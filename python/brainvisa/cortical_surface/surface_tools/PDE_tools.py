# -*- coding: utf-8 -*-
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license version 2 and that you accept its terms.


from __future__ import print_function

from __future__ import absolute_import
import numpy as np
from scipy import sparse
import scipy.stats.stats as sss
from scipy.sparse.linalg import lgmres
from brainvisa.cortical_surface.surface_tools import basic_tools as basicTls
from soma import aims
from six.moves import range
########################
## error tolerance for lgmres solver #
solver_tolerance = 1e-6                      #
########################

####################################################################
#
# distance between the two vertices corresponding to the min and max of the 2d laplacien eigen vector
#
####################################################################
def meshLaplacianEigenVectors(mesh, nbVectors=1):
    from scipy.sparse.linalg import eigsh
    L, LB = computeMeshLaplacian(mesh, lap_type='fem')
    w,v = eigsh(L.tocsr(), nbVectors+1, M=LB.tocsr(), sigma = 0)
    return v[:,1:]

####################################################################
#
# distance between the two vertices corresponding to the min and max of the 2d laplacien eigen vector
#
####################################################################
def meshFiedlerLength(mesh, dist_type='geodesic', fiedler=None):
    if fiedler is None:
        fiedler = meshLaplacianEigenVectors(mesh, 1)
    imin = fiedler.argmin()
    imax = fiedler.argmax()
    vert = np.array(mesh.vertex())

    if dist_type == 'geodesic':
        print('Computing GEODESIC distance between the max and min')
        g = aims.GeodesicPath(mesh, 3, 0)
        dist = g.shortestPath_1_1_len(int(imin), int(imax))

    else:
        print('Computing EUCLIDIAN distance between the max and min')
        min_max = vert[imin, :]-vert[imax, :]
        dist = np.sqrt(np.sum(min_max * min_max, 0))
    return(dist, fiedler)

####################################################################
#
# smoothing the mesh by solving the heat equation using fem Laplacian
#
####################################################################
def laplacianMeshSmoothing(mesh, Niter, dt):
    print('    Smoothing mesh')
    L, B = computeMeshLaplacian(mesh, lap_type='fem')
    avert = np.array(mesh.vertex())
    Mvert = laplacianSmoothing(avert, L, B, Niter, dt)
    vv = aims.vector_POINT3DF()
    for i in range(avert.shape[0]):
        vv.append([Mvert[i, 0], Mvert[i, 1], Mvert[i, 2]])
    smooth = aims.AimsTimeSurface_3_VOID()
    smooth.vertex().assign(vv)
    smooth.polygon().assign(mesh.polygon())
    smooth.updateNormals()
    return smooth

####################################################################
#
# smoothing the texture by solving the heat equation using fem Laplacian
#
####################################################################
def laplacianTextureSmoothing(mesh, tex, Niter, dt):
    print('    Smoothing texture')
    L, B = computeMeshLaplacian(mesh, lap_type='fem')
    return laplacianSmoothing(tex, L, B, Niter, dt)

####################################################################
#
# smoothing the sulcal pits keeping the max at  1
# this is done by solving the heat equation under dirichlet condition at the location of pits
#
####################################################################
def laplacianPitsSmoothing(mesh, tex, Niter, dt):
    print('    Smoothing pits')
    L, B = computeMeshLaplacian(mesh, lap_type='fem')
    inds_pits = np.where(tex == 1)[0]
    #dirichlet condition at the location of pits
    L[inds_pits, :] = 0
    B[inds_pits, :] = 0
    B[inds_pits, inds_pits] = 1
    s_tex = laplacianSmoothing(tex, L, B, Niter, dt)
    # threshold negative values due to the inaccuracy of  matrix inversion to 0
    s_tex[s_tex < 0] = 0
    return s_tex

####################################################################
#
# sub-function for smoothing using fem Laplacian
#
####################################################################
def laplacianSmoothing(Mtex, L, B, Niter, dt):
    mod = 1
    if Niter > 10:
        mod = 10
    if Niter > 100:
        mod = 100
    if Niter > 1000:
        mod = 1000
    #print(tex.shape[0])
    #print(tex.ndim)
    #if tex.ndim < 2:
    #   Mtex = tex.reshape(tex.shape[0],1)
    #else:
     #   Mtex = tex
    # using Implicit scheme
    # B(X^(n+1)-X^n)/dt+L(X^(n+1))=0
    M = B+dt*L
    for i in range(Niter):
        Mtex = B * Mtex
        if Mtex.ndim>1:
            for d in range(Mtex.shape[1]):
                Mtex[:,d], infos = lgmres(M.tocsr(), Mtex[:,d], tol=solver_tolerance)
        else:
            Mtex, infos = lgmres(M.tocsr(), Mtex, tol=solver_tolerance)
        if (i % mod == 0):
            print(i)

    # using Explicit scheme, convergence guaranteed only for dt<1 and not faster than implicit when using fem Laplacian
    # B(X^(n+1)-X^n)/dt+L(X^n)=0
    # M = B-dt*L
    # for i in range(Niter):
    #     Mtex = M * Mtex
    #     Mtex, infos = lgmres(B.tocsr(), Mtex, tol=solver_tolerance)
    #     if (i % mod == 0):
    #         print(i)
    print('    OK')
    return Mtex

####################################################################
# compute_mesh_weight - compute a weight matrix
#
#   W = compute_mesh_weight(vertex,face,type,options);
#
#   W is sparse weight matrix and W(i,j) = 0 is vertex i and vertex j are not
#   connected in the mesh.
#
#   type is either 
#       'combinatorial': W(i,j) = 1 is vertex i is conntected to vertex j.
#       'distance': W(i,j) = 1/d_ij^2 where d_ij is distance between vertex
#           i and j.
#       'conformal': W(i,j) = cot(alpha_ij)+cot(beta_ij) where alpha_ij and
#           beta_ij are the adjacent angle to edge (i,j)
#
#   If options.normalize = 1, the the rows of W are normalize to sum to 1.
#
####################################################################
def computeMeshWeights(mesh, weight_type='conformal', cot_threshold=None, z_threshold=None):
#    cot_threshold=0.00001
 #   print('angle threshold')
    print('    Computing mesh weights')
    vert = np.array(mesh.vertex())
    poly = np.array(mesh.polygon())

    Nbv = vert.shape[0]
    Nbp = poly.shape[0]
    W = sparse.lil_matrix((Nbv, Nbv))
    femB = sparse.lil_matrix((Nbv, Nbv))
    if weight_type == 'conformal' or weight_type == 'fem' :
        threshold = 0.0001 #np.spacing(1)??
        threshold_needed = 0
        for i in range(3):
            i1 = np.mod(i, 3)
            i2 = np.mod(i + 1, 3)
            i3 = np.mod(i + 2, 3)
            pp = vert[poly[:, i2], :] - vert[poly[:, i1], :]
            qq = vert[poly[:, i3], :] - vert[poly[:, i1], :]
            cr  = np.cross(pp,qq)
            area = np.sqrt(np.sum(np.power(cr,2),1))/2
#             nopp = np.apply_along_axis(np.linalg.norm, 1, pp)
#             noqq = np.apply_along_axis(np.linalg.norm, 1, qq)
            noqq = np.sqrt(np.sum(qq * qq, 1))
            nopp = np.sqrt(np.sum(pp * pp, 1))
            thersh_nopp = np.where(nopp<threshold)[0]
            thersh_noqq = np.where(noqq<threshold)[0]
            if len(thersh_nopp) > 0:
                nopp[thersh_nopp] = threshold
                threshold_needed += len(thersh_nopp)
            if len(thersh_noqq) > 0:
                noqq[thersh_noqq] = threshold
                threshold_needed += len(thersh_noqq)
    #        print(np.min(noqq))
            pp = pp / np.vstack((nopp, np.vstack((nopp, nopp)))).transpose()
            qq = qq / np.vstack((noqq, np.vstack((noqq, noqq)))).transpose()
            ang = np.arccos(np.sum(pp * qq, 1))
            ############### preventing infs in weights
            inds_zeros = np.where(ang == 0)[0]
            ang[inds_zeros] = threshold
            threshold_needed_angle = len(inds_zeros)
            ################################
            cot = 1 / np.tan(ang)
            if cot_threshold is not None:
                thresh_inds = cot<0
                cot[thresh_inds]=cot_threshold
                threshold_needed_angle += np.count_nonzero(thresh_inds)
            W = W + sparse.coo_matrix((cot,(poly[:, i2],poly[:, i3])),shape=(Nbv, Nbv))
            W = W + sparse.coo_matrix((cot,(poly[:, i3],poly[:, i2])),shape=(Nbv, Nbv))
            femB = femB + sparse.coo_matrix((area/12,(poly[:, i2],poly[:, i3])),shape=(Nbv, Nbv))
            femB = femB + sparse.coo_matrix((area/12,(poly[:, i3],poly[:, i2])),shape=(Nbv, Nbv))

        # if weight_type == 'fem' :
        #     W.data = W.data/2

        nnz = W.nnz
        if z_threshold is not None:
            z_weights = sss.zscore(W.data)
            inds_out = np.where(np.abs(z_weights) > z_thresh)[0]
            W.data[inds_out] = np.mean(W.data)
            print('    -Zscore threshold needed for ',len(inds_out),' values = ', 100*len(inds_out)/nnz,' %')
        #inds_out_inf = np.where(z_weights < -z_thresh)[0]
        #inds_out_sup = np.where(z_weights > z_thresh)[0]
        #val_inf = np.max(W.data[inds_out_inf])
        #W.data[inds_out_inf] = val_inf
        #val_sup = np.min(W.data[inds_out_sup])
        #W.data[inds_out_sup] = val_sup
        #print('    -Zscore threshold needed for ',len(inds_out_inf)+len(inds_out_sup),' values-')
        print('    -edge length threshold needed for ',threshold_needed,' values = ', 100*threshold_needed/nnz,' %')
        if cot_threshold is not None:
            print('    -cot threshold needed for ',threshold_needed_angle,' values = ', 100*threshold_needed_angle/nnz,' %')

    li = np.hstack(W.data)
    nb_Nan = len(np.where(np.isnan(li))[0])
    nb_neg = len(np.where(li<0)[0])
    print('    -number of Nan in weights: ',nb_Nan ,' = ', 100*nb_Nan/nnz,' %')
    print('    -number of Negative values in weights: ', nb_neg,' = ',100*nb_neg/nnz,' %')

    return W.tocsr(),femB.tocsr()

####################################################################
#
# compute laplacian of a mesh
#
####################################################################
def computeMeshLaplacian(mesh, weights=None, femB=None, lap_type='conformal'):
    print('    Computing Laplacian')
    if weights is None:
        (weights, femB) = computeMeshWeights(mesh, weight_type=lap_type)

    if lap_type == 'fem':
        weights.data = weights.data/2

    N = weights.shape[0]
    sB = femB.sum(axis=0)
    diaB = sparse.dia_matrix((sB, 0), shape=(N, N))
    B = sparse.lil_matrix(diaB + femB)
    s = weights.sum(axis=0)
    dia = sparse.dia_matrix((s, 0), shape=(N, N))
    L = sparse.lil_matrix(dia - weights)

    li = np.hstack(L.data)
    print('    -nb Nan in L : ', len(np.where(np.isnan(li))[0]))
    print('    -nb Inf in L : ', len(np.where(np.isinf(li))[0]))

    return L,B

####################################################################
#
# compute vertex voronoi of a mesh as described in
#    Meyer, M., Desbrun, M., Schröder, P., & Barr, A. (2002). 
#    Discrete differential-geometry operators for triangulated 2-manifolds. 
#    Visualization and Mathematics, 1–26.
#
####################################################################
def vertexVoronoi(mesh, angs=None):
    print('    Computing Vertex Voronoi')
    vert = np.array(mesh.vertex())
    poly = np.array(mesh.polygon())
    if angs is None:
        angs = basicTls.meshPolygonAngles(vert, poly)
    Nbv = vert.shape[0]
    Nbp = poly.shape[0]
    areas = basicTls.meshPolygonArea(vert,poly)
    obt_angs = angs>np.pi/2
    obt_poly = obt_angs[:,0]|obt_angs[:,1]|obt_angs[:,2]
    print('    -percent polygon with obtuse angle ',100.0*len(np.where(obt_poly)[0])/Nbp)
    cot = 1 / np.tan(angs)
    vert_voronoi = np.zeros(Nbv)
    for ind_p,p in enumerate(poly):
        if obt_poly[ind_p]:
            obt_verts = p[obt_angs[ind_p,:]]
            vert_voronoi[obt_verts] = vert_voronoi[obt_verts] + areas[ind_p]/2.0
            non_obt_verts = p[[not x for x in obt_angs[ind_p,:]]]
            vert_voronoi[non_obt_verts] = vert_voronoi[non_obt_verts] + areas[ind_p]/4.0
        else:
            d0 = np.sum( np.power(vert[p[1], :] - vert[p[2], :],2))
            d1 = np.sum( np.power(vert[p[2], :] - vert[p[0], :],2))
            d2 = np.sum( np.power(vert[p[0], :] - vert[p[1], :],2))
            vert_voronoi[p[0]] = vert_voronoi[p[0]] + (d1*cot[ind_p,1] + d2*cot[ind_p,2])/8.0
            vert_voronoi[p[1]] = vert_voronoi[p[1]] + (d2*cot[ind_p,2] + d0*cot[ind_p,0])/8.0
            vert_voronoi[p[2]] = vert_voronoi[p[2]] + (d0*cot[ind_p,0] + d1*cot[ind_p,1])/8.0
# 
#         
#     print(cot.shape)
#     W = sparse.lil_matrix((Nbv, Nbv))
#     W1 = sparse.lil_matrix((Nbv, Nbv))
#     for i in range(3):
#         i1 = np.mod(i, 3)
#         i2 = np.mod(i + 1, 3)
#         i3 = np.mod(i + 2, 3)
#         print((i2,i3))
#         d = np.sum( np.power(vert[poly[:, i2], :] - vert[poly[:, i3], :],2), 1)
#         W1 = W1 + sparse.coo_matrix((d,(poly[:, i2],poly[:, i3])),shape=(Nbv, Nbv))
#         W = W + sparse.coo_matrix((cot[:,i],(poly[:, i2],poly[:, i3])),shape=(Nbv, Nbv))
#         W = W + sparse.coo_matrix((cot[:,i],(poly[:, i3],poly[:, i2])),shape=(Nbv, Nbv))
# 
#     Af = W1.multiply(W).tolil()
#     tAf = sparse.triu(Af,0)
#     vert_voronoi2 = tAf.sum(0)
#     vert_voronoi2 = vert_voronoi2/4
#     vert_voronoi2 = np.array(vert_voronoi2.transpose()).squeeze()
#     diff = vert_voronoi2-vert_voronoi
#     print(len(np.where(diff)[0]))
#     print(vert_voronoi.sum())
#     print(vert_voronoi2.sum())
    return vert_voronoi

####################################################################
#
# compute the depth potential function of a mesh as desribed in
#    Boucher, M., Whitesides, S., & Evans, A. (2009). 
#    Depth potential function for folding pattern representation, registration and analysis. 
#    Medical Image Analysis, 13(2), 203–14. 
#    doi:10.1016/j.media.2008.09.001
#
####################################################################
def depthPotentialFunction(mesh, curvature, alphas):
    L, LB = computeMeshLaplacian(mesh, lap_type='fem')
    Nbv = len(np.array(mesh.vertex()))
    B = -LB * (curvature-( np.sum(curvature*LB.diagonal()) / np.sum(LB.diagonal()) ))

    dpf=[]
    for ind, alpha in enumerate(alphas):
        M = alpha*LB+L/2
        dpf_t, info = lgmres(M.tocsr(), B, tol=solver_tolerance)
        dpf.append(dpf_t)

    ############################
    # old, slower and less accurate implementation using conformal laplacian instead of fem
    ############################
    # vert_voronoi = vertexVoronoi(mesh)
    # L, LB = computeMeshLaplacian(mesh, lap_type='conformal')
    # B = -2 * vert_voronoi * (curvature-( np.sum(curvature*vert_voronoi) / vert_voronoi.sum() ))
    # B=B.squeeze()
    # for ind, alpha in enumerate(alphas):
    #     A = sparse.dia_matrix((alpha*vert_voronoi, 0), shape=(Nbv, Nbv))
    #     M = A+L
    #     dpf_t, info = lgmres(M.tocsr(), B, tol=solver_tolerance)
    #     dpf.append(dpf_t)
    return dpf

####################################################################
#
# compute the mean curvature as detailed in
#    Meyer, M., Desbrun, M., Schröder, P., & Barr, A. (2002). 
#    Discrete differential-geometry operators for triangulated 2-manifolds. 
#    Visualization and Mathematics, 1–26.

#
####################################################################
def meanCurvature(mesh):
    vert = np.array(mesh.vertex())
    poly = np.array(mesh.polygon())
    angs = basicTls.meshPolygonAngles(vert, poly)
    Nbv = len(vert)

    vert_voronoi = vertexVoronoi(mesh,angs)
    mean_curv = np.pi*np.ones(Nbv)
    for ind_p,p in enumerate(poly):
        mean_curv[p[0]] = mean_curv[p[0]] - angs[ind_p,0]
        mean_curv[p[1]] = mean_curv[p[1]] - angs[ind_p,1]
        mean_curv[p[2]] = mean_curv[p[2]] - angs[ind_p,2]
    mean_curv_out = 3*mean_curv / vert_voronoi
#     cot = 1 / np.tan(angs)
#     mean_curv = np.zeros((Nbv,3))
#     for ind_p,p in enumerate(poly):
#         d0 = vert[p[1], :] - vert[p[2], :]
#         d1 = vert[p[2], :] - vert[p[0], :]
#         d2 = vert[p[0], :] - vert[p[1], :]
#         mean_curv[p[0],:] = mean_curv[p[0],:] + d1*cot[ind_p,1] + d2*cot[ind_p,2]
#         mean_curv[p[1],:] = mean_curv[p[1],:] + d2*cot[ind_p,2] + d0*cot[ind_p,0]
#         mean_curv[p[2],:] = mean_curv[p[2],:] + d0*cot[ind_p,0] + d1*cot[ind_p,1]
#     mean_curv_out = np.sqrt(np.sum( np.power(mean_curv,2),1)) / (2*vert_voronoi)
    return mean_curv_out
