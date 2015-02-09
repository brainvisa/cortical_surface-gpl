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


import numpy as np
from scipy import sparse
from scipy.sparse.linalg import lgmres
from brainvisa.cortical_surface.surface_tools import basic_tools as basicTls

####################################################################
#
# EUCLIDIAN distance between the two vertices corresponding to the min and max of the 2d laplacien eigen vector
#
####################################################################
def meshFiedlerLength(mesh):
    from scipy.sparse.linalg import eigsh
    L=computeMeshLaplacian(mesh)
    Lap=0.5*(L+L.transpose())
    print 'Computing fiedler vector'
    w,v=eigsh(Lap, 2, which='LM', sigma = 0)
 
    fiedler=v[:,1]
    print 'Computing EUCLIDIAN distance between the max and min'
    imin=fiedler.argmin()
    imax=fiedler.argmax()
    vert = np.array(mesh.vertex())
    min_max = vert[imin,:]-vert[imax,:]
    dist = np.sqrt(np.sum(min_max * min_max, 0))
    return(dist,fiedler)

####################################################################
#
# laplacian pits smoothing
#
####################################################################
def meshPitsSmoothing(mesh, tex,Niter, dt):
    print '    Smoothing mesh'
    mod = 1
    if Niter > 10:
        mod = 10
    if Niter > 100:
        mod = 100
    if Niter > 1000:
        mod = 1000
    #vert = np.array(mesh.vertex())
    print 'using conformal weights with angular threshold at 0.0001'
    weights = computeMeshWeights(mesh,'conformal',0.0001)
    N = weights.shape[0]
#    s = weights.sum(axis = 1)
    s = weights.sum(axis=0)
#    dia = sparse.lil_matrix((N,N))
#    dia.setdiag(s)
    dia = sparse.dia_matrix((1 / s, 0), shape=(N, N))
#    b = ones(N)
#    W = sparse.linalg.spsolve(dia.tocsr(),b)*weights
    W = dia * weights
    I = sparse.lil_matrix((N, N))
    I.setdiag(np.ones(N))
    tL = I.tocsr() - W
    LI = I.tocsr() - (dt * tL)
    #    LI = dt*L
    #    LI = LI.tocsr()
    #    LI = I.tocsr() + LI
    #    print '    LI.shape = ', LI.shape
    vert=np.array(tex[0]).reshape(N,1)
    Mtex = vert#sparse.lil_matrix(vert).tocsr()
    inds_pits=np.where(vert==1)[0]
    #print inds_pits[:10]
    #print '    Mtex.shape = ', Mtex.shape
    #print inds_pits.shape
    #print Mtex[inds_pits]
    o = np.ones((inds_pits.shape[0],1))

    #LI2 = I.tocsr() + (dt * tL)
    #Mtex = lgmres(LI2.tocsr(), Mtex, tol=solver_tolerance)
    print 'iterative filtering the texture...'
    for i in range(Niter):
        Mtex= LI * Mtex
        Mtex[inds_pits]=o
        if (i % mod == 0):
            print i
    print '    OK'

    return(Mtex)

####################################################################
#
# laplacian smoothing
#
####################################################################
def meshSmoothing(mesh, Niter, dt):
    print '    Smoothing mesh'
    mod = 1
    if Niter > 10:
        mod = 10
    if Niter > 100:
        mod = 100
    if Niter > 1000:
        mod = 1000
    vert = np.array(mesh.vertex())
    Mvert = sparse.lil_matrix(vert).tocsr()
    weights = computeMeshWeights(mesh)
    N = weights.shape[0]
#    s = weights.sum(axis = 1)
    s = weights.sum(axis=0)
#    dia = sparse.lil_matrix((N,N))
#    dia.setdiag(s)
    dia = sparse.dia_matrix((1 / s, 0), shape=(N, N))
#    b = ones(N)
#    W = sparse.linalg.spsolve(dia.tocsr(),b)*weights
    W = dia * weights
    I = sparse.lil_matrix((N, N))
    I.setdiag(np.ones(N))
    tL = I.tocsr() - W
    LI = I.tocsr() - (dt * tL)
    #    LI = dt*L
    #    LI = LI.tocsr()
    #    LI = I.tocsr() + LI
    #    print '    LI.shape = ', LI.shape
    print '    Mvert.shape = ', Mvert.shape
    for i in range(Niter):
        Mvert = LI * Mvert
        if (i % mod == 0):
            print i
    print '    OK'
    vv = aims.vector_POINT3DF()
    for i in range(N):
        vv.append([Mvert[i, 0], Mvert[i, 1], Mvert[i, 2]])

    smooth = aims.AimsTimeSurface_3()
    smooth.vertex().assign(vv)
    smooth.polygon().assign(mesh.polygon())
    smooth.updateNormals()
    return(smooth)

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
def computeMeshWeights(mesh, weight_type=None, angle_threshold=None):
#    angle_threshold=0.00001
 #   print 'angle threshold'
    if weight_type is None:
        weight_type = 'conformal'
    print '    Computing mesh weights'
    vert = np.array(mesh.vertex())
    poly = np.array(mesh.polygon())

    Nbv = vert.shape[0]
    Nbp = poly.shape[0]
    W = sparse.lil_matrix((Nbv, Nbv))
#    W1 = sparse.lil_matrix((Nbv, Nbv))
    if weight_type == 'conformal':
        threshold = 0.0001 #np.spacing(1)??
        threshold_needed = 0
        threshold_needed_angle = 0
        for i in range(3):
            i1 = np.mod(i, 3)
            i2 = np.mod(i + 1, 3)
            i3 = np.mod(i + 2, 3)
            pp = vert[poly[:, i2], :] - vert[poly[:, i1], :]
            qq = vert[poly[:, i3], :] - vert[poly[:, i1], :]
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
    #        print np.min(noqq)
            pp = pp / np.vstack((nopp, np.vstack((nopp, nopp)))).transpose()
            qq = qq / np.vstack((noqq, np.vstack((noqq, noqq)))).transpose()
            ang = np.arccos(np.sum(pp * qq, 1))
            cot = 1 / np.tan(ang)
            if angle_threshold is not None:
                thresh_inds = cot<0
                cot[thresh_inds]=angle_threshold
                threshold_needed_angle += np.count_nonzero(thresh_inds)
            W = W + sparse.coo_matrix((cot,(poly[:, i2],poly[:, i3])),shape=(Nbv, Nbv))
            W = W + sparse.coo_matrix((cot,(poly[:, i3],poly[:, i2])),shape=(Nbv, Nbv))

#             for j in range(Nbp):
#     #            ind1 = poly[j, i1]
#                 ind2 = poly[j, i2]
#                 ind3 = poly[j, i3]
#                 W[ind2, ind3] = W[ind2, ind3] + 1 / np.tan(ang[j])
#                 W[ind3, ind2] = W[ind3, ind2] + 1 / np.tan(ang[j])

#            W[poly[:, i2],poly[:, i3]] = W[poly[:, i2],poly[:, i3]] + cot
#            W[poly[:, i3],poly[:, i2]] = W[poly[:, i3],poly[:, i2]] + cot
           
        print '    -angle threshold needed for ',threshold_needed_angle,' values-'
        print '    -weight threshold needed for ',threshold_needed,' values-'
    li = np.hstack(W.data)
    print '    -percent of Nan in weights: ', 100*len(np.where(np.isnan(li))[0])/len(li)
    print '    -percent of Negative values in weights: ', 100*len(np.where(li<0)[0])/len(li)

    return W


####################################################################
#
# compute laplacian of a mesh
#
####################################################################
def computeMeshLaplacian(mesh, weights=None):
    print '    Computing Laplacian'
    if weights is None:
        weights = computeMeshWeights(mesh)
    N = weights.shape[0]
#    s = weights.sum(axis = 1)
    s = weights.sum(axis=0)
#    dia = sparse.lil_matrix((N,N))
#    dia.setdiag(s)
    dia = sparse.dia_matrix((s, 0), shape=(N, N))

#    print dia - weights
    L = sparse.lil_matrix(dia - weights)
    li = np.hstack(L.data)
    print '    -nb Nan in L : ', len(np.where(np.isnan(li))[0])
    print '    -nb Inf in L : ', len(np.where(np.isinf(li))[0])   

    return L

####################################################################
#
# compute vertex voronoi of a mesh as described in
#    Meyer, M., Desbrun, M., Schröder, P., & Barr, A. (2002). 
#    Discrete differential-geometry operators for triangulated 2-manifolds. 
#    Visualization and Mathematics, 1–26.
#
####################################################################
def vertexVoronoi(mesh, angs=None):
    print '    Computing Vertex Voronoi'
    vert = np.array(mesh.vertex())
    poly = np.array(mesh.polygon())
    if angs is None:
        angs = basicTls.meshPolygonAngles(vert, poly)
    Nbv = vert.shape[0]
    Nbp = poly.shape[0]
    areas = basicTls.meshPolygonArea(vert,poly)
    obt_angs = angs>np.pi/2
    obt_poly = obt_angs[:,0]|obt_angs[:,1]|obt_angs[:,2]
    print '    -percent polygon with obtuse angle ',100.0*len(np.where(obt_poly)[0])/Nbp
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
#     print cot.shape
#     W = sparse.lil_matrix((Nbv, Nbv))
#     W1 = sparse.lil_matrix((Nbv, Nbv))
#     for i in range(3):
#         i1 = np.mod(i, 3)
#         i2 = np.mod(i + 1, 3)
#         i3 = np.mod(i + 2, 3)
#         print (i2,i3)
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
#     print len(np.where(diff)[0])
#     print vert_voronoi.sum()
#     print vert_voronoi2.sum()
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
    vert_voronoi = vertexVoronoi(mesh)
    L = computeMeshLaplacian(mesh)
#     v=np.array(vert_voronoi).squeeze()
#     print 'min',np.min(v)
#     print 'max',np.max(v)
#     print 'sum',np.sum(v)
#     print len(np.where(v<0)[0])
    Nbv = len(np.array(mesh.vertex()))
    solver_tolerance = 1e-6
#    B = sparse.dia_matrix((2 * np.array(vert_voronoi).squeeze() * (k-( np.sum(k*np.array(vert_voronoi).squeeze()) / vert_voronoi.sum() )), 0), shape=(Nbv, Nbv))
#    B = B.tocsr()
    B = -2 * vert_voronoi * (curvature-( np.sum(curvature*vert_voronoi) / vert_voronoi.sum() ))
    B=B.squeeze()
#    B = sparse.csr_matrix(B)
    dpf=[]
    for ind, alpha in enumerate(alphas): 
        A = sparse.dia_matrix((alpha*vert_voronoi, 0), shape=(Nbv, Nbv))
        M = A+L
        M = M.tocsr()
        dpf_t, info = lgmres(M.tocsr(), B, tol=solver_tolerance)
        dpf.append(dpf_t)
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
