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
def computeMeshWeights(mesh, weight_type=None):
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
        for i in range(3):
            i1 = np.mod(i, 3)
            i2 = np.mod(i + 1, 3)
            i3 = np.mod(i + 2, 3)
            print '    ', 3 - i1
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
           
        if threshold_needed > 0:
            print '    -weight threshold needed for ',threshold_needed,' values-'
        print '    OK'
    li = np.hstack(W.data)
    print 'percent of Nan in weights: ', 100*len(np.where(np.isnan(li))[0])/len(li)
    print 'percent of Negative values in weights: ', 100*len(np.where(li<0)[0])/len(li)

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
    print 'nb Nan in L : ', len(np.where(np.isnan(li))[0])
    print 'nb Inf in L : ', len(np.where(np.isinf(li))[0])   
    print '    OK'

    return L

####################################################################
#
# compute vertex voronoi of a mesh as described in
#    Meyer, M., Desbrun, M., Schröder, P., & Barr, A. (2002). 
#    Discrete differential-geometry operators for triangulated 2-manifolds. 
#    Visualization and Mathematics, 1–26.
#
####################################################################
def vertexVoronoi(mesh, weights=None):
    if weights is None:
        weights = computeMeshWeights(mesh)

#     I = np.array([0,0,1,3,1,0,0])
#     J = np.array([0,2,1,3,1,0,0])
#     V = np.array([1,1,1,1,1,1,1])
#     B = sparse.coo_matrix((V,(I,J)),shape=(4,4)).tocsr()
#     print 'B',B
#  
#     I = np.array([0,0,1,3,1,0,0])
#     J = np.array([0,2,1,3,1,0,0])
#     V = np.array([1,1,1,1,1,1,1])
#     B1 = sparse.coo_matrix((V,(I,J)),shape=(4,4)).tocsr()
#     print 'B1*B',B1*B
#     print 'B*B1',B*B1

    print '    Computing Vertex Voronoi'
    vert = np.array(mesh.vertex())
    poly = np.array(mesh.polygon())

    Nbv = vert.shape[0]

    adja = basicTls.meshAdjacencyMatrix(mesh)
    (i,j) = adja.nonzero()
    d = np.sum( np.power(vert[i,:] - vert[j,:],2), 1)
    W1 = sparse.coo_matrix((d, (i, j)),shape=(Nbv, Nbv))  #    W1(i,j) = d_ij^2
    Ae = W1.multiply(weights).tolil()
    print type(Ae)
    (I,J,V) = sparse.find(Ae)
    inds_negs=np.where(V<0)[0]
    V[inds_negs]=0
    Ae2 = sparse.coo_matrix((V,(I,J)),shape=(Nbv, Nbv))
    (I,J,V) = sparse.find(Ae2)
    print len(np.where(V<0)[0])


#    print 'neigh[0].list()',neigh[0].list()
    vert_voronoi = Ae2.sum(0)
    vert_voronoi = vert_voronoi/4
#     for i,nei in enumerate(neigh):
#         ne_i = np.array(neigh[i].list())
#         vert_voronoi[i]=1/2*np.sum(Ae,0)

    return np.array(vert_voronoi.transpose()).squeeze()

####################################################################
#
# compute the depth potential function of a mesh as desribed in
#    Boucher, M., Whitesides, S., & Evans, A. (2009). 
#    Depth potential function for folding pattern representation, registration and analysis. 
#    Medical Image Analysis, 13(2), 203–14. 
#    doi:10.1016/j.media.2008.09.001
#
####################################################################
def depthPotentialFunction(mesh, curvature, alpha):
    weights = computeMeshWeights(mesh)
    vert_voronoi = vertexVoronoi(mesh, weights)
    L = computeMeshLaplacian(mesh, weights)
    print vert_voronoi.shape
    print type(vert_voronoi)
    print vert_voronoi.sum()
#     v=np.array(vert_voronoi).squeeze()
#     print 'min',np.min(v)
#     print 'max',np.max(v)
#     print 'sum',np.sum(v)
#     print len(np.where(v<0)[0])
    Nbv = len(np.array(mesh.vertex()))
    solver_tolerance = 1e-6
    A = sparse.dia_matrix((alpha*vert_voronoi, 0), shape=(Nbv, Nbv))
    M = A+L
    M = M.tocsr()
#    B = sparse.dia_matrix((2 * np.array(vert_voronoi).squeeze() * (k-( np.sum(k*np.array(vert_voronoi).squeeze()) / vert_voronoi.sum() )), 0), shape=(Nbv, Nbv))
#    B = B.tocsr()
    B = 2 * vert_voronoi * (curvature-( np.sum(curvature*vert_voronoi) / vert_voronoi.sum() ))
    print B.shape
    B=B.squeeze()
    print B.shape
#    B = sparse.csr_matrix(B)
    print M.shape
    print B.shape
    dpf, info = lgmres(M.tocsr(), B, tol=solver_tolerance)
    
    return dpf

