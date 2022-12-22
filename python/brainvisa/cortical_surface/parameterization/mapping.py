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
from brainvisa.cortical_surface.parameterization import sulcalLine as sln
from brainvisa.cortical_surface.parameterization import sulcalLinesSet as slSet
from brainvisa.cortical_surface.parameterization import model as md
from brainvisa.cortical_surface.surface_tools import basic_tools as basicTls
from brainvisa.cortical_surface.surface_tools import PDE_tools as pdeTls
from scipy import sparse
import scipy
import six
from six.moves import range


#########################################################
# tolerance parameter for the linear system          ####
# solver in scipy                                    ####
# smaller value makes solver slower but more precise ####
solver_tolerance = 1e-6                              ####
#########################################################
ver = [1]#[ int(x) for x in scipy.__version__.split( '.' ) ]
if ver < [ 0, 9 ]:
#    print('HIP-HOP :: scipy is too old, using gmres for solving linear systems (will be slower), with tolerance = ',solver_tolerance)
    from scipy.sparse.linalg import gmres
    solver = 'gmres'
else:
#    print('HIP-HOP :: using lgmres for solving linear systems, with tolerance = ',solver_tolerance)
    from scipy.sparse.linalg import lgmres
    solver = 'lgmres'  
import numpy as np
from soma import aims, aimsalgo
import time


#########################################################
# tolerance parameter for the linear system          ####
# solver in scipy                                    ####
# smaller value makes solver slower but more precise ####
solver_tolerance = 1e-6                              ####
#########################################################



#class Mapping(object):
#    '''
#    classdocs
#    '''
#
#    def __init__(self, mesh_orig=None):
#        self.mesh_orig = mesh_orig
#        self.mesh_mapped = []
#        self.induced_coords = []
#        self.regularization = 'conformal'
#        self.reference_space = 'sphere'
####################################################################
#
# compute comformal mapping of a mesh to a sphere
#
####################################################################
def sphereConformalMapping(mesh):

    print('    Spherical mapping')
    L, LB = computeMeshLaplacian(mesh, lap_type='conformal')
    #print('Laplacian : ', L)

    Nv = np.array(mesh.vertex()).shape[0]

    Nor = np.array(mesh.normal())

    Rx = sparse.lil_matrix(Nor[:, 0]).tocsr()
    Ry = sparse.lil_matrix(Nor[:, 1]).tocsr()
    Rz = sparse.lil_matrix(Nor[:, 2]).tocsr()

    print('    Solving linear system')
    L = L.tocsr()
    x = spsolve(L, Rx)
    y = spsolve(L, Ry)
    z = spsolve(L, Rz)

    print('    OK')
    vv = aims.vector_POINT3DF()
    for i in range(Nv):
        vv.append([x[i], y[i], z[i]])

    disk = aims.AimsTimeSurface_3_VOID()
    disk.vertex().assign(vv)
    disk.polygon().assign(mesh.polygon())
    disk.updateNormals()
    return(disk)


####################################################################
#
# compute comformal mapping of a mesh to a disk
#
####################################################################
def diskConformalMapping(mesh, boundary=None, boundary_coords=None):
    if boundary is None:
        boundary = basicTls.meshBoundary(mesh)
        #print('Boundary: ', boundary)
    boundary = np.array(boundary)
    if boundary_coords is None:
        p = boundary.size
        t = np.arange(0, 2 * np.math.pi, (2 * np.math.pi / p))
        boundary_coords = np.array([np.cos(t), np.sin(t)])
    L, LB = pdeTls.computeMeshLaplacian(mesh, lap_type='conformal')
    #print('Laplacian : ', L)
    #vert = np.array(mesh.vertex())
    Nv = len(mesh.vertex())  # np.array(mesh.vertex()).shape[0]

    print('Boundary Size:', boundary.shape)
    print('Laplacian Size:', L.shape)
    for i in boundary:
        L[i, :] = 0
        L[i, i] = 1
    L = L.tocsr()

    Rx = np.zeros(Nv)
    Ry = np.zeros(Nv)
    Rx[boundary] = boundary_coords[0, :]
    Ry[boundary] = boundary_coords[1, :]

#    Rx = sparse.lil_matrix(Rx).tocsr()
#    Ry = sparse.lil_matrix(Ry).tocsr()
    if solver == 'lgmres':
        x, info = lgmres(L, Rx, tol=solver_tolerance)
        y, info = lgmres(L, Ry, tol=solver_tolerance)
    else:
        x, info = gmres(L, Rx, tol=solver_tolerance)
        y, info = gmres(L, Ry, tol=solver_tolerance)

#    y = spsolve(L, Ry)
    z = np.zeros(Nv)

    vv = aims.vector_POINT3DF()
    for i in range(Nv):
        vv.append([x[i], y[i], z[i]])

    disk = aims.AimsTimeSurface_3_VOID()
    disk.vertex().assign(vv)
    disk.polygon().assign(mesh.polygon())
    disk.updateNormals()
    return(disk)


####################################################################
#
# compute comformal mapping of a mesh to a rectangle
#
# the mesh is already topologically 'rectangular'
# boundaries is an array wit 4 lines, each of them being an ordered
# list of indexes.
# ratio is the length/width ratio (integer)
#
####################################################################
def rectConformalMapping(mesh, boundary, length, width, fixed_boundary=0):
    "boundaries (see path2Boundary for details:"
    "boundary[0] == insula_boundary"
    "boundary[1] == neocortex_poles_path always from insula to cingular pole"
    "boundary[2] == cingular_boundary"
    "boundary[3] == new vertices always from insula to cingular pole"
    print('mapping to the rectangle  ', length, ' x ', width, 'with fixed_boundary = ', fixed_boundary)
    #print('Laplacian : ', L)
    Nbv = np.array(mesh.vertex()).shape[0]

    Rx = np.zeros(Nbv)
    Ry = np.zeros(Nbv)

    Lx, LB = pdeTls.computeMeshLaplacian(mesh, lap_type='conformal')

    if fixed_boundary:
        for i in boundary[3]:
            Lx[i, :] = 0
            Lx[i, i] = 1
        for i in boundary[1]:
            Lx[i, :] = 0
            Lx[i, i] = 1
        for i in boundary[0]:
            Lx[i, :] = 0
            Lx[i, i] = 1
        for i in boundary[2]:
            Lx[i, :] = 0
            Lx[i, i] = 1
        vert = np.array(mesh.vertex())

        right = np.zeros(len(boundary[3]))
        for v in range(len(boundary[3])):
            right[v] = right[v - 1] + np.linalg.norm(vert[boundary[3][v], :] - vert[boundary[3][v - 1], :])
        right = width * right / max(right)
        left = np.zeros(len(boundary[1]))
        for v in range(1, len(boundary[1])):
            left[v] = left[v - 1] + np.linalg.norm(vert[boundary[1][v], :] - vert[boundary[1][v - 1], :])
        left = width * left / max(left)
        bottom = np.zeros(len(boundary[0]))  # insula boundary
        for v in range(1, len(boundary[0])):
            bottom[v] = bottom[v - 1] + np.linalg.norm(vert[boundary[0][v], :] - vert[boundary[0][v - 1], :])
        bottom = length * bottom / max(bottom)
        top = np.zeros(len(boundary[2]))  # cingular boundary
        for v in range(1, len(boundary[2])):
            top[v] = top[v - 1] + np.linalg.norm(vert[boundary[2][v], :] - vert[boundary[2][v - 1], :])
        top = length * top / max(top)
        Rx[boundary[3]] = length * np.ones(len(boundary[3]))
        Ry[boundary[3]] = right
        Rx[boundary[1]] = np.zeros(len(boundary[1]))
        Ry[boundary[1]] = left
        Rx[boundary[0]] = bottom
        Ry[boundary[0]] = np.zeros(len(boundary[0]))
        Rx[boundary[2]] = top
        Ry[boundary[2]] = width * np.ones(len(boundary[1]))

        print('solve the linear system')
        Lx = Lx.tocsr()
        Rx = sparse.csr_matrix(Rx)
        Ry = sparse.csr_matrix(Ry)
#        mtx1 = mtx.astype(np.float32)
#        print('using gmres')
#        result, info = gmres(Lx, Rx, tol=1e-3)
        print('using umfpack')
        x = spsolve(Lx, Rx, use_umfpack=True)
        print('without umfpack')
        x = spsolve(Lx, Rx)
        print('solve for y')
        y = spsolve(Lx, Ry, use_umfpack=True)

    else:
        Ly = Lx.copy()
        for i in boundary[3]:
            Lx[i, :] = 0
            Lx[i, i] = 1
        for i in boundary[1]:
            Lx[i, :] = 0
            Lx[i, i] = 1
        for i in boundary[0]:
            Ly[i, :] = 0
            Ly[i, i] = 1
        for i in boundary[2]:
            Ly[i, :] = 0
            Ly[i, i] = 1

        right = length * np.ones(len(boundary[3]))
        left = np.zeros(len(boundary[1]))
        bottom = np.zeros(len(boundary[0]))  # insula boundary
        top = width * np.ones(len(boundary[2]))  # cingular boundary
        Rx[boundary[3]] = right
        Rx[boundary[1]] = left
        Ry[boundary[0]] = bottom
        Ry[boundary[2]] = top

        print('solve the linear system')
        Lx = Lx.tocsc()
        Ly = Ly.tocsc()
#        print(Rx.shape)
#        Rxt=Rx.transpose()
#        print(Rxt.shape)
#        Rx = sparse.lil_matrix(Rx).tocsr()
#        Ry = sparse.lil_matrix(Ry).tocsr()
        print(Rx.shape)
#        print('using spsolve')
#        t0 = time.perf_counter()
#        x_ex = spsolve(Lx, Rx)
#        print(time.perf_counter() - t0, "seconds process time")
#        print('using splu')
#        t0 = time.perf_counter()
#        print('    precond')
#        P = splu(Lx, drop_tol=1e-1)
#        print(time.perf_counter() - t0, "seconds process time")
#        print('    solve')
#        t0 = time.perf_counter()
#        x0 = P.solve(Rx)
#        print(time.perf_counter() - t0, "seconds process time")
#        print('error:',np.linalg.norm(x_ex-x0))
#        M_x = lambda x: P.solve(x)
#        M = scipy.sparse.linalg.LinearOperator((n * m, n * m), M_x)
#        result = scipy.sparse.linalg.lgmres(matrix, b, tol=1e-4, M=M)[0]
#        print('using bicgstab')
#        t0 = time.perf_counter()
#        x1, info = bicgstab(Lx, Rx)#, tol=solver_tolerance)
##        print(info)
#        print(time.perf_counter() - t0, "seconds process time")
#        print('error:',np.linalg.norm(x_ex-x1))
        t0 = time.perf_counter()
        if solver == 'lgmres':
            x, info = lgmres(Lx, Rx, tol=solver_tolerance)
            y, info = lgmres(Ly, Ry, tol=solver_tolerance)#spsolve(Ly, Ry)
            print('using lgmres')
        else:
            print('using gmres')
            x, info = gmres(Lx, Rx, tol=solver_tolerance)
            y, info = gmres(Ly, Ry, tol=solver_tolerance)
#        print(info)
        print(time.perf_counter() - t0, "seconds process time")
    print('matrix inverted')
    z = np.zeros(Nbv)
    vv = aims.vector_POINT3DF()
    for i in range(Nbv):
        vv.append([x[i], y[i], z[i]])

    rect = aims.AimsTimeSurface_3_VOID()
    rect.vertex().assign(vv)
    rect.polygon().assign(mesh.polygon())
    rect.updateNormals()
    return(rect)


####################################################################
#
# compute constrained comformal mapping of a mesh in a rectangle
#
# the mesh is already topologically 'rectangular'
# boundaries is an array wit 4 lines, each of them being an ordered
# list of indexes.
# ratio is the length/width ratio (integer)
#
####################################################################
def cstrRectConformalMapping(Lx, modele, mesh, boundary, sulcalCstr, cstrBalance):
    "boundaries (see path2Boundary for details):"
    "boundary[0] == insula_boundary"
    "boundary[1] == neocortex_poles_path always from insula to cingular pole"
    "boundary[2] == cingular_boundary"
    "boundary[3] == new vertices always from insula to cingular pole"
    print('cstr mapping in the rectangle  with cstrBalance = ', cstrBalance)
    vert = np.array(mesh.vertex())
    Nbv = vert.shape[0]

    Rx = np.zeros(Nbv)
    Ry = np.zeros(Nbv)

    Ly = Lx.copy()
    for i in boundary[3]:
        Lx[i, :] = 0
        Lx[i, i] = 1
    for i in boundary[1]:
        Lx[i, :] = 0
        Lx[i, i] = 1
    for i in boundary[0]:
        Ly[i, :] = 0
        Ly[i, i] = 1
    for i in boundary[2]:
        Ly[i, :] = 0
        Ly[i, i] = 1

    Rx[boundary[3]] = vert[boundary[3], 0]
    Rx[boundary[1]] = vert[boundary[1], 0]
    Ry[boundary[0]] = vert[boundary[0], 1]  # insula boundary
    Ry[boundary[2]] = vert[boundary[2], 1]  # cingular boundary

    #---------------build the matrices A and C in the paper
    lon_weights = np.zeros(Nbv)
    C_lon = np.zeros(Nbv)#sparse.lil_matrix(Nbv, 1)
    A_lon_diag = np.zeros(Nbv)
    for lon_cstr_ind in sulcalCstr.longitudeCstrIndex:
        if sulcalCstr.sulcalLines[lon_cstr_ind].axisID in modele.longitudeAxisID:
            if modele.longitudeAxisCoord[modele.longitudeAxisID.index(sulcalCstr.sulcalLines[lon_cstr_ind].axisID)] is not None:
                print(sulcalCstr.sulcalLines[lon_cstr_ind].name+' weight in HOP energy : '+str(sulcalCstr.sulcalLines[lon_cstr_ind].weight))
                lon_weights[sulcalCstr.sulcalLines[lon_cstr_ind].indices] = sulcalCstr.sulcalLines[lon_cstr_ind].vertexWeight * sulcalCstr.sulcalLines[lon_cstr_ind].weight
                C_lon[sulcalCstr.sulcalLines[lon_cstr_ind].indices] = modele.longitudeAxisCoord[modele.longitudeAxisID.index(sulcalCstr.sulcalLines[lon_cstr_ind].axisID)] * lon_weights[sulcalCstr.sulcalLines[lon_cstr_ind].indices]
                A_lon_diag[sulcalCstr.sulcalLines[lon_cstr_ind].indices] = -lon_weights[sulcalCstr.sulcalLines[lon_cstr_ind].indices]
    A_lon = sparse.dia_matrix((A_lon_diag, 0), (Nbv, Nbv))

    lat_weights = np.zeros(Nbv)
    C_lat = np.zeros(Nbv)
    A_lat_diag = np.zeros(Nbv)
    for lat_cstr_ind in sulcalCstr.latitudeCstrIndex:
        if sulcalCstr.sulcalLines[lat_cstr_ind].axisID in modele.latitudeAxisID:
            if modele.latitudeAxisCoord[modele.latitudeAxisID.index(sulcalCstr.sulcalLines[lat_cstr_ind].axisID)] is not None:
                print(sulcalCstr.sulcalLines[lat_cstr_ind].name+' weight in HOP energy : '+str(sulcalCstr.sulcalLines[lat_cstr_ind].weight))
                lat_weights[sulcalCstr.sulcalLines[lat_cstr_ind].indices] = sulcalCstr.sulcalLines[lat_cstr_ind].vertexWeight * sulcalCstr.sulcalLines[lat_cstr_ind].weight
                C_lat[sulcalCstr.sulcalLines[lat_cstr_ind].indices] = modele.latitudeAxisCoord[modele.latitudeAxisID.index(sulcalCstr.sulcalLines[lat_cstr_ind].axisID)] * lat_weights[sulcalCstr.sulcalLines[lat_cstr_ind].indices]
                A_lat_diag[sulcalCstr.sulcalLines[lat_cstr_ind].indices] = -lat_weights[sulcalCstr.sulcalLines[lat_cstr_ind].indices]
    A_lat = sparse.dia_matrix((A_lat_diag, 0), (Nbv, Nbv))

##    t1 = timeit.default_timer()
##    Rx = sparse.csr_matrix(Rx)
##    t2 = timeit.default_timer()
##    print("%.2f sec"  %(t2-t1))
##    t1 = timeit.default_timer()
##    Rx = sparse.lil_matrix(Rx).tocsr()
##    t2 = timeit.default_timer()
##    print("%.2f sec"  %(t2-t1))

#     li = np.hstack(Lx.data)
#     print('nb Nan in Lx : ', len(np.where(np.isnan(li))[0]))
#     print('nb Inf in Lx : ', len(np.where(np.isinf(li))[0])   )
#     print('nb Nan in Rx : ', len(np.where(np.isnan(Rx))[0]))
#     print('nb Inf in Rx : ', len(np.where(np.isinf(Rx))[0]))
#     li = np.hstack(A_lon.data)
#     print('nb Nan in A_lon : ', len(np.where(np.isnan(li))[0]) )
#     print('nb Inf in A_lon : ', len(np.where(np.isinf(li))[0])  )
#     print('nb Nan in C_lon : ', len(np.where(np.isnan(C_lon))[0])  )
#     print('nb Inf in C_lon : ', len(np.where(np.isinf(C_lon))[0])  )
#      
#     li = np.hstack(Ly.data)
#     print('nb Nan in Ly : ', len(np.where(np.isnan(li))[0]))
#     print('nb Inf in Ly : ', len(np.where(np.isinf(li))[0]))
#     print('nb Nan in Ry : ', len(np.where(np.isnan(Ry))[0]))
#     print('nb Inf in Ry : ', len(np.where(np.isinf(Ry))[0])  )
#     li = np.hstack(A_lat.data)
#     print('nb Nan in A_lat : ', len(np.where(np.isnan(li))[0]) )
#     print('nb Inf in A_lat : ', len(np.where(np.isinf(li))[0])    )
#     print('nb Nan in C_lat : ', len(np.where(np.isnan(C_lat))[0])   )
#     print('nb Inf in C_lat : ', len(np.where(np.isinf(C_lat))[0])   )

#
    print('solve the linear system')
    Lx = Lx.tocsr()
    Ly = Ly.tocsr()
#    Rx = sparse.csr_matrix(Rx)
#    Ry = sparse.csr_matrix(Ry)
#    A_lon = A_lon.tocsr()
#    C_lon = sparse.csr_matrix(C_lon)
#    A_lat = A_lat.tocsr()
#    C_lat = sparse.csr_matrix(C_lat)
#    print(sparse.issparse(A_lat))
#    print(sparse.issparse(C_lat))

#    x = spsolve(Lx - cstrBalance * A_lon, cstrBalance * C_lon + Rx)
#    y = spsolve(Ly - cstrBalance * A_lat, cstrBalance * C_lat + Ry)
    
    



    t0 = time.perf_counter()
    if solver == 'lgmres':
        x, info = lgmres(Lx - cstrBalance * A_lon, cstrBalance * C_lon + Rx, tol=solver_tolerance)
        print(time.perf_counter() - t0, "seconds process time for x")
        t0 = time.perf_counter()
        y, info = lgmres(Ly - cstrBalance * A_lat, cstrBalance * C_lat + Ry, tol=solver_tolerance)
        print(time.perf_counter() - t0, "seconds process time for y")
    else:        
        x, info = gmres(Lx - cstrBalance * A_lon, cstrBalance * C_lon + Rx, tol=solver_tolerance)
        print(time.perf_counter() - t0, "seconds process time for x")
        t0 = time.perf_counter()
        y, info = gmres(Ly - cstrBalance * A_lat, cstrBalance * C_lat + Ry, tol=solver_tolerance)
        print(time.perf_counter() - t0, "seconds process time for y")
    print('matrix inverted')
    z = np.zeros(Nbv)
    vv = aims.vector_POINT3DF()
    for i in range(Nbv):
        vv.append([x[i], y[i], z[i]])

    rect = aims.AimsTimeSurface_3_VOID()
    rect.vertex().assign(vv)
    rect.polygon().assign(mesh.polygon())
    rect.updateNormals()
    return(rect)


####################################################################
#
# USELESS
#
#  load a rectangular mesh and extract the boundaries from coordinates
#
####################################################################
# def readRectangularMesh(mesh_filename,tex_filename):
#     re = aims.Reader()
#     mesh = re.read(mesh_filename)
#     vert = np.array(mesh.vertex())
#     tex = re.read(tex_filename)
#     "tex_filename[0] gives neocortex indices"
#     neocortex_indices = np.where(tex[0].arraydata() > 0)[0]
#     boundary = []
#     boundary.append(indsToROI(neocortex_indices, tex[1].arraydata()))  # insula boundary
#     boundary.append(indsToROI(neocortex_indices, tex[2].arraydata()))
#     boundary.append(indsToROI(neocortex_indices, tex[3].arraydata()))  # cingular boundary
#     "last boundary can be obtained from boundary[1]"
#     nb_new_verts = len(boundary[1])
#     new_verts_inds = range(vert.shape[0] - nb_new_verts, vert.shape[0])
#     boundary.append(new_verts_inds)
# 
# #     m = np.amin(vert, axis = 0)
# #     M = np.amax(vert, axis = 0)
# #     inds_x_min = np.where(vert[:, 0] == m[0])[0]
# #     inds_x_max = np.where(vert[:, 0] == M[0])[0]
# #     inds_y_min = np.where(vert[0, :] == m[1])[0]
# #     inds_y_max = np.where(vert[0, :] == M[1])[0]
# #     boundary = []
# #     boundary.append(inds_y_max)  # insula boundary
# #     boundary.append(inds_x_min)
# #     boundary.append(inds_y_min)  # cingular boundary
# #     boundary.append(inds_x_max)
#     print('boundary not good')
#     return(mesh, boundary, neocortex_indices)


####################################################################
#
# compute the coordinate textures for a mapped mesh
# from the coordinates of the vertices in the parameterized space (disk, rectangle or sphere)
#
#
####################################################################
def coordinatesFromRect(vert, poles_lat_insula, poles_lat_cingular):
    #vert = np.array(mesh.vertex())
    # mesh vertices to square [0 1] [0 1]
    m = vert.min(0)
    M = vert.max(0)

    m_lon = m[0]
    M_lon = M[0]
    m_lat = m[1]
    M_lat = M[1]

    lon = vert[:, 0] / (M_lon - m_lon)
    lat = vert[:, 1] / (M_lat - m_lat)

    lon = 360 * lon
    lat = (180 - poles_lat_insula - poles_lat_cingular) * lat + poles_lat_insula
    negs = np.where(lon < 0)
    lon[negs] = lon[negs] + 360
    return(lon, lat)


####################################################################
#
# compute the coordinate textures for a mapped mesh
# from the coordinates of the vertices in the parameterized disk
# lat is [0, 1]
# lon is [0, 360]
#
####################################################################
def coordinatesFromDisk(mesh):
    vert = np.array(mesh.vertex())
    r = np.sqrt(vert[:, 0] * vert[:, 0] + vert[:, 1] * vert[:, 1])
    lat = r / np.max(r) ## lat is always in [0, 1]
    lon = 2 * np.arctan(vert[:, 1] / (vert[:, 0] + r)) * 180 / np.pi + 180
    return(lon, lat)


####################################################################
#
# compute the coordinate textures for a mapped mesh
#
#
####################################################################
def computeCoordinates(nb_vert_full_mesh, neocortex_indices, rect_mesh, rect_mesh_boundary, poles_lat_insula, poles_lat_cingular, insula_mesh=None, insula_boundary=None, insula_indices=None, cingular_mesh=None, cingular_boundary=None, cingular_indices=None):

    rect_lon, rect_lat = coordinatesFromRect(np.array(rect_mesh.vertex()), poles_lat_insula, poles_lat_cingular)
    "rect_mesh_boundary[3] == new vertices always from insula to cingular pole"
    rect_lon = np.delete(rect_lon, rect_mesh_boundary[3], None)
    rect_lat = np.delete(rect_lat, rect_mesh_boundary[3], None)
    ar_neocortex_indices = np.array(neocortex_indices)
    lon = np.zeros(nb_vert_full_mesh)
    lon[ar_neocortex_indices] = rect_lon
    lat = np.zeros(nb_vert_full_mesh)
    lat[ar_neocortex_indices] = rect_lat
    return(lon, lat)

####################################################################
#
# compute shortest geodesic path between two array of indexes
# on a mesh
#
####################################################################
def getShortestPath(mesh, ind1, ind2):
    gp = aims.GeodesicPath(mesh, 3, 0)
    v = aims.vector_U32(ind2)
    lmin = 1000
    imin = 0
    jmin = 1
    for i in aims.vector_U32(ind1):
        j, l = gp.shortestPath_1_N_ind(i, v)
        if (l < lmin):
            lmin = l
            imin = i
            jmin = j
    return(gp.shortestPath_1_1_ind(imin, jmin).arraydata())




def path2Boundary(neoCortex_mesh, neoCortex_boundary, neocortex_poles_path, neigh=None):
    if neigh is None:
        neigh = aims.SurfaceManip.surfaceNeighbours(neoCortex_mesh)
    vert = np.array(neoCortex_mesh.vertex())

    neoCortex_open_boundary = boundaryReordering(neoCortex_boundary, neocortex_poles_path, vert)

    neigh_verts = set(np.hstack(list(neigh[v]) for v in neoCortex_open_boundary[1]))
    other_verts = neigh_verts.difference(neoCortex_open_boundary[1])
#     tex_out = aims.TimeTexture_S16()
#     tex_out[0].reserve(neoCortex_mesh.vertex().size())  # pre-allocates memory
#     for i in six.moves.xrange(neoCortex_mesh.vertex().size()):
#         if i in other_verts:
#             tex_out[0].append(1)
#         else:
#             tex_out[0].append(0)
#     ws = aims.Writer()
#     ws.write(tex_out, '/home/toz/ammon_Lwhite_neocortex_other_verts.tex')

    "group other verts into two connected sets, one on each side of the poles_path"
    l_other_verts = list(other_verts)
    tag = np.zeros(len(other_verts))
    nb_tagged = 1
    nb_tagged_o = 0
    cluster1 = [l_other_verts[0]]
    tag[l_other_verts.index(cluster1[0])] = 1
    neigh_cluster1 = set(neigh[cluster1[0]])
    while nb_tagged > nb_tagged_o:
        nb_tagged_o = nb_tagged
        for v in other_verts:
            if (v in neigh_cluster1) & (tag[l_other_verts.index(v)] == 0):
                cluster1.append(v)
                neigh_cluster1 = neigh_cluster1.union(neigh[v])
                tag[l_other_verts.index(v)] = 1
        nb_tagged = tag.sum()
        test = nb_tagged > nb_tagged_o
    "identify the anterior bank of the cut :: first vertex of the insula boundary is anterior while last one is posterior"   
    inter_bound0 = set(cluster1).intersection(neoCortex_open_boundary[0])
    #print('inter_bound0 ',inter_bound0)
    if len(inter_bound0) > 0:
        if neoCortex_open_boundary[0].index(list(inter_bound0)[0]) < (len(neoCortex_open_boundary[0]) / 2):
            posterior_cluster = other_verts.difference(cluster1)
        else:
            posterior_cluster = cluster1
    else:
        print('problem: cluster1.intersection(insula) is empty')

    "add the new vertices and poly to the mesh"
    poly = np.array(neoCortex_mesh.polygon())

    posterior_cluster_ind_poly = np.hstack(np.where(poly == i)[0] for i in posterior_cluster)
    poles_path_ind_poly = np.hstack(np.where(poly == i)[0] for i in neoCortex_open_boundary[1])

    added_poly_inds = set(posterior_cluster_ind_poly).intersection(poles_path_ind_poly)
    added_poly = poly[list(added_poly_inds), :]

    nb_new_verts = len(neoCortex_open_boundary[1])
    new_verts_inds = list(range(vert.shape[0], vert.shape[0] + nb_new_verts))
    for i in range(nb_new_verts):
        places = np.where(added_poly == neoCortex_open_boundary[1][i])
        added_poly[places[0], places[1]] = new_verts_inds[i]
    #poly = np.vstack([poly,added_poly])
    poly[list(added_poly_inds), :] = added_poly
    pp = aims.vector_AimsVector_U32_3()
    pp.assign([aims.AimsVector_U32_3([poly[i, 0], poly[i, 1], poly[i, 2]]) for i in range(poly.shape[0])])

    vert = np.vstack([vert, vert[neoCortex_open_boundary[1], :]])
    vv = aims.vector_POINT3DF()
    for i in range(vert.shape[0]):
        vv.append([vert[i, 0], vert[i, 1], vert[i, 2]])

    neoCortex_open_mesh = aims.AimsTimeSurface_3_VOID()
    neoCortex_open_mesh.vertex().assign(vv)
    neoCortex_open_mesh.polygon().assign(pp)
    neoCortex_open_mesh.updateNormals()

    "add the new boundary and connect it"
    "neoCortex_open_boundary[3] == new verices always from insula to cingular pole as for neoCortex_open_boundary[1]"
    neoCortex_open_boundary.append(new_verts_inds)
    neoCortex_open_boundary[0].append(new_verts_inds[0])
    neoCortex_open_boundary[2].append(new_verts_inds[-1])
#
#    vv = aims.vector_POINT3DF()
#    ee = aims.vector_AimsVector_U32_3()
#        for i in boundary[bound_ind]:
#            vv.append([vert[i,0], vert[i,1],vert[i,2]])
#        for i in range(len(boundary[bound_ind])-1):
#            ee.append(np.array([i,i+1],np.uint32))
#        bound_mesh.vertex(bound_ind).assign(vv)
#        bound_mesh.polygon(bound_ind).assign(ee)
#        pol = np.vstack( (np.zeros( len(boundary[bound_ind])-2, dtype = np.int32 ),boundary[bound_ind][0:-2],boundary[bound_ind][1:-1] )).transpose()
#        print(pol.dtype#astype(np.float32).dtype)
#        bound_mesh.polygon(bound_ind).assign([ aims.AimsVector(x,'U32') for x in pol ])

    return (neoCortex_open_mesh, neoCortex_open_boundary)

####################################################################
#
#    ensures neoCortex_boundary[0] == insula_boundary and boundary[2] == cingular_boundary
#    neoCortex_boundary[1] == neocortex_poles_path
#    insula and cingular boundary are always oriented in the same way
#
#    neocortex_poles_path given as parameter is assumed to be always oriented from insula to cingular pole
#    it is the orientation of neocortex_poles_path that drives the orientation of insula and cingular boundary
####################################################################
def boundaryReordering(neoCortex_boundary, neocortex_poles_path, vert):

    path_bound0 = list(set(neocortex_poles_path).intersection(set(neoCortex_boundary[0])))
    path_bound1 = list(set(neocortex_poles_path).intersection(set(neoCortex_boundary[1])))
    if not path_bound0:
        print("problem: no intersection between boundary 0 and poles_path")
    if not path_bound1:
        print("problem: no intersection between boundary 1 and poles_path")


    "ensures a single vertex intersection between poles_path and boundaries"
    "note that neocortex_poles_path is the ordered list of vertices, "
    "i.e if the intersection consists of more than 1 vertex, the supplementary vertices to suppress are at the begining or at the end but not in the middle"
    if len(path_bound0) > 1:
        path_length = len(neocortex_poles_path)
        inds_path_bound0 = [neocortex_poles_path.index(i) for i in path_bound0]
        m_ind = min(inds_path_bound0)
        M_ind = max(inds_path_bound0)
        if M_ind < (path_length/2.): # delete the first items
            removeset = set(range(M_ind))
        else:# delete the last items
            removeset = set(range(m_ind+1,path_length))
        neocortex_poles_path = [v for i, v in enumerate(neocortex_poles_path) if i not in removeset]
        path_bound0 = list(set(neocortex_poles_path).intersection(set(neoCortex_boundary[0])))
    if len(path_bound1) > 1:
        path_length = len(neocortex_poles_path)
        inds_path_bound1 = [neocortex_poles_path.index(i) for i in path_bound1]
        m_ind = min(inds_path_bound1)
        M_ind = max(inds_path_bound1)
        if M_ind < (path_length/2.): # delete the first items
            removeset = set(range(M_ind))
        else:# delete the last items
            removeset = set(range(m_ind+1,path_length))
        neocortex_poles_path = [v for i, v in enumerate(neocortex_poles_path) if i not in removeset]
        path_bound1 = list(set(neocortex_poles_path).intersection(set(neoCortex_boundary[1])))
    ind_path_bound0 = neocortex_poles_path.index(path_bound0[0])
    ind_path_bound1 = neocortex_poles_path.index(path_bound1[0])
    "order the boundaries and rotate insula and cingular boundary if necessary"
    rotate = lambda l, n :l[n:] + l[:n]  # list rotation fct
    rt_bound0 = neoCortex_boundary[0].index(path_bound0[0])
    rt_bound1 = neoCortex_boundary[1].index(path_bound1[0])

    ordered_neoCortex_boundary = list()
    if ind_path_bound0 < ind_path_bound1:  # neoCortex_boundary[0] = =  insula boundary
        ordered_neoCortex_boundary.append(rotate(neoCortex_boundary[0], rt_bound0))
        ordered_neoCortex_boundary.append(neocortex_poles_path)
        ordered_neoCortex_boundary.append(rotate(neoCortex_boundary[1], rt_bound1))
    else:  # neoCortex_boundary[1] = =  insula boundary
        ordered_neoCortex_boundary.append(rotate(neoCortex_boundary[1], rt_bound1))
        ordered_neoCortex_boundary.append(neocortex_poles_path)
        ordered_neoCortex_boundary.append(rotate(neoCortex_boundary[0], rt_bound0))

    ###########################################
    # inversion of insula or cingular boundary: not good to use coordinates of the vertices here
    ###########################################
    "reverse insula and cingular boundary if necessary"
#    sec_coord_insula = vert[ordered_neoCortex_boundary[0],1]
#    print(sec_coord_insula)
#    print((np.argmin(sec_coord_insula),sec_coord_insula[np.argmin(sec_coord_insula)]))

    if np.argmin(vert[ordered_neoCortex_boundary[0], 1]) < (len(ordered_neoCortex_boundary[0]) / 2):
        ordered_neoCortex_boundary[0].reverse()
        ordered_neoCortex_boundary[0] = rotate(ordered_neoCortex_boundary[0], -1)
        print('reversing insula boundary')
    if np.argmin(vert[ordered_neoCortex_boundary[2], 1]) < (len(ordered_neoCortex_boundary[2]) / 2):
        ordered_neoCortex_boundary[2].reverse()
        ordered_neoCortex_boundary[2] = rotate(ordered_neoCortex_boundary[2], -1)
        print('reversing cingular boundary')
    return ordered_neoCortex_boundary
#    if len(neoCortex_boundary[0]) == len(insula_boundary[0]):
#        if len(neoCortex_boundary[1]) == len(cingular_boundary[0]):
#            print('bounrary correspondance validated')
#        elif len(neoCortex_boundary[1]) == len(insula_boundary[0]):
#            print('all boundaries have the same length.....')
#        else:
#            print('problem in the cinular boundary')
#    elif len(neoCortex_boundary[0]) == len(cingular_boundary[0]):
#        if len(neoCortex_boundary[1]) == len(insula_boundary[0]):
#            print('bounrary correspondance inverted -> reversing boundaries')
#            tmp_bound = neoCortex_boundary
#            neoCortex_boundary[0] = tmp_bound[1]
#            neoCortex_boundary[1] = tmp_bound[0]
#        else:
#            print('problem in the insular boundary')
#    else:
#        print('problem in both insular and cingular boundary')
#    return neoCortex_boundary


def indsToROI(neocortex_indices, inds):
    neocortex_inds = [neocortex_indices.index(i) for i in inds]
    return neocortex_inds


def catConstraints(sulci_list):
    union_constraints = sulci_list[0]
    if len(sulci_list) > 1:
        for sulc in sulci_list[1:-1]:
            union_constraints.cat(sulc)
    return union_constraints


def texture2ROI(tex_cstr, neocortex_indices):
    if tex_cstr is np.ndarray:
        return tex_cstr[np.array(neocortex_indices)]
    else:
        return np.array(tex_cstr)[np.array(neocortex_indices)]


####################################################################
#
# compute the spherical coordinates of vertices from lon and lat textures
#
####################################################################
def sphericalMeshFromCoords(tex_lat, tex_lon, ray):
    tex_lat = tex_lat * np.pi / 180
    tex_lon = tex_lon * np.pi / 180
#    print(tex_lon)
    spherical_verts = np.ndarray((tex_lat.shape[0], 3))
    spherical_verts[:, 1] = ray * np.cos(tex_lon) * np.sin(tex_lat)
    spherical_verts[:, 0] = ray * np.sin(tex_lon) * np.sin(tex_lat)
    spherical_verts[:, 2] = ray * np.cos(tex_lat)
#    print(spherical_verts[:, 2])
    return spherical_verts



####################################################################
#
# identify the inverted triangles in a parameterized mesh
#
####################################################################
def invertedPolygon(mesh, shape=None):
    if shape is None:
        shape = 'square'
#    norms = np.array(mesh.normal())
    norms = basicTls.meshPolygonNormal(mesh)
    if shape is 'sphere':
        print('sphere not available yet')
#N(1,:)=FaceNormal(1,:)+FaceCenter(1,:);
#N(2,:)=FaceNormal(2,:)+FaceCenter(2,:);
#N(3,:)=FaceNormal(3,:)+FaceCenter(3,:);
#nb_inward=0;
#ind=1;
#inward=[];
#    for i=1:length(N(1,:))
#        if norm(N(:,i))<1
#            nb_inward=nb_inward+1;
#            inward(ind)=i;
#            ind=ind+1;
#        end
#    end
    else:
        inward = np.where(norms[:, 2] > 0)[0]
        nb_inward = len(inward)
        if nb_inward > (len(mesh.polygon()) / 2):
            inward = np.where(norms[:, 2] < 0)[0]
            nb_inward = len(inward)
    return(nb_inward, inward)

def polygonNormals(vertices, faces):
    print('vertices.shape',vertices.shape)
    print('faces.shape',faces.shape)
    print(faces[:,2].shape)
    print(vertices[faces[:,2],:].shape)
    normalf = crossp( vertices[faces[:,1],:]-vertices[faces[:,0],:], vertices[faces[:,2],:]-vertices[faces[:,0],:] )
    print(normalf.shape)
    d = np.sqrt(np.sum(np.square(normalf),1))
    d[d<0.0000001]=1
    print(d.shape)
    rep =  np.tile( d, (3,1) ).transpose()
    print(rep.shape)
    normalf = normalf / rep
    return normalf

def crossp(x,y):
    # x and y are (m,3) dimensional
    z = x.copy()
    z[:,0] = x[:,1]*y[:,2] - x[:,2]*y[:,1]
    z[:,1] = x[:,2]*y[:,0] - x[:,0]*y[:,2]
    z[:,2] = x[:,0]*y[:,1] - x[:,1]*y[:,0]
    return z

####################################################################
#
# reverse properly the 4 boundaries of a rectangular mesh
#
####################################################################
def recantgleBoundaryFlip(neoCortex_boundary):
    rotate = lambda l, n :l[n:] + l[:n]  # list rotation fct
    fliped_neoCortex_boundary = []
    fliped_neoCortex_boundary.append(neoCortex_boundary[0])
    fliped_neoCortex_boundary.append(neoCortex_boundary[3])
    fliped_neoCortex_boundary.append(neoCortex_boundary[2])
    fliped_neoCortex_boundary.append(neoCortex_boundary[1])
    
    fliped_neoCortex_boundary[0].reverse()
    fliped_neoCortex_boundary[0] = rotate(fliped_neoCortex_boundary[0], -1)
    fliped_neoCortex_boundary[2].reverse()
    fliped_neoCortex_boundary[2] = rotate(fliped_neoCortex_boundary[2], -1)
     
    return fliped_neoCortex_boundary

####################################################################
#
# flip a rectangular mesh following the first axis (vert[:,0])
#
####################################################################
def recantgleFlip(rect_mesh):
    vert = np.array(rect_mesh.vertex())
    # mesh vertices to square [0 1] [0 1]
    m = vert.min(0)
    M = vert.max(0)

    vert[:, 0] = -vert[:,0] + (M[0] - m[0])
    m2 = vert.min(0)
    M2 = vert.max(0)
    vv = aims.vector_POINT3DF()
    for x in vert:
         vv.append(x)
    rect_mesh.vertex().assign(vv)
    return rect_mesh

def parcelsFromCoordinates(template_lat,template_lon,model,
                           parcellation_type=None,
                           between_poles_parcel_anterior_coord=130,
                           between_poles_parcel_posterior_coord=250):
    if parcellation_type is None:
        parcellation_type == 'model' # default resolution
    
    (longitude_axis_coords, latitude_axis_coords) = model.axisCoordToDegree()
#    print(longitude_axis_coords)
#    print(latitude_axis_coords)
    
    nb_vert = template_lat.shape[0]
    tex_parcels = np.zeros(nb_vert)
    lab_parcel = 2 #1 will correspond to the cingular pole
    sort_axes_lon = [360]
    for f in longitude_axis_coords[:-1]:#[360 - f for f in model.longitudeAxisID]
        if f is not None:
            sort_axes_lon.append(f)
    # additional longitude axes that delimitate the 'trash' region between
    # the insular and cingular poles
    '''
    old fashioned: a central value and a  width parameter
    between_poles_parcel_central_value = 190
    between_poles_parcel_width = 120
    sort_axes_lon.append(between_poles_parcel_central_value - between_poles_parcel_width / 2)
    sort_axes_lon.append(between_poles_parcel_central_value + between_poles_parcel_width / 2)
    '''
    sort_axes_lon.append(between_poles_parcel_anterior_coord)
    sort_axes_lon.append(between_poles_parcel_posterior_coord)

    sort_axes_lat = [0, model.insularPoleBoundaryCoord]
    for f in latitude_axis_coords:
        if f is not  None:
           sort_axes_lat.append(f) 
         
    sort_axes_lon.sort()
    sort_axes_lat.sort()
    sort_axes_lat.append(180-model.cingularPoleBoundaryCoord)

    # add suplementary axes to the model <=> subdivise parcels
    if parcellation_type == 'marsAtlas' or parcellation_type == 'model' :
        # antero-posterior subdivision of the prefrontal lobe
        sort_axes_lon.append(sort_axes_lon[1] + (sort_axes_lon[2] - sort_axes_lon[1]) / 2)
        # antero-posterior subdivision of the temporal lobe
        sort_axes_lon.append(sort_axes_lon[5] + 3 * (sort_axes_lon[6] - sort_axes_lon[5]) / 4)
#        sort_axes_lon.append(sort_axes_lon[5] + temporal_pole_parcel_width)
        sort_axes_lon.sort()

    elif parcellation_type == 'heschl' :
        # antero-posterior subdivision of the prefrontal lobe
        sort_axes_lon.append(sort_axes_lon[1] + (sort_axes_lon[2] - sort_axes_lon[1]) / 2)
        sort_axes_lon.sort()

    elif parcellation_type == 'lobes':
        # antero-posterior subdivision of the temporal lobe
        sort_axes_lon.append(sort_axes_lon[5] + 7 * (sort_axes_lon[6] - sort_axes_lon[5]) / 8)
#        sort_axes_lon.append(sort_axes_lon[5] + temporal_pole_parcel_width)
        del sort_axes_lon[8]
        #del sort_axes_lon[6]
        del sort_axes_lon[3]
        del sort_axes_lon[2]
        del sort_axes_lon[1]
        sort_axes_lon.sort()

        del sort_axes_lat[5]
        #del sort_axes_lat[4]
        del sort_axes_lat[3]
        del sort_axes_lat[2]

    elif parcellation_type == 'gyri':
        del sort_axes_lat[3]

    for t_lon in range(len(sort_axes_lon)-1):
#        print(sort_axes_lon[t_lon])
        inds_lon = np.where((template_lon >= sort_axes_lon[t_lon])&(template_lon<=sort_axes_lon[t_lon+1]))[0]
        for t_lat in range(len(sort_axes_lat)-1):
#            print(sort_axes_lat[t_lat])
            inds_lat = np.where((template_lat[inds_lon]>=sort_axes_lat[t_lat])&(template_lat[inds_lon]<=sort_axes_lat[t_lat+1]))[0]
            tex_parcels[inds_lon[inds_lat]]=lab_parcel
#            print('lab_parcel', lab_parcel)
            lab_parcel = lab_parcel+1



    if parcellation_type == 'model' or parcellation_type == 'heschl' :
        # parcels merging
        # INSULA
        tex_parcels[tex_parcels == 9] = 44
        tex_parcels[tex_parcels == 2] = 44
        tex_parcels[tex_parcels == 16] = 44
        tex_parcels[tex_parcels == 23] = 44
        tex_parcels[tex_parcels == 30] = 44
        tex_parcels[tex_parcels == 65] = 44
        tex_parcels[tex_parcels == 72] = 44
        tex_parcels[tex_parcels == 51] = 44
        tex_parcels[tex_parcels == 58] = 44
        tex_parcels[tex_parcels == 37] = 44
        # arround the path between the poles
        tex_parcels[tex_parcels == 38] = 255
        tex_parcels[tex_parcels == 39] = 255
        tex_parcels[tex_parcels == 40] = 255
        tex_parcels[tex_parcels == 41] = 255
        tex_parcels[tex_parcels == 42] = 255
        tex_parcels[tex_parcels == 43] = 255
        # cingular pole
        tex_parcels[tex_parcels == 0] = 1

        (tex_parcels, nb_parcels) = reorganize_parcels(tex_parcels)
        # ensure path between poles parcel has a label of 255
        m = np.max(tex_parcels)
        tex_parcels[tex_parcels == m] = 255

        ## old: original model, without additional axes
        # # INSULA
        # tex_parcels[tex_parcels == 9] = 2
        # tex_parcels[tex_parcels == 16] = 2
        # tex_parcels[tex_parcels == 23] = 2
        # tex_parcels[tex_parcels == 37] = 2
        # tex_parcels[tex_parcels == 44] = 2
        # tex_parcels[tex_parcels == 51] = 2
        # tex_parcels[tex_parcels == 58] = 2
        # # arround the path between the poles
        # tex_parcels[tex_parcels == 31] = 30
        # tex_parcels[tex_parcels == 32] = 30
        # tex_parcels[tex_parcels == 33] = 30
        # tex_parcels[tex_parcels == 34] = 30
        # tex_parcels[tex_parcels == 35] = 30
        # tex_parcels[tex_parcels == 36] = 30
        # # cingular pole
        # tex_parcels[tex_parcels == 0] = 1
        #  (tex_parcels, nb_parcels) = reorganize_parcels(tex_parcels)

    elif parcellation_type == 'marsAtlas':
        # parcels merging, parcel names correspond to the nomenclature given in the paper
        # INSULA sup ant
        # tex_parcels[tex_parcels == 16] = 9
        # tex_parcels[tex_parcels == 23] = 9
        # tex_parcels[tex_parcels == 30] = 9
        # INSULA sup post
        # tex_parcels[tex_parcels == 65] = 2
        # tex_parcels[tex_parcels == 72] = 2
        # INSULA inf
        # tex_parcels[tex_parcels == 51] = 44
        # tex_parcels[tex_parcels == 58] = 44
        # INSULA
        # tex_parcels[tex_parcels == 2] = 44  # 3] = 69
        # tex_parcels[tex_parcels == 9] = 44  # 17] = 69

        tex_parcels_tmp = tex_parcels.copy()
        # INSULA
        tex_parcels_tmp[tex_parcels == 9] = 44
        tex_parcels_tmp[tex_parcels == 2] = 44
        tex_parcels_tmp[tex_parcels == 16] = 44
        tex_parcels_tmp[tex_parcels == 23] = 44
        tex_parcels_tmp[tex_parcels == 30] = 44
        tex_parcels_tmp[tex_parcels == 65] = 44
        tex_parcels_tmp[tex_parcels == 72] = 44
        tex_parcels_tmp[tex_parcels == 51] = 44
        tex_parcels_tmp[tex_parcels == 58] = 44
        tex_parcels_tmp[tex_parcels == 37] = 44
        # arround the path between the poles
        tex_parcels_tmp[tex_parcels == 39] = 38
        tex_parcels_tmp[tex_parcels == 40] = 38
        tex_parcels_tmp[tex_parcels == 41] = 38
        tex_parcels_tmp[tex_parcels == 42] = 38
        tex_parcels_tmp[tex_parcels == 43] = 38
        # ITCm
        tex_parcels_tmp[tex_parcels == 48] = 49  # 77] = 79
        # VCs
        tex_parcels_tmp[tex_parcels == 60] = 61  # 97] = 99
        # Cu
        tex_parcels_tmp[tex_parcels == 62] = 63  # 101] = 103
        # ICC
        tex_parcels_tmp[tex_parcels == 57] = 64  # 81] = 105
        tex_parcels_tmp[tex_parcels == 50] = 64  # 93] = 105
        # PCC
        tex_parcels_tmp[tex_parcels == 71] = 78  # 117] = 129
        # MCC
        tex_parcels_tmp[tex_parcels == 15] = 22  # 15] = 41
        tex_parcels_tmp[tex_parcels == 8] = 22  # 29] = 41
        # IPCv
        tex_parcels_tmp[tex_parcels == 59] = 66  # 95] = 107
        # Sdl
        tex_parcels_tmp[tex_parcels == 74] = 75  # 121] = 123
        # Sdm
        tex_parcels_tmp[tex_parcels == 76] = 77  # 125] = 127
        # PFC vm
        tex_parcels_tmp[tex_parcels == 34] = 35  # 61] = 63
        # Mdl
        tex_parcels_tmp[tex_parcels == 4] = 5  # 7] = 9
        # Mdm
        tex_parcels_tmp[tex_parcels == 6] = 7  # 11] = 13
        # PMrV
        tex_parcels_tmp[tex_parcels == 10] = 17  # 19] = 31
        # PMdl
        tex_parcels_tmp[tex_parcels == 11] = 12  # 21] = 23
        # PMdm
        tex_parcels_tmp[tex_parcels == 13] = 14  # 25] = 27
        # PFcdl
        tex_parcels_tmp[tex_parcels == 18] = 19  # 33] = 35
        # PFcdm
        tex_parcels_tmp[tex_parcels == 20] = 21  # 37] = 39
        # ACC
        tex_parcels_tmp[tex_parcels == 29] = 36  # 53] = 65

        tex_parcels = np.zeros(tex_parcels_tmp.shape)
        # cingular pole = 0
        # tex_parcels_tmp[tex_parcels == 0] = 1
        # ensure path between poles parcel has a label of 255
        tex_parcels[tex_parcels_tmp == 38] = 255

        # reorganization of the labels to match the nomenclature of the paper
        #	VCcm
        tex_parcels[tex_parcels_tmp == 55] = 1
        #	VCl
        tex_parcels[tex_parcels_tmp == 54] = 2
        # VCs
        tex_parcels[tex_parcels_tmp == 61] = 3
        # Cu
        tex_parcels[tex_parcels_tmp == 63] = 4
        #	VCrm
        tex_parcels[tex_parcels_tmp == 56] = 5
        # ITCm
        tex_parcels[tex_parcels_tmp == 49] = 6
        #	ITCr
        tex_parcels[tex_parcels_tmp == 47] = 7
        #	MTCc
        tex_parcels[tex_parcels_tmp == 53] = 8
        #	STCc
        tex_parcels[tex_parcels_tmp == 52] = 9
        #	STCr
        tex_parcels[tex_parcels_tmp == 45] = 10
        #	MTCr
        tex_parcels[tex_parcels_tmp == 46] = 11
        # ICC
        tex_parcels[tex_parcels_tmp == 64] = 12
        # IPCv
        tex_parcels[tex_parcels_tmp == 66] = 13
        #	IPCd
        tex_parcels[tex_parcels_tmp == 67] = 14
        #	SPC
        tex_parcels[tex_parcels_tmp == 68] = 15
        #	SPCm
        tex_parcels[tex_parcels_tmp == 69] = 16
        #	PCm
        tex_parcels[tex_parcels_tmp == 70] = 17
        # PCC
        tex_parcels[tex_parcels_tmp == 78] = 18
        #	Sv
        tex_parcels[tex_parcels_tmp == 73] = 19
        # Sdl
        tex_parcels[tex_parcels_tmp == 75] = 20
        # Sdm
        tex_parcels[tex_parcels_tmp == 77] = 21
        #	Mv
        tex_parcels[tex_parcels_tmp == 3] = 22
        # Mdl
        tex_parcels[tex_parcels_tmp == 5] = 23
        # Mdm
        tex_parcels[tex_parcels_tmp == 7] = 24
        # PMrV
        tex_parcels[tex_parcels_tmp == 17] = 25
        # PMdl
        tex_parcels[tex_parcels_tmp == 12] = 26
        # PMdm
        tex_parcels[tex_parcels_tmp == 14] = 27
        # PFcdl
        tex_parcels[tex_parcels_tmp == 19] = 28
        # PFcdm
        tex_parcels[tex_parcels_tmp == 21] = 29
        # MCC
        tex_parcels[tex_parcels_tmp == 22] = 30
        #PFrvl
        tex_parcels[tex_parcels_tmp == 24] = 31
        #	Pfrdli
        tex_parcels[tex_parcels_tmp == 25] = 32
        #	Pfrdls
        tex_parcels[tex_parcels_tmp == 26] = 33
        #	PFrd
        tex_parcels[tex_parcels_tmp == 27] = 34
        #	PFrm
        tex_parcels[tex_parcels_tmp == 28] = 35
        #	OFCvl
        tex_parcels[tex_parcels_tmp == 31] = 36
        #	OFCv
        tex_parcels[tex_parcels_tmp == 32] = 37
        #	OFCvm
        tex_parcels[tex_parcels_tmp == 33] = 38
        # PFC vm
        tex_parcels[tex_parcels_tmp == 35] = 39
        # ACC
        tex_parcels[tex_parcels_tmp == 36] = 40
        # insula
        tex_parcels[tex_parcels_tmp == 44] = 41

        nb_parcels = np.unique(tex_parcels).shape[0]

    elif parcellation_type == 'lobes' :
        # parcels merging
        # INSULA
        tex_parcels[tex_parcels == 22] = 2
        tex_parcels[tex_parcels == 18] = 2
        tex_parcels[tex_parcels == 14] = 2
        tex_parcels[tex_parcels == 10] = 2
        tex_parcels[tex_parcels == 6] = 2
        # arround the path between the poles
        tex_parcels[tex_parcels == 7] = 255
        tex_parcels[tex_parcels == 8] = 255
        tex_parcels[tex_parcels == 9] = 255
        # limbic
        tex_parcels[tex_parcels == 25] = 5
        tex_parcels[tex_parcels == 21] = 5
        tex_parcels[tex_parcels == 17] = 5
        tex_parcels[tex_parcels == 13] = 5
        # frontal
        tex_parcels[tex_parcels == 4] = 3
        # parietal
        tex_parcels[tex_parcels == 24] = 19
        tex_parcels[tex_parcels == 23] = 19
        # temporal
        tex_parcels[tex_parcels == 15] = 11
        tex_parcels[tex_parcels == 12] = 11
        # occipital
        tex_parcels[tex_parcels == 20] = 16
        # cingular pole
        tex_parcels[tex_parcels == 0] = 1
        nb_parcels = np.unique(tex_parcels).shape[0]
        #(tex_parcels, nb_parcels) = reorganize_parcels(tex_parcels)
        # ensure path between poles parcel has a label of 255
        m = np.max(tex_parcels)
        tex_parcels[tex_parcels == m] = 255

    elif parcellation_type == 'gyri' :
        # parcels merging
        # INSULA
        tex_parcels[tex_parcels == 14] = 2
        tex_parcels[tex_parcels == 20] = 2
        tex_parcels[tex_parcels == 8] = 2
        tex_parcels[tex_parcels == 50] = 2
        tex_parcels[tex_parcels == 44] = 2
        tex_parcels[tex_parcels == 32] = 2
        tex_parcels[tex_parcels == 26] = 2
        tex_parcels[tex_parcels == 38] = 2
        # arround the path between the poles
        tex_parcels[tex_parcels == 27] = 255
        tex_parcels[tex_parcels == 28] = 255
        tex_parcels[tex_parcels == 29] = 255
        tex_parcels[tex_parcels == 30] = 255
        tex_parcels[tex_parcels == 31] = 255
        # precentral
        tex_parcels[tex_parcels == 4] = 3
        tex_parcels[tex_parcels == 5] = 3
        # postcentral
        tex_parcels[tex_parcels == 52] = 51
        tex_parcels[tex_parcels == 53] = 51
        # superior frontal
        tex_parcels[tex_parcels == 17] = 11
        tex_parcels[tex_parcels == 12] = 11
        tex_parcels[tex_parcels == 18] = 11
        # medial orbitofrontal
        tex_parcels[tex_parcels == 24] = 23
        # paracentral
        tex_parcels[tex_parcels == 54] = 6
        # caudal anterior cingulate
        tex_parcels[tex_parcels == 19] = 13
        # posterior cingulate
        tex_parcels[tex_parcels == 55] = 7
        # isthmus cingulate
        tex_parcels[tex_parcels == 49] = 43
        # parahipocampal
        tex_parcels[tex_parcels == 36] = 37
        # cuneus
        tex_parcels[tex_parcels == 42] = 41
        # supramarginal
        tex_parcels[tex_parcels == 45] = 39
        # inferior parietal
        tex_parcels[tex_parcels == 46] = 40

        (tex_parcels, nb_parcels) = reorganize_parcels(tex_parcels)
        # ensure path between poles parcel has a label of 255
        m = np.max(tex_parcels)
        tex_parcels[tex_parcels == m] = 255

    elif parcellation_type == 'model_foetus':
        # concatenate some parcels
        # INSULA
        tex_parcels[tex_parcels == 7] = 2
        tex_parcels[tex_parcels == 12] = 2
        tex_parcels[tex_parcels == 17] = 2
        tex_parcels[tex_parcels == 22] = 2
        tex_parcels[tex_parcels == 27] = 2
        tex_parcels[tex_parcels == 32] = 2
        tex_parcels[tex_parcels == 37] = 2
        # arround the path between the poles
        tex_parcels[tex_parcels == 19] = 18
        tex_parcels[tex_parcels == 20] = 18
        tex_parcels[tex_parcels == 21] = 18
        # post sup temporal gyrus
        tex_parcels[tex_parcels == 28] = 23
        # post inf temporal gyrus
        tex_parcels[tex_parcels == 29] = 24
        # int frontal gyrus
        tex_parcels[tex_parcels == 16] = 11
        # sup frontal gyrus
        tex_parcels[tex_parcels == 15] = 10
        # med frontal gyrus
        tex_parcels[tex_parcels == 14] = 9
        # inf frontal gyrus
        tex_parcels[tex_parcels == 13] = 8
        # precentral gyrus
        tex_parcels[tex_parcels == 4] = 3
        tex_parcels[tex_parcels == 5] = 3
        # postcentral gyrus
        tex_parcels[tex_parcels == 39] = 38
        tex_parcels[tex_parcels == 40] = 38
        # med cingular gyrus
        tex_parcels[tex_parcels == 41] = 6
        # precuneus
        tex_parcels[tex_parcels == 31] = 36
        # cingular pole
        tex_parcels[tex_parcels == 0] = 1

        (tex_parcels, nb_parcels) = reorganize_parcels(tex_parcels)

    else:  # unknown parcellation type
        nb_parcels = np.unique(tex_parcels).shape[0]

    return (tex_parcels, nb_parcels)


####################################################################
#
# reorganization of the labels in parcels texture in order to have values that follows each other (1 to nb_parcels)
#
####################################################################
def reorganize_parcels(tex_parcels):
    uparcells = np.unique(tex_parcels)
    tex_parcels_tmp = tex_parcels.copy()
    for ind,u_parc in enumerate(uparcells): # note that ind>=1
        tex_parcels[tex_parcels_tmp == u_parc] = ind
    return (tex_parcels, uparcells.shape[0])




####################################################################
#
# identify the inverted triangles in a parameterized mesh
#
####################################################################
def solveInvertedPolygon(mesh, boundary, nb_it_smooth, neigh=None):
    max_tot_count = 30
    max_count = 5
    if neigh is None:
        neigh = aims.SurfaceManip.surfaceNeighbours(mesh)
    poly = np.array(mesh.polygon())
    (nb_inward, inward) = invertedPolygon(mesh)
    nb_inward_evol = [nb_inward]
    inward_evol = [inward]
    count = 0
    tot_count = 0
    while nb_inward_evol[-1] > 0:
        # group connected polyons into patchs
#        t0 = time.perf_counter()
        bary = np.unique(poly[inward, :])
#         print(time.perf_counter() - t0, "seconds process time")
#         t0 = time.perf_counter()
#         adj_list = []
#         for b in bary:
#             adj_list.append(neigh[b].list())
#         print(time.perf_counter() - t0, "seconds process time")
        vert = np.array(mesh.vertex())
#        t0 = time.perf_counter()
        vert_out = vert.copy()
#        print(time.perf_counter() - t0, "seconds process time")
#        t0 = time.perf_counter()
        for it in range(nb_it_smooth):
            for ind_bary,b in enumerate(bary):#l_neigh_bary in enumerate(adj_list):
                l_neigh_bary = neigh[b].list()
                for neigh_bary in l_neigh_bary:
                    neig = vert[neigh[neigh_bary].list(), :]
                    vert_out[neigh_bary, :] = np.mean(neig, axis = 0)
                #g_neig = vert_out[l_neigh_bary, :]
                vert_out[bary[ind_bary], :] = np.mean(vert_out[l_neigh_bary, :], axis = 0)#g_neig, axis = 0)
            ## re projection of bounday vertices onto the bounday
            vert_out[boundary[0], 1] = vert[boundary[0], 1]
            vert_out[boundary[1], 0] = vert[boundary[1], 0]
            vert_out[boundary[2], 1] = vert[boundary[2], 1]
            vert_out[boundary[3], 0] = vert[boundary[3], 0]
            ###########################
            vert = vert_out.copy()
#        print(time.perf_counter() - t0, "seconds process time")
#        t0 = time.perf_counter()
        vv = aims.vector_POINT3DF()
        for x in vert:
            vv.append(x)
        mesh.vertex().assign(vv)
#        print(time.perf_counter() - t0, "seconds process time")
#        t0 = time.perf_counter()
        (nb_inward, inward) = invertedPolygon(mesh)
#        print(time.perf_counter() - t0, "seconds process time")
        nb_inward_evol.append(nb_inward)
        inward_evol.append(inward)
        print(nb_inward_evol)
        tot_count +=1
        if tot_count > max_tot_count:
            print('unable to solve the inverted faces')
            break
        if len(nb_inward_evol)>3:
            if nb_inward_evol[-1] == nb_inward_evol[-2] or nb_inward_evol[-1] == nb_inward_evol[-3]:
                count += 1
            else:
                count = 0
        if count > max_count:
            print('unable to solve the inverted faces')
            break
        mesh.updateNormals()
    return (mesh, nb_inward_evol, inward_evol)




####################################################################
#
# HIP with multiple poles_path
#
####################################################################
def hip_multi_path(mesh, insula_tex_clean, cingular_tex_clean, length, width,Npath,step,path_displacement):
#    square_ratio = 4.5
    neocortex_tex_value = 0
    insula_tex_value = 180
    cingular_tex_value = 1
    #write_all_steps_to_disk = 0
    cingular_inds_in = np.where(cingular_tex_clean == cingular_tex_value)[0]
    insula_inds_in = np.where(insula_tex_clean == insula_tex_value)[0]
    # test that there is no intersection between the two poles
    inter = set(cingular_inds_in).intersection(insula_inds_in)
    if inter:
        raise Exception('Cingular and Insular poles are connected ! You should run the poles textures sanity check')
        return

    tex_poles_clean = np.zeros(mesh.vertex().size())
    tex_poles_clean[cingular_inds_in] = cingular_tex_value
    tex_poles_clean[insula_inds_in] = insula_tex_value

    neigh = aims.SurfaceManip.surfaceNeighbours(mesh)
    print('------------------CutMesh')
    (sub_meshes, labels, sub_indexes) = basicTls.cutMesh(mesh, tex_poles_clean)
    print('labels found in the texture ', labels)
    neo_ind = labels.index(neocortex_tex_value)
    neoCortex_mesh = sub_meshes[neo_ind]
    neoCortex_boundary = basicTls.meshBoundary(sub_meshes[neo_ind])
    if len(neoCortex_boundary)>2:
        print(('problem: more than 2 boundaries in the neoCortex, the pole textures have topological defects'))
        raise Exception('more than 2 boundaries in the neoCortex, the pole textures have topological defects')
        return
    neocortex_indices = sub_indexes[neo_ind]
    ins_ind = labels.index(insula_tex_value)
    insula_mesh = sub_meshes[ins_ind]
#    insula_boundary = basicTls.meshBoundary(sub_meshes[ins_ind])
    insula_indices = sub_indexes[ins_ind]
    cing_ind = labels.index(cingular_tex_value)
    cingular_mesh = sub_meshes[cing_ind]
#    cingular_boundary = basicTls.meshBoundary(sub_meshes[cing_ind])
    cingular_indices = sub_indexes[cing_ind]
    "------------------computing the poles path on the original mesh, always from the boundary of insula pole to the boundary of cingular pole"
    cing_tex_boundary = basicTls.textureBoundary(mesh, cingular_tex_clean, cingular_tex_value, neigh)
    ins_tex_boundary = basicTls.textureBoundary(mesh, insula_tex_clean, insula_tex_value, neigh)
    #poles_path = getShortestPath(mesh, ins_tex_boundary[-1], cing_tex_boundary[-1])

    #------------------- compute Npath other possible poles_path using neighbors of the extremities -------------------------
    #----- reorganize the poles_path and boundaries to distinguish anterior and posterior parts of the boundaries arround the poles_path
#    neoCortex_boundary_on_mesh =  [ins_tex_boundary[-1], cing_tex_boundary[-1]]
#    neoCortex_boundary_on_mesh = boundaryReordering(neoCortex_boundary_on_mesh, poles_path.tolist(), np.array(mesh.vertex()))
    "neoCortex_boundary_on_mesh[0] == insula_boundary"
    " boundary[2] == cingular_boundary"
    "neoCortex_boundary[1] == neocortex_poles_path always from insula to cingular pole"
    "the Npath neighbors of poles_path extremities in ins_tex_boundary and cing_tex_boundary are then :"
    "neoCortex_boundary_on_mesh[0][0:Npath-1]"
    "neoCortex_boundary_on_mesh[2][-1:-Npath]"
    ins_tex_boundary_o =  indsToROI(neocortex_indices, ins_tex_boundary[-1])
    cing_tex_boundary_o = indsToROI(neocortex_indices, cing_tex_boundary[-1])
    poles_path = getShortestPath(neoCortex_mesh, ins_tex_boundary_o, cing_tex_boundary_o)
    neoCortex_boundary_o =  [ins_tex_boundary_o, cing_tex_boundary_o]
    neoCortex_boundary_o = boundaryReordering(neoCortex_boundary_o, poles_path.tolist(), np.array(neoCortex_mesh.vertex()))
    # compute other possible poles_path
    all_poles_path = []
    gp = aims.GeodesicPath(neoCortex_mesh, 3, 0)
    if path_displacement == 'backward' or path_displacement == 'both':
        for p in range(1, Npath+1):
            ind_vert_insula = int(neoCortex_boundary_o[0][-p*step])
            ind_vert_cingular = int(neoCortex_boundary_o[2][p*step])#(neoCortex_boundary_o[2][-p*step-1])
            all_poles_path.append( gp.shortestPath_1_1_ind( ind_vert_insula, ind_vert_cingular).arraydata() )
    all_poles_path.append(poles_path)
    if path_displacement == 'forward' or path_displacement == 'both':
        for p in range(1, Npath+1):
            ind_vert_insula = int(neoCortex_boundary_o[0][p*step])
            ind_vert_cingular = int(neoCortex_boundary_o[2][-p*step])#(neoCortex_boundary_o[2][-p*step-1])
            all_poles_path.append( gp.shortestPath_1_1_ind( ind_vert_insula, ind_vert_cingular).arraydata() )
#    path_meshes = basicTls.meshBoundaryMesh(neoCortex_mesh, all_poles_path)
#    aims.write( path_meshes, '/home/toz/paths_all_sav.gii')

#     ws = aims.Writer()
#     tex_out = aims.TimeTexture_S16()
#     tex_out[0].reserve(mesh.vertex().size())  # pre-allocates memory
#     for i in six.moves.xrange(mesh.vertex().size()):
#         if i in poles_path:
#              tex_out[0].append(1)
#         else:
#              tex_out[0].append(0)
#     tex_out[1].reserve(mesh.vertex().size())  # pre-allocates memory
#     for i in six.moves.xrange(mesh.vertex().size()):
#         if i in cing_tex_boundary[-1]:
#              tex_out[1].append(1)
#         else:
#              tex_out[1].append(0)
#     tex_out[2].reserve(mesh.vertex().size())  # pre-allocates memory
#     for i in six.moves.xrange(mesh.vertex().size()):
#         if i in ins_tex_boundary[-1]:
#              tex_out[2].append(1)
#         else:
#              tex_out[2].append(0)
#     ws.write(tex_out, '/home/toz/poles_link.tex')
    neoCortex_open_mesh = []
    neoCortex_open_boundary = []
    neoCortex_square = []
    tested_poles_path = []
    for path in all_poles_path:
        '''poles_path to neocortex'''
        #neocortex_poles_path = indsToROI(neocortex_indices, path)
        #print('neocortex_poles_path',neocortex_poles_path)
        print('------------------path2Boundary')
        #(neoCortex_open_mesh_p, neoCortex_open_boundary_p) = path2Boundary(neoCortex_mesh,neoCortex_boundary,neocortex_poles_path)
        (neoCortex_open_mesh_p, neoCortex_open_boundary_p) = path2Boundary(neoCortex_mesh,neoCortex_boundary, path.tolist())
        l = len(neoCortex_open_boundary_p[1])
        done = False
        for p in tested_poles_path:
            if len(set(neoCortex_open_boundary_p[1]).intersection(set(p))) == l:
                print('pole_path already processed')
                done = True
        if not done:
            tested_poles_path.append(neoCortex_open_boundary_p[1])
            #print('neoCortex_open_boundary',neoCortex_open_boundary)
            #vert = np.array(neoCortex_open_mesh.vertex())
            print('------------------rectConformalMapping')
            neoCortex_square_p = rectConformalMapping(neoCortex_open_mesh_p, neoCortex_open_boundary_p, length, width, 0)

            neoCortex_open_mesh.append(neoCortex_open_mesh_p)
            neoCortex_open_boundary.append(neoCortex_open_boundary_p)
            neoCortex_square .append(neoCortex_square_p)

    path_meshes = basicTls.meshBoundaryMesh(neoCortex_mesh, tested_poles_path)

    return (neoCortex_square, neoCortex_open_boundary, neocortex_indices, insula_indices, cingular_indices, insula_mesh, cingular_mesh, neoCortex_open_mesh, path_meshes)



####################################################################
#
# HIP
#
####################################################################
def hip(mesh, insula_tex_clean, cingular_tex_clean, length, width):
#    square_ratio = 4.5
    neocortex_tex_value = 0
    insula_tex_value = 180
    cingular_tex_value = 1
    #write_all_steps_to_disk = 0
    cingular_inds_in = np.where(cingular_tex_clean == cingular_tex_value)[0]
    insula_inds_in = np.where(insula_tex_clean == insula_tex_value)[0]
    # test that there is no intersection between the two poles
    inter = set(cingular_inds_in).intersection(insula_inds_in)
    if inter:
        raise Exception('Cingular and Insular poles are connected ! You should run the poles textures sanity check')
        return

    tex_poles_clean = np.zeros(mesh.vertex().size())
    tex_poles_clean[cingular_inds_in] = cingular_tex_value
    tex_poles_clean[insula_inds_in] = insula_tex_value

    neigh = aims.SurfaceManip.surfaceNeighbours(mesh)
    print('------------------CutMesh')
    (sub_meshes, labels, sub_indexes) = basicTls.cutMesh(mesh, tex_poles_clean)
    print('labels found in the texture ', labels)
    neo_ind = labels.index(neocortex_tex_value)
    neoCortex_mesh = sub_meshes[neo_ind]
    neoCortex_boundary = basicTls.meshBoundary(sub_meshes[neo_ind])
    if len(neoCortex_boundary)>2:
        print(('problem: more than 2 boundaries in the neoCortex, the pole textures have topological defects'))
        raise Exception('more than 2 boundaries in the neoCortex, the pole textures have topological defects')
        return
    neocortex_indices = sub_indexes[neo_ind]
    ins_ind = labels.index(insula_tex_value)
    insula_mesh = sub_meshes[ins_ind]
#    insula_boundary = basicTls.meshBoundary(sub_meshes[ins_ind])
    insula_indices = sub_indexes[ins_ind]
    cing_ind = labels.index(cingular_tex_value)
    cingular_mesh = sub_meshes[cing_ind]
#    cingular_boundary = basicTls.meshBoundary(sub_meshes[cing_ind])
    cingular_indices = sub_indexes[cing_ind]
    print('------------------poles path, always from insula to cingular pole')
    cing_tex_boundary = basicTls.textureBoundary(mesh, cingular_tex_clean, cingular_tex_value, neigh)
    ins_tex_boundary = basicTls.textureBoundary(mesh, insula_tex_clean, insula_tex_value, neigh)
    poles_path = getShortestPath(mesh, ins_tex_boundary[-1], cing_tex_boundary[-1])
#     ws = aims.Writer()
#     tex_out = aims.TimeTexture_S16()
#     tex_out[0].reserve(mesh.vertex().size())  # pre-allocates memory
#     for i in six.moves.xrange(mesh.vertex().size()):
#         if i in poles_path:
#              tex_out[0].append(1)
#         else:
#              tex_out[0].append(0)
#     tex_out[1].reserve(mesh.vertex().size())  # pre-allocates memory
#     for i in six.moves.xrange(mesh.vertex().size()):
#         if i in cing_tex_boundary[-1]:
#              tex_out[1].append(1)
#         else:
#              tex_out[1].append(0)
#     tex_out[2].reserve(mesh.vertex().size())  # pre-allocates memory
#     for i in six.moves.xrange(mesh.vertex().size()):
#         if i in ins_tex_boundary[-1]:
#              tex_out[2].append(1)
#         else:
#              tex_out[2].append(0)
#     ws.write(tex_out, '/home/toz/poles_link.tex')

    '''poles_path to neocortex'''
    neocortex_poles_path = indsToROI(neocortex_indices, poles_path)
    #print('neocortex_poles_path',neocortex_poles_path)
    print('------------------path2Boundary')
    (neoCortex_open_mesh, neoCortex_open_boundary) = path2Boundary(neoCortex_mesh,neoCortex_boundary,neocortex_poles_path)
    #print('neoCortex_open_boundary',neoCortex_open_boundary)
    #vert = np.array(neoCortex_open_mesh.vertex())
    print('------------------rectConformalMapping')
    neoCortex_square = rectConformalMapping(neoCortex_open_mesh, neoCortex_open_boundary, length, width, 0)
################################
# return neoCortex_open_mesh instead of neoCortex_mesh
################################
    return (neoCortex_square, neoCortex_open_boundary, neocortex_indices, insula_indices, cingular_indices, insula_mesh, cingular_mesh, neoCortex_open_mesh)
#     if write_all_steps_to_disk:
#         print('------------------textureBoundary')
#         ws.write(neoCortex_mesh, '/home/toz/ammon_Lwhite_neocortex_cut_mesh.mesh')
#         ws.write(insula_mesh, '/home/toz/ammon_Lwhite_insula_cut_mesh.mesh')
#         ws.write(cingular_mesh, '/home/toz/ammon_Lwhite_cingular_cut_mesh.mesh')
#         ws.write(meshBoundaryMesh(mesh, cing_tex_boundary), '/home/toz/ammon_Lwhite_decim_cing_boundary.mesh' )
#         ws.write(meshBoundaryMesh(mesh, ins_tex_boundary), '/home/toz/ammon_Lwhite_decim_ins_boundary.mesh' )
#         print(poles_path)
#         print(neocortex_poles_path)
#         tex_out = aims.TimeTexture_S16()
#         tex_out[0].reserve(neoCortex_mesh.vertex().size())  # pre-allocates memory
#         for i in six.moves.xrange(neoCortex_mesh.vertex().size()):
#             if i in neocortex_poles_path:
#                 tex_out[0].append(1)
#             else:
#                 tex_out[0].append(0)
#         ws.write(tex_out, '/home/toz/ammon_Lwhite_neocortex_poles_link.tex')
#         ws.write(neoCortex_square, '/home/toz/ammon_Lwhite_square.mesh')
#     "open_neocortex_indices = neocortex_indices....................plus path???"

####################################################################
#
# HOP
#
####################################################################
def hop(cstrBalance, neoCortex_square, neoCortex_open_boundary, neoCortex_open_mesh,full_sulci, side, model=None):
    vert = np.array(neoCortex_square.vertex())


    
    '''
    realign the S.C. to 0, should already be the case!
    '''
    try :
        SC_ind = full_sulci.names.index(('S.C._'+side))   
        SC_label = full_sulci.labels[SC_ind]
        print('SC_label: ', SC_label)
#    full_sulci.sulcalLines[SC_ind].printArgs()
        translation = -full_sulci.sulcalLines[SC_ind].barycenter[0]
    except ValueError:
        print('----------------------------------------------------------------')
        print('Central sulcus is missing!')
        translation = 0
        print('no translation applied')
##        if model is not None:
##            translation = model.left
##            print('translation is driven by the boundaries')
##        else:
##            translation = 0
##            print('no translation applied')
        print('----------------------------------------------------------------')

    vert[:, 0] = vert[:, 0] + translation # * np.ones(vert.shape[0])
    neoCortex_square.vertex().assign([aims.Point3df(x) for x in vert])
    full_sulci.updateVertices(vert)

    if model is None:
        model = md.Model()

    full_sulci.sulcalLine2SulcalConstraint(model)
#    full_sulci.sulcalLines[SC_ind].printArgs()
    if model is None:
        model.setBoundary(vert[neoCortex_open_boundary[0][0], 0], vert[neoCortex_open_boundary[0][-1], 0], vert[neoCortex_open_boundary[2][0], 1], vert[neoCortex_open_boundary[0][0], 1])
        model.setAxisCoord(full_sulci)


    Lx, LB = pdeTls.computeMeshLaplacian(neoCortex_open_mesh, lap_type='conformal')#neoCortex_square)#

    neoCortex_square_cstr = cstrRectConformalMapping(Lx, model, neoCortex_square, neoCortex_open_boundary, full_sulci, cstrBalance)
    return neoCortex_square_cstr

####################################################################
#
# compute comformal mapping of a mesh to a disk
#
####################################################################
def mesh2Disk(insula_mesh, insula_boundary, insula_lon):
    insula_bound_rad = np.pi * (insula_lon[insula_boundary] - 180) / 180 
    circle = np.array([np.cos(insula_bound_rad), np.sin(insula_bound_rad)])
    insula_disk = diskConformalMapping(insula_mesh, insula_boundary, circle)
    insula_lon, insula_lat = coordinatesFromDisk(insula_disk)
    return (insula_lon, insula_lat, insula_disk)
    
    
####################################################################
#
# HIP-HOP
#
####################################################################
def hipHop(mesh, insula_tex_clean, cingular_tex_clean, texture_sulci, side, model=None):
#    square_ratio = 4.5
    length = 4.5
    width = 1
    cstrBalance = 200
    poles_lat_insula = 30
    poles_lat_cingular = 30
    neocortex_tex_value = 0
    insula_tex_value = 180
    cingular_tex_value = 1
    write_all_steps_to_disk = 0
    if side == 'right':
        SC_label = 44#25
    else:
        SC_label = 43#25       
    print('max(cingular_tex_clean) : ', np.max(cingular_tex_clean))
    print('max(insula_tex_clean) : ', np.max(insula_tex_clean))
    neigh = aims.SurfaceManip.surfaceNeighbours(mesh)
    #    cingular_tex_clean, cing_tex_boundary = poleTextureClean(mesh, texture_poles, cingular_tex_value)    #    insula_tex_clean, ins_tex_boundary = poleTextureClean(mesh, texture_poles, insula_tex_value)
    tex_poles_clean = np.zeros(cingular_tex_clean.size)
    tex_poles_clean[np.where(cingular_tex_clean == cingular_tex_value)[0]] = cingular_tex_value
    tex_poles_clean[np.where(insula_tex_clean == insula_tex_value)[0]] = insula_tex_value
    print('------------------CutMesh')
    (sub_meshes, labels, sub_indexes) = basicTls.cutMesh(mesh, tex_poles_clean)
    print(labels)
    neo_ind = labels.index(neocortex_tex_value)
    neoCortex_mesh = sub_meshes[neo_ind]
    neoCortex_boundary = basicTls.meshBoundary(sub_meshes[neo_ind])
    neocortex_indices = sub_indexes[neo_ind]
    ins_ind = labels.index(insula_tex_value)
    insula_mesh = sub_meshes[ins_ind]
    insula_boundary = basicTls.meshBoundary(sub_meshes[ins_ind])
    insula_indices = sub_indexes[ins_ind]
    cing_ind = labels.index(cingular_tex_value)
    cingular_mesh = sub_meshes[cing_ind]
    cingular_boundary = basicTls.meshBoundary(sub_meshes[cing_ind])
    cingular_indices = sub_indexes[cing_ind]
    print('------------------poles path, always from insula to cingular pole')
    cing_tex_boundary = basicTls.textureBoundary(mesh, cingular_tex_clean, cingular_tex_value, neigh)
    ins_tex_boundary = basicTls.textureBoundary(mesh, insula_tex_clean, insula_tex_value, neigh)
    poles_path = getShortestPath(mesh, ins_tex_boundary[0], cing_tex_boundary[0])
    "poles_path to neocortex"
    neocortex_poles_path = indsToROI(neocortex_indices, poles_path)
    print('------------------path2Boundary')
    (neoCortex_open_mesh, neoCortex_open_boundary) = path2Boundary(neoCortex_mesh,neoCortex_boundary,neocortex_poles_path)
    vert = np.array(neoCortex_open_mesh.vertex())
    print('------------------rectConformalMapping')
    neoCortex_square = rectConformalMapping(neoCortex_open_mesh, neoCortex_open_boundary, length, width, 0)
    print('------------------solveInvertedPolygon')
    (neoCortex_square, nb_inward_evol) = solveInvertedPolygon(neoCortex_square, neoCortex_open_boundary, 100)
    print(nb_inward_evol)

    if write_all_steps_to_disk:
        print('------------------textureBoundary')
        ws = aims.Writer()
        ws.write(neoCortex_mesh, '/home/toz/ammon_Lwhite_neocortex_cut_mesh.mesh')
        ws.write(insula_mesh, '/home/toz/ammon_Lwhite_insula_cut_mesh.mesh')
        ws.write(cingular_mesh, '/home/toz/ammon_Lwhite_cingular_cut_mesh.mesh')
        ws.write(meshBoundaryMesh(mesh, cing_tex_boundary), '/home/toz/ammon_Lwhite_decim_cing_boundary.mesh' )
        ws.write(meshBoundaryMesh(mesh, ins_tex_boundary), '/home/toz/ammon_Lwhite_decim_ins_boundary.mesh' )
        print(poles_path)
        tex_out = aims.TimeTexture_S16()
        tex_out[0].reserve(mesh.vertex().size())  # pre-allocates memory
        for i in six.moves.xrange(mesh.vertex().size()):
            if i in poles_path:
                tex_out[0].append(1)
            else:
                tex_out[0].append(0)
        ws.write(tex_out, '/home/toz/ammon_Lwhite_decim_poles_link.tex')
        print(neocortex_poles_path)
        tex_out = aims.TimeTexture_S16()
        tex_out[0].reserve(neoCortex_mesh.vertex().size())  # pre-allocates memory
        for i in six.moves.xrange(neoCortex_mesh.vertex().size()):
            if i in neocortex_poles_path:
                tex_out[0].append(1)
            else:
                tex_out[0].append(0)
        ws.write(tex_out, '/home/toz/ammon_Lwhite_neocortex_poles_link.tex')
        ws.write(neoCortex_square, '/home/toz/ammon_Lwhite_square.mesh')
    "open_neocortex_indices = neocortex_indices....................plus path???"
#    texture_sulci[np.where(cingular_tex_clean == cingular_tex_value)[0]]=0
#    texture_sulci[np.where(insula_tex_clean == insula_tex_value)[0]]=0
    tex_cstr_square = texture2ROI(texture_sulci, neocortex_indices)

    full_sulci = slSet.SulcalLinesSet()
    full_sulci.extractFromTexture(texture_sulci, mesh)
    full_sulci.updateIndices(neocortex_indices)
    vert = np.array(neoCortex_square.vertex())
    full_sulci.updateVertices(vert)
    print('SC_label: ', SC_label)
    SC_ind = full_sulci.labels.index(SC_label)
    full_sulci.sulcalLines[SC_ind].printArgs()
    translation = -full_sulci.sulcalLines[SC_ind].barycenter[0]
    vert[:, 0] = vert[:, 0] + translation # * np.ones(vert.shape[0])
    neoCortex_square.vertex().assign([aims.Point3df(x) for x in vert])
    neoCortex_square.updateNormals()
    full_sulci.updateVertices(vert)

    model = md.Model()
    model.printArgs()
    full_sulci.sulcalLine2SulcalConstraint(model)
    full_sulci.printArgs()

    model.setBoundary(vert[neoCortex_open_boundary[0][0], 0], vert[neoCortex_open_boundary[0][-1], 0], vert[neoCortex_open_boundary[2][0], 1], vert[neoCortex_open_boundary[0][0], 1])
    model.printArgs()
    model.setAxisCoord(full_sulci)
    model.printArgs()
    model.saveToFile('/home/toz/model_current.txt')

    Lx, LB = pdeTls.computeMeshLaplacian(neoCortex_square, lap_type='conformal')#neoCortex_open_mesh)

    neoCortex_square_cstr = cstrRectConformalMapping(Lx, model, neoCortex_square, neoCortex_open_boundary, full_sulci, cstrBalance)
    (neoCortex_square_cstr, nb_inward_cstr_evol) = solveInvertedPolygon(neoCortex_square_cstr, neoCortex_open_boundary, 100)
    print('nb_inward_cstr_evol', nb_inward_cstr_evol)

    lon, lat = computeCoordinates(mesh, neocortex_indices, neoCortex_square_cstr, neoCortex_open_boundary, poles_lat_insula, poles_lat_cingular)

    print('param de l insula')
    insula_lon = texture2ROI(lon, insula_indices)
    insula_bound_rad = np.pi * (insula_lon[insula_boundary] - 180) / 180 
    circle = np.array([np.cos(insula_bound_rad), np.sin(insula_bound_rad)])
    insula_disk = diskConformalMapping(insula_mesh, insula_boundary, circle)
    insula_lon, insula_lat = coordinatesFromDisk(insula_disk, poles_lat_insula)
    lon[insula_indices] = insula_lon
    lat[insula_indices] = insula_lat

    cingular_lon = texture2ROI(lon, cingular_indices)
    cingular_bound_rad = np.pi * (cingular_lon[cingular_boundary] - 180) / 180 
    circle = np.array([np.cos(cingular_bound_rad), np.sin(cingular_bound_rad)])
    cingular_disk = diskConformalMapping(cingular_mesh, cingular_boundary, circle)
    cingular_lon, cingular_lat = coordinatesFromDisk(cingular_disk, poles_lat_cingular)
    lon[cingular_indices] = cingular_lon
    lat[cingular_indices] = 180 - cingular_lat
    if write_all_steps_to_disk:
        s_tex_cstr_square = aims.TimeTexture_S16()
        s_tex_cstr_square[0].assign(tex_cstr_square)
        for pt in neocortex_poles_path:
            s_tex_cstr_square[0].append(0)
        ws.write(s_tex_cstr_square, '/home/toz/ammon_Lwhite_neocortex_cstr.tex')
        full_sulci.save('/home/toz/ammon_Lwhite_square_full_sulci')
        model.save('/home/toz/model.mesh')
        ws.write(neoCortex_square_cstr, '/home/toz/ammon_Lwhite_square_cstr_'+str(cstrBalance)+'.mesh')
        ws.write(insula_disk, '/home/toz/ammon_Lwhite_insula_disk.mesh')
        ws.write(cingular_disk, '/home/toz/ammon_Lwhite_cingular_disk.mesh')
    print('hip-hop done!')
    return (lon, lat)

