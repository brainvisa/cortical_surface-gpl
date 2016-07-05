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
from soma import aims
from scipy import sparse
np_ver = [1,6]#[ int(x) for x in np.__version__.split( '.' ) ]

def ismember(ar1, ar2):
    if np_ver < [1, 6]:
        (uni, inds) = np.unique1d(ar1, False, True)## deprecated since numpy.__version__ > 1.6
        I = np.setmember1d(uni, ar2)
    else:
        (uni, inds) = np.unique(ar1, False, True)
        I = np.in1d(uni, ar2)
    return np.reshape(I[inds], ar1.shape)

####################################################################
#
# compute the area of each triangle in a mesh
#
####################################################################
def meshPolygonArea(vert,poly):
    
    pp = vert[poly[:, 1], :] - vert[poly[:, 0], :]
    qq = vert[poly[:, 2], :] - vert[poly[:, 0], :]
    cr  = np.cross(pp,qq)
    area = np.sqrt(np.sum(np.power(cr,2),1))/2

    return area

####################################################################
#
# compute the 3 angles of each triangle in a mesh
#
####################################################################
def meshPolygonAngles(vert, poly):
    angles_out = np.zeros(poly.shape)
    for i in range(3):
        i1 = np.mod(i, 3)
        i2 = np.mod(i + 1, 3)
        i3 = np.mod(i + 2, 3)
        pp = vert[poly[:, i2], :] - vert[poly[:, i1], :]
        qq = vert[poly[:, i3], :] - vert[poly[:, i1], :]      
        noqq = np.sqrt(np.sum(qq * qq, 1))
        nopp = np.sqrt(np.sum(pp * pp, 1))

        pp = pp / np.vstack((nopp, np.vstack((nopp, nopp)))).transpose()
        qq = qq / np.vstack((noqq, np.vstack((noqq, noqq)))).transpose()
        angles_out[:,i] = np.arccos(np.sum(pp * qq, 1))
    return angles_out


####################################################################
#
# compute the distortions between the two meshes given
# the type of distortions is given as parameter
# the two meshes must have the same number of vertices
#
####################################################################
def meshDistortions(mesh1,mesh2,type):

    if type == 'angle':
        print('computing angular distortions')
        angles1 = meshPolygonAngles(np.array(mesh1.vertex()), np.array(mesh1.polygon()))
        angles2 = meshPolygonAngles(np.array(mesh2.vertex()), np.array(mesh2.polygon()))
        diff_angles = angles1 - angles2
        distortions_out = diff_angles#np.sum(np.abs(diff_angles))
#             tet1=mesh_face_angles(FV1);
#             tet2=mesh_face_angles(FV2);       
#             tet_1{1,length(FV1.vertices(:,1))}=[];
#             f_tet=abs(tet1-tet2);
#             for f=1:length(FV1.faces(:,1))
#                 for j=1:3
#                     tet_1{FV1.faces(f,j)}=[tet_1{FV1.faces(f,j)},f_tet(j,f)];
#                 end
#             end
#             for v=1:length(FV1.vertices(:,1))
#                 tex1(v)=sum(tet_1{v}(:))/length(tet_1{v}(:));
#             end
    else:
        print('available types of distortions :: angle')
#         case 'distance'
#             VertConn = compute_vertex_ring(FV1.faces');
#             for v=1:length(FV1.vertices(:,1))
#                 for n=1:length(VertConn{v})
#                     tex1(v)=tex1(v)+(euclidian_dist(FV2.vertices(v,:),FV2.vertices(VertConn{v}(n),:))-euclidian_dist(FV1.vertices(v,:),FV1.vertices(VertConn{v}(n),:)))^2;
#                 end
#                 tex1(v)=tex1(v)/(2*length(VertConn{v}));%normalisation par rapport au nb de voisins (x2 car chaque arete est parcourue 2 fois)
#             end
#         case 'area'
#             tet1=mesh_face_area(FV1);
#             tet2=mesh_face_area(FV2);       
#             tet_1{1,length(FV1.vertices(:,1))}=[];
#             f_tet=abs(tet1-tet2);
#             for f=1:length(FV1.faces(:,1))
#                 for j=1:3
#                     tet_1{FV1.faces(f,j)}=[tet_1{FV1.faces(f,j)},f_tet(f)];
#                 end
#             end
#             for v=1:length(FV1.vertices(:,1))
#                 tex1(v)=sum(tet_1{v}(:))/length(tet_1{v}(:));
#             end

    return distortions_out


    
####################################################################
#
# compute unit normals to the faces
#
####################################################################
def meshPolygonNormal(mesh):
    vert = np.array(mesh.vertex())
    poly = np.array(mesh.polygon())

    normalf = np.cross(vert[poly[:, 1], :] - vert[poly[:, 0], :], vert[poly[:, 2], :]-vert[poly[:, 0], :])
    d = np.sqrt(np.sum(np.power(normalf, 2), 1))
    d[d < 0.00000001] = 1 #problem with eps??
    normalf = normalf / np.transpose(np.tile(d, (3, 1)))
#if nargout>1
#    % unit normal to the vertex
#    normal = zeros(nvert,3);
#    for i=1:nface
#        f = faces(i,:);
#        for j=1:3
#            normal(f(j),:) = normal(f(j),:) + normalf(i,:);
#        end
#    end
#    % normalize
#    d = sqrt( sum(normal.^2,2) ); d(d<eps)=1;
#    normal = normal ./ repmat( d, 1,3 );
#end
    return normalf


####################################################################
#
# compute the adjacency matrix of a mesh
# adja(i,j) = 2 if the edge (i,j) is in two faces
# adja(i,j) = 1 means that i and j are on the boundary of the mesh
# adja(i,j) = 0 elsewhere
#
####################################################################
def meshAdjacencyMatrix(mesh):
    Nv = np.array(mesh.vertex()).shape[0]
    poly = np.array(mesh.polygon())
    adja = sparse.lil_matrix((Nv, Nv))

    for tri in poly:
        tri = np.sort(tri)
        adja[tri[0], tri[1]] = adja[tri[0], tri[1]] + 1
        adja[tri[0], tri[2]] = adja[tri[0], tri[2]] + 1
        adja[tri[1], tri[2]] = adja[tri[1], tri[2]] + 1

    return adja


####################################################################
#
# build the boundary by traversing edges
# return list of connected components ORDERED ACCORDING TO THEIR LENGTH, i.e. THE FIRST THE SHORTEST
# complex boundary corresponds to multiple holes in the surface or bad shaped boundary
####################################################################
def edges2Boundary(li, lj):
    tag = np.zeros(li.size)
    boundary = [[]]
    bound_ind = 0
    curr_edge_i = 0
    boundary[bound_ind].extend([li[curr_edge_i], lj[curr_edge_i]])
    tag[curr_edge_i] = 1

    reverse = 1
    while (np.where(tag == 1)[0].size != tag.size):
#        print boundary
        p = boundary[bound_ind][-1]
        curr_edge_i = np.where((li == p) & (tag == 0))[0]
        if (curr_edge_i.size == 0):
            curr_edge_j = np.where((lj == p) & (tag == 0))[0]
            if (curr_edge_j.size == 0):
                if reverse:
                    boundary[bound_ind].reverse()
                    reverse = 0
                else:
#                    if boundary[bound_ind][0] == boundary[bound_ind][-1]:
#                        boundary[bound_ind].pop()
                    bound_ind += 1
#                    print str(bound_ind+1)+' boundaries in this mesh'
                    reverse = 1
                    new_first = np.where((tag == 0))[0][0]
#                    print new_first
#                    print [li[new_first], lj[new_first]]
                    boundary.append([li[new_first], lj[new_first]])
                    tag[new_first] = 1
            else:
                boundary[bound_ind].append(li[curr_edge_j[0]])
                tag[curr_edge_j[0]] = 1
        else:
            boundary[bound_ind].append(lj[curr_edge_i[0]])
            tag[curr_edge_i[0]] = 1
    "concatenate boundary pieces"
    bound_conn = boundariesIntersection(boundary)
    while len(bound_conn) > 0:  # while np.array(bound_conn)>0:
#        print 'concatenate boundaries'
        cat_bound = catBoundary(boundary[bound_conn[0][0]], boundary[bound_conn[0][1]], bound_conn[0][2])
#        print  cat_bound
        boundary[bound_conn[0][0]] = cat_bound
        boundary.pop(bound_conn[0][1])
        bound_conn = boundariesIntersection(boundary)
        
    for bound in boundary:
        if bound[0] == bound[-1]:
            bound.pop()

    ''' cut complex boundaries
    '''
    for b_ind,bound in enumerate(boundary):
        occurence = listCount(bound)
        if max(occurence.keys()) > 1:
            print 'complex boundary --> cut into simple parts'
            while max(occurence.keys()) > 1:
#                 '''find the vertex that is taken more than one time in the boundary'''
#                print bound
#                print len(bound)
#                print 'occurence',occurence
                ite = occurence[max(occurence.keys())]
                first_occ = bound.index(ite)
                sec_occ = first_occ + 1 + bound[first_occ + 1:].index(ite)
#                print 'first_occ',first_occ
#                print 'sec_occ',sec_occ
                '''create a new boundary that corresponds to the loop '''
#                print '[len(b) for b in boundary]',[len(b) for b in boundary]
                boundary.append(bound[first_occ:sec_occ])
                '''remove the loop from current boundary'''
                bound[first_occ:sec_occ] = []
#                print bound
                occurence = listCount(bound)
            boundary[b_ind] = bound
#            print '[len(b) for b in boundary]',[len(b) for b in boundary]

            
    "sort the boundaries the first the longest"
    boundaries_len = [len(bound) for bound in boundary]
    inx = range(len(boundaries_len))
    inx.sort(lambda x, y: boundaries_len[x] - boundaries_len[y])
# < = >    inx = np.array(boundaries_len).argsort()
    sort_boundary = [boundary[i] for i in inx]
#    boundary.sort()
#    print 'in boundary boundaries_len = ',[len(bound) for bound in sort_boundary]
    return sort_boundary

#    boundary_cat = {}
#    nb_bound = bound_ind
#    curr_bound = 0
#    boundary_cat[curr_bound] = boundary[curr_bound]
#    print boundary
#    tag_bound = np.zeros((nb_bound,1))
#    tag_bound[curr_bound] = 1
#    while ( np.where(tag_bound == 1)[0].size ! =  tag_bound.size):
#        print curr_bound
#        print boundary[curr_bound]
#        for bound_ind in range(curr_bound+1,nb_bound):
#            common = set(boundary_cat[curr_bound]).intersection(set(boundary[bound_ind]))
#            if common:
#                print common
#                print boundary_cat[curr_bound]
#                print boundary[bound_ind]
#                boundary_cat[curr_bound].extend(boundary[bound_ind])
#                tag_bound[bound_ind] = 1
#                print bound_ind
#        print boundary_cat
#        print boundary
#        curr_bound+ = 1
#        print tag_bound
#        if np.where(tag_bound == 1)[0].size ! =  tag_bound.size:
#            print np.where((tag_bound == 0))
#            next = np.where((tag_bound == 0))[0][0]
#            print next
#            boundary_cat[curr_bound] = boundary[next]
#            tag_bound[next] = 1

#    "find complex point"
#    hist = np.bincount(np.hstack((li,lj)))
#    c_pts = np.where(hist>2)[0]
#    tag_c_pts = hist[c_pts]
#    print c_pts
#    print tag_c_pts
#    print np.vstack((li,lj)).transpose()
#    tag = np.zeros(li.size)
#    boundary = {}
#    bound_ind = 0
#    if c_pts[0]:#bad shaped boundary
#        if c_pts[0] in list(li):
#            curr_edge_i = list(li).index(c_pts[0])
#            boundary[bound_ind] = [li[curr_edge_i], lj[curr_edge_i]]
#        else:
#            curr_edge_i = list(lj).index(c_pts[0])
#            boundary[bound_ind] = [lj[curr_edge_i],li[curr_edge_i]]
#    else:
#        curr_edge_i = 0
#        boundary[bound_ind] = [li[curr_edge_i], lj[curr_edge_i]]
#    tag[curr_edge_i] = 1
#    print list(c_pts)
##    reverse = 1
#    c_pt = 0
#    while ( np.where(tag == 1)[0].size ! =  tag.size):
##        print boundary
#        p = boundary[bound_ind][-1]
#        if p in list(c_pts):
#            print p
#            c_pt = 1
#        curr_edge_i = np.where((li == p) & (tag == 0))[0]
#        if (curr_edge_i.size == 0):
#            curr_edge_j = np.where((lj == p) & (tag == 0))[0]
#            if (curr_edge_j.size == 0):
##                if reverse:
##                    boundary[bound_ind].reverse()
##                    reverse = 0
##                else:# multiple boundaries in this mesh
#                    # no complex point remains but multiple holes
#                    if boundary[bound_ind][0] == boundary[bound_ind][-1]:
#                        boundary[bound_ind].pop()
#                    bound_ind+ = 1
##                    print str(bound_ind+1)+' boundaries in this mesh'
#                    reverse = 1
#                    new_first = np.where((tag == 0))[0][0]
##                    print new_first
##                    print [li[new_first], lj[new_first]]
#                    boundary[bound_ind] = [li[new_first], lj[new_first]]
#                    tag[new_first] = 1
#            else:
#                if c_pt:
#                    bound_ind+ = 1
#                    boundary[bound_ind] = [li[curr_edge_j[0]], lj[curr_edge_j[0]]]
#                    c_pt = 0
#                else:
#                    boundary[bound_ind].append(li[curr_edge_j[0]])
#                tag[curr_edge_j[0]] = 1
#        else:
#            if c_pt:
#                bound_ind+ = 1
#                boundary[bound_ind] = [lj[curr_edge_i[0]],li[curr_edge_i[0]]]
#                c_pt = 0
#            else:
#                boundary[bound_ind].append(lj[curr_edge_i[0]])
#            tag[curr_edge_i[0]] = 1



#    "cut complex boundaries"
#    diff_list = lambda l1,l2: [x for x in l1 if l2.count(x)>1]
#    nb_bound = bound_ind
#    bound_ind = 1
#    while bound_ind<nb_bound:
#        u_boundary = list(set(boundary[bound_ind]))
#        u_boundary.sort()
#        dif_len = len(boundary[bound_ind])-len(u_boundary)
#        while dif_len>0:
#            print dif_len
#            pb_pt = diff_list(u_boundary,boundary[bound_ind])[0]
#            print u_boundary
#            print boundary[bound_ind]
#            print pb_pt
#            pb_pt_ind1 = boundary[bound_ind].index(pb_pt)
#            pb_pt_ind2 = boundary[bound_ind].index(pb_pt,pb_pt_ind1+1)-1
#            print pb_pt_ind1
#            print pb_pt_ind2
#            nb_bound = nb_bound+1
#            boundary[nb_bound] = boundary[bound_ind][pb_pt_ind2:len(boundary[bound_ind])];
#            print boundary[nb_bound]
#            for i in range(pb_pt_ind2-pb_pt_ind1):
#                boundary[bound_ind].pop();
#            u_boundary = list(set(boundary[bound_ind]))
#            u_boundary.sort()
#            dif_len = len(boundary[bound_ind])-len(u_boundary)
#        bound_ind = bound_ind+1
#    return boundary
#    boundary = [li[0], lj[0]]
#    tag[0] = 1
#    p = lj[0]
#    while ( np.where(tag == 1)[0].size ! =  tag.size):
#        wi = np.where((li == p) & (tag == 0))[0]
#        if (wi.size == 0):
#            wj = np.where((lj == p) & (tag == 0))[0][0]
#            boundary.append(li[wj])
#            p = li[wj]
#            tag[wj] = 1
#        else:
#            boundary.append(lj[wi[0]])
#            p = lj[wi[0]]
#            tag[wi[0]] = 1
#    boundary.pop()
#    return np.array(boundary)

####################################################################
#
# count the occurrences of all items in a list and return a dictionary
# that is of the form {nb_occurence:list_item}, which is the opposite of
# standard implementation usually found on the web
# -----------------------------------------------------------------
# in python >= 2.7, collections may be used, see example below
# >>> from collections import Counter
# >>> z = ['blue', 'red', 'blue', 'yellow', 'blue', 'red']
# >>> Counter(z)
# Counter({'blue': 3, 'red': 2, 'yellow': 1})
####################################################################
def listCount(l):
    return dict((l.count(it),it) for it in l)


def catBoundary(bound1, bound2, common_pts):
#    print common_pts
#    ind_common_pts = []
#    for x in common_pts:
#        ind_common_pts.append([bound1.index(x),bound2.index(x)])
#    print ind_common_pts
#    print bound1
#    print bound1[ind_common_pts[0][0]]
    for x in bound2:
        bound1.append(x)
    return bound1


def boundariesIntersection(boundary):
    bound_conn = []
    for bound_ind1 in range(len(boundary) - 1):
        for bound_ind2 in range(bound_ind1 + 1, len(boundary)):
            common = set(boundary[bound_ind1]).intersection(set(boundary[bound_ind2]))
            if common:
                bound_conn.append([bound_ind1, bound_ind2, list(common)])
    return bound_conn


def edges2SimpleBoundary(li, lj):
    tag = np.zeros(li.size)
    boundary = {}
    bound_ind = 0
    curr_edge_i = 0
    boundary[bound_ind] = [li[curr_edge_i], lj[curr_edge_i]]
    tag[curr_edge_i] = 1

    reverse = 1
    while (np.where(tag == 1)[0].size != tag.size):
#        print boundary
        p = boundary[bound_ind][-1]
        curr_edge_i = np.where((li == p) & (tag == 0))[0]
        if (curr_edge_i.size == 0):
            curr_edge_j = np.where((lj == p) & (tag == 0))[0]
            if (curr_edge_j.size == 0):
                if reverse:
                    boundary[bound_ind].reverse()
                    reverse = 0
                else:  # multiple boundaries in this mesh
                    if boundary[bound_ind][0] == boundary[bound_ind][-1]:
                        boundary[bound_ind].pop()
                    bound_ind += 1
#                    print str(bound_ind+1)+' boundaries in this mesh'
                    reverse = 1
                    new_first = np.where((tag == 0))[0][0]
#                    print new_first
#                    print [li[new_first], lj[new_first]]
                    boundary[bound_ind] = [li[new_first], lj[new_first]]
                    tag[new_first] = 1
            else:
                boundary[bound_ind].append(li[curr_edge_j[0]])
                tag[curr_edge_j[0]] = 1
        else:
            boundary[bound_ind].append(lj[curr_edge_i[0]])
            tag[curr_edge_i[0]] = 1
    return boundary


def meshBoundaryMesh(mesh, boundary):
    vert = np.array(mesh.vertex())
    bound_mesh = aims.AimsTimeSurface_2_VOID()

    for bound_ind in range(len(boundary)):
        vv = aims.vector_POINT3DF()
        ee = aims.vector_AimsVector_U32_2()
        for i in boundary[bound_ind]:
            vv.append([vert[i, 0], vert[i, 1], vert[i, 2]])
        for i in range(len(boundary[bound_ind]) - 1):
            ee.append(np.array([i, i + 1], np.uint32))
        bound_mesh.vertex(bound_ind).assign(vv)
        bound_mesh.polygon(bound_ind).assign(ee)
#        pol = np.vstack( (np.zeros( len(boundary[bound_ind])-2, dtype = np.int32 ),boundary[bound_ind][0:-2],boundary[bound_ind][1:-1] )).transpose()
#        print pol.dtype#astype(np.float32).dtype
#        bound_mesh.polygon(bound_ind).assign([ aims.AimsVector(x,'U32') for x in pol ])
    bound_mesh.updateNormals()
    return bound_mesh


####################################################################
#
# compute borders of a mesh
#
####################################################################
def meshBoundary(mesh):
    adja = meshAdjacencyMatrix(mesh)
    r = sparse.extract.find(adja)
    li = r[0][np.where(r[2] == 1)]
    lj = r[1][np.where(r[2] == 1)]

    if (li.size == 0):
        print 'No holes in the surface !!!!'
        return np.array()
    else:
        return edges2Boundary(li, lj)

####################################################################
#
# compute indices that are the boundary of a region defined by value
# in a texture, without any topological or ordering constraint
#
####################################################################
def textureSimpleBoundary(mesh, atex, val, neigh=None):
    tex_val_indices = np.where(atex == val)[0]
    if not tex_val_indices.size:
        print 'no value ' + str(val) + ' in the input texture!!'
        return list()
    else:
        ####################################################################
        # print 'the vertices on the boundary have the same texture value (boundary inside the patch)'
        ####################################################################
        if neigh is None:
            neigh = aims.SurfaceManip.surfaceNeighbours(mesh)

        '''identify the vertices that are on the boundary,
        i.e that have at least one neigbor that has not the same value in the texture '''
        bound_verts = list()
        for i in tex_val_indices:
            ne_i = np.array(neigh[i].list())
            #print ne_i.size
            #print np.intersect1d_nu(ne_i, tex_val_indices).size
            if np_ver < [ 1, 6 ]:
                inters_size = np.intersect1d_nu(ne_i, tex_val_indices).size
            else:
                inters_size = np.intersect1d(ne_i, tex_val_indices).size
            if (inters_size != ne_i.size):
                bound_verts.append(i)
        return bound_verts

####################################################################
#
# compute indexes that are the boundary of a region defined by value
# in a texture
#
####################################################################
def textureBoundary(mesh, atex, val, neigh=None):
    tex_val_indices = np.where(atex == val)[0]
    if not tex_val_indices.size:
        print 'no value ' + str(val) + ' in the input texture!!'
        return list()
    else:
        bound_verts =  textureSimpleBoundary(mesh, atex, val, neigh)
        ''' select the edges that are on the boundary in the polygons
        '''
        adja = meshAdjacencyMatrix(mesh)
        r = sparse.extract.find(adja)
#        print bound_verts
#        print bound_verts[0] in r[0]
#        print r
#        inter = set(r[0]).intersection(set(bound_verts))
#        print inter
#        print set(r[0])-inter
#        print np.intersect1d(np.intersect1d(r[0],bound_verts),np.intersect1d(r[1],bound_verts))
#        print r[0] == bound_verts[0]
#        print bound_verts[0:1]
#        print np.intersect1d(r[0],bound_verts[0:1]).shape
#        print np.intersect1d(np.intersect1d(r[0],bound_verts),np.intersect1d(r[1],bound_verts))    
        inr0 = []
        inr1 = []
        for v in bound_verts:
            inr0.extend(np.where(r[0] == v)[0])
            inr1.extend(np.where(r[1] == v)[0])
        r[2][inr0] = r[2][inr0] + 1
        r[2][inr1] = r[2][inr1] + 1
        li = r[0][np.where(r[2] == 4)]
        lj = r[1][np.where(r[2] == 4)]
#        print 'li lj = '
#        print li
#        print lj
        return edges2Boundary(li, lj)

####################################################################
#
# return the indices of vertex corresponding to the list of polgon given
#
####################################################################
def meshPolygon2vertex(mesh, poly_indices):
    poly = np.array(mesh.polygon())
    return np.unique(poly[poly_indices, :])

####################################################################
#
# cut a hole in a mesh at nodes defined by value in texture
# returns two meshes of hole and mesh-hole
# the hole border belongs to both meshes
#
####################################################################
def cutMesh(mesh, atex):
    atex2 = atex.copy()
    labels = np.around(np.unique(atex))
    labels = labels.tolist()
    labels.reverse()
    sub_meshes = list()
    sub_indexes = list()
    last_label = labels[-1]
    for label_ind in range(len(labels) - 1):
        (sub_mesh, sub_index) = subCutMesh(mesh, atex, labels[label_ind])
        sub_meshes.append(sub_mesh)
        sub_indexes.append(sub_index.tolist())
        boundary = textureBoundary(mesh, atex, labels[label_ind])
        atex2[boundary] = last_label
    (sub_mesh, sub_index) = subCutMesh(mesh, atex2, last_label)
    sub_meshes.append(sub_mesh)
    sub_indexes.append(sub_index.tolist())
    return (sub_meshes, labels, sub_indexes)


def subCutMesh(mesh, atex, val):
    # [FV_cut,vert_ind_cut,boundary_cut,FV_comp,corresp_ind_comp,boundary_comp]=cut_mesh(FV,texture,FIG)
    poly = np.array(mesh.polygon())
    vert = np.array(mesh.vertex())
    tex_val_indices = np.where(atex == val)[0]
    I = ismember(poly, tex_val_indices)
    poly_set = poly[I[:, 0] & I[:, 1] & I[:, 2], :]
#    print tex_val_indices
    if np_ver < [1, 6]:
        (uni, inds) = np.unique1d(poly_set, False, True)
    else:
        (uni, inds) = np.unique(poly_set, False, True)
    submesh = aims.AimsTimeSurface_3_VOID()
    vv = aims.vector_POINT3DF()
    for i in vert[uni, :]:
        vv.append(i)
    submesh.vertex().assign(vv)
    pp = aims.vector_AimsVector_U32_3()
    for i in np.reshape(inds, poly_set.shape):
        pp.append(i)
    submesh.polygon().assign(pp)
    submesh.updateNormals()
    return (submesh, tex_val_indices)
#    [unqVertIds, ~, newVertIndices] = unique(setF);
#    FV_cut.faces = reshape(newVertIndices,size(setF));
#    FV_cut.vertices = FV.vertices(unqVertIds,:);
#
#    boundary_cut = compute_boundary(FV_cut.faces');
#
#    v_bound=vert_ind_cut(boundary_cut{1});
#    texture(v_bound)=0;
#    vert_ind_comp=find(~texture);
#    I=(ismember(FV.faces,vert_ind_comp));
#    setF = FV.faces(I(:,1)&I(:,2)&I(:,3),:);
#    on supprime les faces dont les 3 sommets sont sur la frontiere
#    J=(ismember(setF,v_bound));
#    if ~isempty(find(J(:,1)&J(:,2)&J(:,3)))
#    disp('pb cut_mesh')
#    setF(J(:,1)&J(:,2)&J(:,3),:)=[];
#---------------------------------------------------------
#    [unqVertIds, corresp, newVertIndices] = unique(setF);
#    corresp_ind_comp=setF(corresp);
#
#    FV_comp.faces = reshape(newVertIndices,size(setF));
#    FV_comp.vertices = FV.vertices(unqVertIds,:);
#    boundary_comp=compute_boundary(FV_comp.faces');

####################################################################
#
# compute the edges  between a list of vertices
# the indices in the resulting edges are within [0 length of the list]
#
####################################################################

def vertsIndicesToEdges(mesh, indices, neigh=None):
    if neigh is None:
        neigh = aims.SurfaceManip.surfaceNeighbours(mesh)
    "build the segments that link the vertices"
    Nv = len(indices)
    segm = []
    C = sparse.lil_matrix((Nv, Nv))
    for v in indices:
        ne_i = np.array(neigh[v].list())
        if np_ver < [ 1, 6 ]:
            intersect = np.intersect1d_nu(ne_i, indices)
        else:
            intersect = np.intersect1d(ne_i, indices)
        if intersect is not None:
            v_index = indices.index(v)
            for inter in intersect:
                inter_index = indices.index(inter)
                if C[v_index, inter_index] == 0:
                    segm.append([v_index, inter_index])
                    C[v_index, inter_index] = 1
                    C[inter_index, v_index] = 1
    return segm


####################################################################
#
# compute iso-parameter lines on the mesh
#
####################################################################

def meshIsoLine(mesh, tex, vals):
    #print 'Looking for isoLine'
    points=np.array(mesh.vertex())
    values=np.array(tex[0])
    #print('isoLine: points:', points.shape)
    #print('isoLine: values:', values.shape)
    line=aims.AimsTimeSurface_2_VOID()
    isoV=aims.vector_POINT3DF()
    isoP=aims.vector_AimsVector_U32_2()
    triangles=np.array(mesh.polygon())

    for val in vals:
         sign=np.zeros(values.size)
         sign[values < val]=10
         sign[values >= val]=20

         for tr in triangles:
              count=sign[tr].sum()
              if (count==40):
                   if (sign[tr[0]]==20):
                        v1=interpolateVertices(tr[0], tr[1], points, values, val)
                        v2=interpolateVertices(tr[0], tr[2], points, values, val)
                   elif (sign[tr[1]]==20):
                        v1=interpolateVertices(tr[1], tr[0], points, values, val)
                        v2=interpolateVertices(tr[1], tr[2], points, values, val)
                   elif (sign[tr[2]]==20):
                        v1=interpolateVertices(tr[2], tr[0], points, values, val)
                        v2=interpolateVertices(tr[2], tr[1], points, values, val)
                   isoV, isoP = addSegment(v1, v2, isoV, isoP)
              elif (count==50):
                   if (sign[tr[0]]==10):
                        v1=interpolateVertices(tr[0], tr[1], points, values, val)
                        v2=interpolateVertices(tr[0], tr[2], points, values, val)
                   elif (sign[tr[1]]==10):
                        v1=interpolateVertices(tr[1], tr[0], points, values, val)
                        v2=interpolateVertices(tr[1], tr[2], points, values, val)
                   elif (sign[tr[2]]==10):
                        v1=interpolateVertices(tr[2], tr[0], points, values, val)
                        v2=interpolateVertices(tr[2], tr[1], points, values, val)
                   isoV, isoP=addSegment(v1, v2, isoV, isoP)

    line.vertex().assign(isoV)
    line.polygon().assign(isoP)
    return line

####################################################################
# 
# THIS ONE DOES NOT INTERPOLATE ON EDGES: 
# it sends back a list of vertex indices closest to the iso-line
####################################################################  

def meshAlmostIsoLine(mesh, tex, val):
     #print 'Looking for isoLine'
     points=np.array(mesh.vertex())
     values=np.array(tex[0])
     #print('isoLine: points:', points.shape)
     #print('isoLine: values:', values.shape)
     sign=np.zeros(values.size)
     sign[values < val]=10
     sign[values >= val]=20
     
     si=set()
     
     triangles=np.array(mesh.polygon())
     for tr in triangles:
          count=sign[tr].sum()
          if (count==40):
               if (sign[tr[0]]==20):
                    i1=bestVertex(tr[0], tr[1], values, val)
                    i2=bestVertex(tr[0], tr[2], values, val)
               elif (sign[tr[1]]==20):
                    i1=bestVertex(tr[1], tr[0], values, val)
                    i2=bestVertex(tr[1], tr[2], values, val)
               elif (sign[tr[2]]==20):
                    i1=bestVertex(tr[2], tr[0], values, val)
                    i2=bestVertex(tr[2], tr[1], values, val)
               si.add(i1)
               si.add(i2)
          elif (count==50):
               if (sign[tr[0]]==10):
                    i1=bestVertex(tr[0], tr[1], values, val)
                    i2=bestVertex(tr[0], tr[2], values, val)
               elif (sign[tr[1]]==10):
                    i1=bestVertex(tr[1], tr[0], values, val)
                    i2=bestVertex(tr[1], tr[2], values, val)
               elif (sign[tr[2]]==10):
                    i1=bestVertex(tr[2], tr[0], values, val)
                    i2=bestVertex(tr[2], tr[1], values, val)
               si.add(i1)
               si.add(i2)
     return np.array(list(si))
               
def interpolateVertices(v1, v2, points, texture, valeur):
     t=np.abs(valeur - texture[v1])/np.abs(texture[v1] - texture[v2])
     new=(1-t)*points[v1] + t*points[v2]
     return new
     
def bestVertex(v1, v2, texture, valeur):
     t=np.abs(valeur - texture[v1])/np.abs(texture[v1] - texture[v2])
     if t<0.5:
          return v1
     else:
          return v2
     
def addSegment(v1, v2, vert, seg):
     nv=vert.size()
     a1=0
     a2=0
     
     for i in range(nv):
          v=vert[i]
          d1=v-v1
          d2=v-v2
          if (np.sqrt(d1[0]*d1[0] + d1[1]*d1[1] + d1[2]*d1[2]))<0.0001:
               i1=i
               a1=1
          if (np.sqrt(d2[0]*d2[0] + d2[1]*d2[1] + d2[2]*d2[2]))<0.0001:
               i2=i
               a2=1
 
     if (a1==0):
          vert.append(v1)
          i1=vert.size()-1
     if (a2==0):
          vert.append(v2)
          i2=vert.size()-1
     seg.append([i1, i2])
     return vert, seg 


####################################################################
#
# create cylindric meshes from a list of vertices
#
####################################################################

def linesToTubes(lines, diam=0.25):
     gen = aims.SurfaceGenerator()
     tubes = aims.AimsSurfaceTriangle()

     seg = np.array(lines.polygon())
     vert = np.array(lines.vertex())
     for s in seg:
          tubes += gen.cylinder(vert[np.uint32(s[0])], vert[np.uint32(s[1])], diam, diam, 10, 0)
     return tubes
     
####################################################################
# 
# find the vertex of a mesh closest to a point p
#
####################################################################

def closerVert(p, mesh):
     diff=p-mesh
     return(dot(diff,diff.transpose()).diagonal().argmin())

#################################################################
