'''
Created on Jan 15, 2013

@author: toz
'''
import numpy as np
from soma import aims
from soma import aimsalgo
from scipy import sparse


def ismember(ar1, ar2):
    if np.__version__<1.6:
        (uni, inds) = np.unique1d(ar1, False, True)## deprecated since numpy.__version__ > 1.6
        I = np.setmember1d(uni, ar2)
    else:
        (uni, inds) = np.unique(ar1, False, True)
        I = np.in1d(uni, ar2)
    return np.reshape(I[inds], ar1.shape)

####################################################################
#
# compute the 3 angles of each triangle in a mesh
#
####################################################################
def meshPolygonAngles(vertex, polygon):
    angles_out = np.zeros(polygon.shape)
    print angles_out.shape
    for ind,p in enumerate(polygon):
        tet = polygonAngles(vertex[p,:])
        angles_out[ind,:] = tet
    return angles_out

####################################################################
#
# compute the 3 angles between the 3 vertices given in parameter
# for too small or flat triangles, return [0, 0, 0]
#
####################################################################
def polygonAngles(vertices):
    tet = np.zeros(3)
    a = fastNorm(vertices[1,:] - vertices[0,:])
    b = fastNorm(vertices[2,:] - vertices[0,:])
    c = fastNorm(vertices[2,:] - vertices[1,:])
#     print vertices[2,:] - vertices[0,:]
#     print b
    if a>0 and b>0 and c>0: #else: flat triangle
        a2 = a**2
        b2 = b**2
        c2 = c**2
        tet[0] = np.arccos((a2+b2-c2)/(2*a*b))
        tet[1] = np.arccos((a2+c2-b2)/(2*a*c))
        tet[2] = np.arccos((c2+b2-a2)/(2*c*b))
    return tet

####################################################################
#
# compute the norm distortions between the two meshes given following
# the type of distortions is given in parameter
# the two meshes must have the same number of vertex
#
####################################################################
def meshDdistortions(FV1,FV2,type):
     distortions_out = np.zeros(FV1.vertex.shape[0])
#     
#     switch type
#         case 'distance'
#             VertConn = compute_vertex_ring(FV1.faces');
#             for v=1:length(FV1.vertices(:,1))
#                 for n=1:length(VertConn{v})
#                     tex1(v)=tex1(v)+(euclidian_dist(FV2.vertices(v,:),FV2.vertices(VertConn{v}(n),:))-euclidian_dist(FV1.vertices(v,:),FV1.vertices(VertConn{v}(n),:)))^2;
#                 end
#                 tex1(v)=tex1(v)/(2*length(VertConn{v}));%normalisation par rapport au nb de voisins (x2 car chaque arete est parcourue 2 fois)
#             end
#         case 'angle'
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
#         otherwise
#             error('type of distortion not valid')
     return distortions_out

####################################################################
#
# compute the norm of a vector very fast
#
####################################################################
def fastNorm(x):
    return np.sqrt(x.dot(x))

####################################################################
#
# read sulcus label translation file and return a dictionnary
#
####################################################################
def readSulcusLabelTranslationFile(sulcus_label_file):
    sulc_labels = []
    with open(sulcus_label_file,'r') as inf:
        for line in inf:
            sulc_labels.append(line.split())    
    sulc_labels_dict = dict((int(value), key) for (key, value) in sulc_labels)
    return sulc_labels_dict
    
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
    print 'Creating new mesh'
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
def computeMeshWeights(mesh):
    threshold = 0.00001 #np.spacing(1)??
    print '    Computing mesh weights'
    vert = np.array(mesh.vertex())
    poly = np.array(mesh.polygon())

    Nbv = vert.shape[0]
    Nbp = poly.shape[0]
    W = sparse.lil_matrix((Nbv, Nbv))
    #W = zeros((Nv,Nv))
    # this old numpy array representation cannot handle big meshes in memory
    threshold_needed = 0
    for i in range(3):
        i1 = np.mod(i, 3)
        i2 = np.mod(i + 1, 3)
        i3 = np.mod(i + 2, 3)
        print '    ', 3 - i1
        pp = vert[poly[:, i2], :] - vert[poly[:, i1], :]
        qq = vert[poly[:, i3], :] - vert[poly[:, i1], :]
        nopp = np.apply_along_axis(np.linalg.norm, 1, pp)
        noqq = np.apply_along_axis(np.linalg.norm, 1, qq)
        thersh_nopp = np.where(nopp<threshold)[0]
        thersh_noqq = np.where(noqq<threshold)[0]
        if len(thersh_nopp) > 0:
            nopp[thersh_nopp] = threshold
            threshold_needed = 1
        if len(thersh_noqq) > 0:
            noqq[thersh_noqq] = threshold
            threshold_needed = 1
#        print np.min(noqq)
        pp = pp / np.vstack((nopp, np.vstack((nopp, nopp)))).transpose()
        qq = qq / np.vstack((noqq, np.vstack((noqq, noqq)))).transpose()
        ang = np.arccos(np.sum(pp * qq, 1))
        for j in range(Nbp):
#            ind1 = poly[j, i1]
            ind2 = poly[j, i2]
            ind3 = poly[j, i3]
            W[ind2, ind3] = W[ind2, ind3] + 1 / np.tan(ang[j])
            W[ind3, ind2] = W[ind3, ind2] + 1 / np.tan(ang[j])
    if threshold_needed:
        print '    -weight threshold needed-'
    print '    OK'

    li = np.hstack(W.data)
    print 'nb of Nan in weights: ', len(np.where(np.isnan(li))[0])

    return W


####################################################################
#
# compute laplacian of a mesh
#
####################################################################
def computeMeshLaplacian(mesh):
    print '    Computing Laplacian'
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
    print '    OK'

    return L


####################################################################
#
# ensure the texture corresponding to the value tex_val has only one connex component with simple boundary
#
####################################################################
def textureTopologicalCorrection(mesh, atex, tex_val, background_val=0, neigh=None):
    tex_val_indices = np.where(atex == tex_val)[0]
    if not tex_val_indices.size:
        print 'no value ' + str(tex_val) + ' in the input texture!!'
        return list()
    else:
        print str(tex_val_indices.size) + ' vertices have the texture value ' + str(tex_val)
        if neigh is None:
            neigh = aims.SurfaceManip.surfaceNeighbours(mesh)
        labels = np.unique(atex)
        labels = labels.tolist()
        #----------------begin:: ensure a single connex component
        tex2 = aims.TimeTexture_S16()
        tex2[0].reserve(atex.size)
        for i in range(atex.size):
            if i in tex_val_indices:
                tex2[0].push_back(tex_val)
            else:
                tex2[0].push_back(0)
        conn_tex = aimsalgo.AimsMeshLabelConnectedComponent(mesh, tex2, 0, 0)  # aimsalgo.AimsMeshLabelConnectedComponent(mesh, tex2,1,0)
#        aconn_tex = conn_tex[0].arraydata()
        conn_labels = set(conn_tex[0].arraydata())
        conn_labels = conn_labels.difference(set([-1]))  # value -1 corresponds to background
        if len(conn_labels) > 1:
            print len(conn_labels), " connex component(s) for this value in the texture, keeping only the largest component"
            lab_val_indices_nb = [len(np.where(conn_tex[0].arraydata() == lab)[0]) for lab in conn_labels]#.difference(set([-1]))]
            #print lab_val_indices_nb#conn_labels#.difference(set([-1])),lab_val_indices_nb
            index_val_largest_comp = lab_val_indices_nb.index(max(lab_val_indices_nb))
          #  index_val_largest_comp = index_val_largest_comp.tolist()
            adeleted_conn_comp_indices = np.where(conn_tex[0].arraydata() == -1)[0]
            deleted_conn_comp_indices = adeleted_conn_comp_indices.tolist()
            for conn_labels_index in set(range(len(conn_labels))) - set([index_val_largest_comp]):
                    deleted_conn_comp_indices.extend(np.where(conn_tex[0].arraydata() == list(conn_labels)[conn_labels_index])[0])
            atex[deleted_conn_comp_indices] = background_val
#            conn_tex.erase()
#            for i in range(atex.size):
#                if i in deleted_conn_comp_indices:
#                    conn_tex[0].push_back(0)
#                else:
#                    conn_tex[0].push_back(atex[i])
        #----------------end:: ensure a single connex component
        #----------------begin:: fill any hole in this single connex component
        boundary = textureBoundary(mesh, atex, tex_val, neigh)
        if len(boundary) > 1:
            print "filling holes in the largest connex component"
#            ws=aims.Writer()
#            tex_out = aims.TimeTexture_S16()
#            i=0
            while len(boundary) > 1:
                print 'length of boundaries : ',[len(bound) for bound in boundary]
#                i+=1
#                nm= 's12158_Rhippo_cleaned'+str(i)+'.tex'
#                tex_out[0].assign(atex)
#                ws.write(tex_out, nm)
                for l in range(len(boundary) - 1):
                    atex = dilateTexture(atex, neigh, boundary[l])
                boundary = textureBoundary(mesh, atex, tex_val, neigh)
        print 'length of boundaries : ',[len(bound) for bound in boundary]
        #----------------end:: fill any hole in this single connex component
        #----------------begin:: ensure that the boundary do not contain the 3 vertices of a triangle
        print 'cleaning the boundary of the texture'
        atex, boundary = cleanTextureBoundary(mesh, atex, tex_val, boundary[0], background_val, neigh)
    return (atex, boundary)

####################################################################
# ensure that the boundary do not contain the 3 vertices of a triangle
####################################################################
def cleanTextureBoundary(mesh, tex, tex_val, bound, background_val=0, neigh=None):
    if neigh is None:
        neigh = aims.SurfaceManip.surfaceNeighbours(mesh)
    '''check if the 3 vertices of any polygon is on the boundary
    if it is the case, the vertex that can be removed is selected and removed'''    
    poly = np.array(mesh.polygon())
    I = ismember(poly, bound)
    poly_set = poly[I[:, 0] & I[:, 1] & I[:, 2], :]
    print poly_set.shape[0]
    u_bound = set(bound)
    count_list = [bound.count(x) for x in u_bound]
    lu_bound = list(u_bound)
    pb_pt = np.array(lu_bound)[np.where(np.array(count_list) > 1)[0]]
    print pb_pt
    pts_to_remove = []
    for pb_poly in poly_set:
        pts_to_remove.append(pb_poly[np.where(ismember(pb_poly, pb_pt) == False)[0]])
    print pts_to_remove
    if len(pts_to_remove) > 0:
        tex[np.array(pts_to_remove)] = background_val
    boundary = textureBoundary(mesh, tex, tex_val, neigh)
    if len(boundary) > 1:
        print 'complex boundary after cleanTextureBoundary'
    return (tex, boundary)


def dilateTexture(tex, neigh, inds):
    for i in inds:
        tex[np.array(neigh[i].list())] = tex[i]
    return tex


####################################################################
#
# compute indexes that are the boundary of a region defined by value
# in a texture
#
####################################################################
def textureBoundary(mesh, atex, val, neigh=0):
    tex_val_indices = np.where(atex == val)[0]
    if not tex_val_indices.size:
        print 'no value ' + str(val) + ' in the input texture!!'
        return list()
    else:
        ####################################################################
        # print 'the vertices on the boundary have the same texture value (boundary inside the patch)'
        ####################################################################
        if not neigh:
            neigh = aims.SurfaceManip.surfaceNeighbours(mesh)
            
        '''identify the vertices that are on the boundary,         
        i.e that have at least one neigbor that has not the same value in the texture '''
        bound_verts = list()
        for i in tex_val_indices:
            ne_i = np.array(neigh[i].list())
            #print ne_i.size
            #print np.intersect1d_nu(ne_i, tex_val_indices).size
            np_ver = [ int(x) for x in np.__version__.split( '.' ) ]
            if np_ver < [ 1, 6 ]:
                inters_size = np.intersect1d_nu(ne_i, tex_val_indices).size
            else:
                inters_size = np.intersect1d(ne_i, tex_val_indices).size
            if (inters_size != ne_i.size):
                bound_verts.append(i)
                
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
# return list of connected components ORDERED, THE FIRST THE LONGEST
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
                print '[len(b) for b in boundary]',[len(b) for b in boundary]
                boundary.append(bound[first_occ:sec_occ])
                '''remove the loop from current boundary'''
                bound[first_occ:sec_occ] = []
#                print bound
                occurence = listCount(bound)
            boundary[b_ind] = bound
            print '[len(b) for b in boundary]',[len(b) for b in boundary]

            
    "sort the boundaries the first the longest"
    boundaries_len = [len(bound) for bound in boundary]
    inx = range(len(boundaries_len))
    inx.sort(lambda x, y: boundaries_len[x] - boundaries_len[y])
# < = >    inx = np.array(boundaries_len).argsort()
    sort_boundary = [boundary[i] for i in inx]
#    boundary.sort()

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
    bound_mesh = aims.AimsTimeSurface_2()

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
    print 'labels:' + str(labels)
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
    if np.__version__<1.6:
        (uni, inds) = np.unique1d(poly_set, False, True)
    else:
        (uni, inds) = np.unique(poly_set, False, True)        
    submesh = aims.AimsTimeSurface_3()
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



