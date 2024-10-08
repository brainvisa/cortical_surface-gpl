'''
Created on 2 august 2012

@author: toz
'''
import numpy as np
from soma import aims
from scipy import sparse
from brainvisa.cortical_surface.surface_tools.basic_tools import vertsIndicesToEdges
from six.moves import range

np_ver = [1,6]#[ int(x) for x in np.__version__.split( '.' ) ]


class SulcalLine(object):
    '''
    classdocs
    '''

    def __init__(self, label=None, indices=None, vertices=None, segm=None, sulc_labels_dict=None, name=None):
        '''
        Constructor
        '''
        if label is None:
            self.label = []
        else:
            self.label = label
        if indices is None:
            self.indices = np.array([], np.uint32)
        else:
            self.indices = indices
        if vertices is None:
            self.vertices = np.array([])
        else:
            self.vertices = vertices
        if segm is None:
            self.segm = np.array([], np.uint32)
        else:
            self.segm = segm
        self.color = label
        self.nbVertices = self.vertices.shape[0]
        if sulc_labels_dict is None:
            self.name = 'unknown'
        else:
            self.name = sulc_labels_dict[self.label]
        if name is not None:
            self.name = name
        self.computeAttributes()

    def printArgs(self):
        print('label = ', self.label)
        print('name = ', self.name)
        print('indices = ', self.indices)
        print('vertices = ', self.vertices)
        print('nbVertices = ', self.nbVertices)
        print('segm = ', self.segm)
        print('barycenter = ', self.barycenter)
        print('length = ', self.length)
        print('color = ', self.color)

    def cat(self, sl):
        self.segm = np.vstack((self.segm, sl.segm+len(self.vertices)))
        self.vertices = np.vstack((self.vertices, sl.vertices))
        self.nbVertices = self.vertices.shape[0]
#        self.indices
        self.computeAttributes()

    def updateName(self, sulc_labels_dict):
        self.name = sulc_labels_dict[self.label]

    def updateVertices(self, vertices=None):
        if vertices is None:
            self.vertices = np.array([])
        else:
            self.vertices = vertices
            self.computeAttributes()

    def updateSegm(self, segm=None):
        if segm is None:
            self.segm = np.array([], np.uint32)
        else:
            self.segm = segm

#     def updateIndices(self, indices=None):
#         if indices is None:
#             self.indices = np.array([], np.uint32)
#         else:
#             # remove the inds that are = -1 from indices and self.segm
#             deleted_inds = np.where(indices == -1)[0]
#             nb_deleted_inds = deleted_inds.shape[0]
#             if nb_deleted_inds == 0:
#                 self.indices = np.array(indices, np.uint32)
#             else:
# #                print('sulcalLine indices : ',indices)
#                 # remove the indices that are = -1...
# #                print('deleted_inds : ', deleted_inds)
#                 indices = np.delete(indices, deleted_inds)
#                 self.indices = np.array(indices, np.uint32)
#                 print(indices)
#                 # ...and the corresponding segment in self.segm
#                 m_segm = np.max(self.segm)+1
#                 print('m_segm : ', m_segm)
#                 deleted_segm_tmp = np.array([np.where(self.segm == ind)[0] for ind in deleted_inds])
#                 deleted_segm = np.unique(deleted_segm_tmp.flatten())
# #                print('deleted_segm :', deleted_segm)
# #                print('segm :', self.segm)
#                 self.segm = np.delete(self.segm, deleted_segm, 0)
#                 print('segm :', self.segm)
#                 for el in range(m_segm - nb_deleted_inds, m_segm):
#                     self.segm[np.where(self.segm == el)] = el - nb_deleted_inds
#                 print('segm :', self.segm)

    def computeAttributes(self):
        if self.vertices.shape[0] > 0:
            self.barycenter = np.mean(self.vertices, 0)
#            for s in self.segm:
#                print(s)
#                print(s[0])
#            print(self.vertices.shape)
#            print([self.vertices[s[0], :] for s in self.segm])
            vert_segm = np.array([self.vertices[s[1], :] - self.vertices[s[0], :] for s in self.segm])
            # -np.array(self.vertices)[self.segm[:,0],:]
            nn = np.apply_along_axis(np.linalg.norm, 1, vert_segm)
            self.length = np.sum(nn)
        else:
            self.barycenter = np.array([])
            self.length = np.array([])

    def computePCA(self):
        if self.vertices.shape[0] > 0:
            #compute the PCA
            #----center the data
            data = self.vertices - self.vertices.mean(axis=0)
            #----compute the SVD
            U, s, Vt = np.linalg.svd(data[:,0:2], full_matrices=False)
            self.PCA = Vt
        else:
            self.PCA = np.array([])


    def __iter__(self):
        return self

    def next(self):
        return self

    def extractFromTexture(self, label, tex, mesh, sulc_labels_dict=None, neigh=None):
        atex = np.around(tex)
        tex_val_indices = list(np.where(atex == label)[0])
        Nv = len(tex_val_indices)
        if Nv == 0:
            print('no value ' + str(label) + ' in the input texture, return empty sulcalLine!!')
        else:
            segm = vertsIndicesToEdges(mesh, tex_val_indices, neigh)
            # if neigh is None:
            #     neigh = aims.SurfaceManip.surfaceNeighbours(mesh)
            # "build the segments that link the vertices"
            # segm = []
            # C = sparse.lil_matrix((Nv, Nv))
            # for v in tex_val_indices:
            #     ne_i = np.array(neigh[v].list())
            #     if np_ver < [ 1, 6 ]:
            #         intersect = np.intersect1d_nu(ne_i, tex_val_indices)
            #     else:
            #         intersect = np.intersect1d(ne_i, tex_val_indices)
            #     if intersect is not None:
            #         v_index = tex_val_indices.index(v)
            #         for inter in intersect:
            #             inter_index = tex_val_indices.index(inter)
            #             if C[v_index, inter_index] == 0:
            #                 segm.append([v_index, inter_index])
            #                 C[v_index, inter_index] = 1
            #                 C[inter_index, v_index] = 1
            verts = np.array(mesh.vertex())
            return self.__init__(label, np.array(tex_val_indices, np.uint32), verts[tex_val_indices, :], np.array(segm, np.uint32), sulc_labels_dict)
            
    def toMesh(self):
        out_mesh = aims.AimsTimeSurface_2_VOID()
        verts = aims.vector_POINT3DF()
        poly = aims.vector_AimsVector_U32_2()
        for v in self.vertices:
            verts.append(v)
        for s in self.segm:
            poly.append(np.array(s, np.uint32))
        out_mesh.vertex().assign(verts)
        out_mesh.polygon().assign(poly)
#        out_mesh.vertex().assign(self.verts)
#        out_mesh.polygon().assign(self.segm)
        out_mesh.updateNormals()
        return out_mesh

    def toTex(self):
        out_tex = aims.TimeTexture_S16()
        for i in np.ones(self.vertices.shape[0], np.int16) * self.label:
                out_tex[0].append(int(i))
        return out_tex

    def clean(self):
        "clean the line using geodesic shortest path"
        print('TO DO : clean the line using geodesic shortest path!!!!')
        return self

    def plot(self, plt, color=None, modele=None):
        if modele is None:
            plot_target = False
        else:
            plot_target = True
        if color is None:
            color = ['b']
        for s in self.segm:
            plt.plot(self.vertices[s, 0], self.vertices[s, 1], color, linewidth=3.0)
        if plot_target & isinstance(self, SulcalConstraint):
            self.plotTarget(plt, modele)


class SulcalConstraint(SulcalLine):
    def __init__(self, sl=None, axisID=None, isLon=None, isLat=None, weight=None, vertexWeight=None):
        '''
        Constructor
        '''
        if sl is None:
            super(SulcalConstraint, self).__init__()
        else:
            super(SulcalConstraint, self).__init__(label=sl.label, indices=sl.indices, vertices=sl.vertices, segm=sl.segm,name=sl.name)
        if isLat is None:
            self.isLat = False
        else:
            self.isLat = isLat
        if isLon is None:
            self.isLon = False
        else:
            self.isLon = isLon
        if axisID is None:
            self.axisID = None
        else:
            self.axisID = axisID

        ######################################################################
        if weight is None:
            self.weight = 1#self.nbVertices #self.length  # sum(self.vertexWeight)
        else:
            self.weight = weight
        if vertexWeight is None:
            self.vertexWeight = np.ones(self.nbVertices)# should be used for weighting between vertices inside a sulcalConstraint
        else:
            self.vertexWeight = vertexWeight
        ######################################################################

    def plotTarget(self, plt, modele):
        if self.isLat:
            coord = modele.latitudeAxisCoord[modele.latitudeAxisID.index(self.axisID)]
            for v in range(self.nbVertices):
                plt.plot([self.vertices[v, 0], self.vertices[v, 0]], [self.vertices[v, 1],coord],'g')
        if self.isLon:
            coord = modele.longitudeAxisCoord[modele.longitudeAxisID.index(self.axisID)]
            for v in range(self.nbVertices):
                plt.plot([self.vertices[v, 0], coord], [self.vertices[v, 1],self.vertices[v, 1]],'r')

    def printArgs(self):
        print('isLat = ', self.isLat)
        print('isLon = ', self.isLon)
        print('axisID = ', self.axisID)
        print('vertexWeight = ', self.vertexWeight)
        print('weight = ', self.weight)
        super(SulcalConstraint, self).printArgs()

    def cat(self, sc):
        super(SulcalConstraint, self).cat(sc)
    #     self.setVertexWeight()  # should be used for weighting between vertices inside a sulcaConstraint
    #     ######################################################################
    #     self.weight = 1#self.nbVertices#self.length  # sum(self.verticesWeight)
    #     ######################################################################
    #
    # def setVertexWeight(self, input_weights=None):
    #     if input_weights is None:
    #         self.vertexWeight = np.ones(self.nbVertices)
    #     else:
    #         self.vertexWeight = input_weights
