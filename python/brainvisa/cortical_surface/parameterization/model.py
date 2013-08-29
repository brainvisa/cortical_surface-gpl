'''
Created on Nov 16, 2012

@author: toz
'''
import numpy as np
from soma import aims


class Model(object):
    '''
    classdocs
    '''

    def __init__(self, sulci=None, left=None, right=None, top=None, bottom=None):
        '''
        Constructor
        '''

        self.longitudeAxis = [16, 40, 61, 281, 297, 339, 360]
        self.latitudeAxis = [1, 55, 56, 81, 92, 106]
        if sulci is None:
            self.longitudeAxisCoord = []
            self.latitudeAxisCoord = []
        else:
            self.setAxisCoord(sulci)
        if left == None:
            self.left = []
        else:
            self.left = left
        if right == None:
            self.right = []
        else:
            self.right = right
        if top == None:
            self.top = []
        else:
            self.top = top
        if bottom == None:
            self.bottom = []
        else:
            self.bottom = bottom

    def printArgs(self):
        print 'model ::'
        print '    left = ', self.left
        print '    right = ', self.right
        print '    top = ', self.top
        print '    bottom = ', self.bottom
        print '    longitudeAxis = ', self.longitudeAxis
        print '    longitudeAxisCoord = ', self.longitudeAxisCoord
        print '    latitudeAxis = ', self.latitudeAxis
        print '    latitudeAxisCoord = ', self.latitudeAxisCoord

    def setBoundary(self, left, right, top, bottom):
        self.left = left
        self.right = right
        self.top = top
        self.bottom = bottom

    def saveToMesh(self, fileName):
        ws = aims.Writer()
        ws.write(self.toMesh(), fileName)

    def saveToFile(self, fileName):
        f = open(fileName, 'w')
        txt = 'model ::\n'
#         for col in data[index,:]:
#             txt = txt+' '+str(col)
        txt = txt + 'left = ' + str(self.left) + '\n'
        txt = txt + 'right = ' + str(self.right) + '\n'
        txt = txt + 'top = ' + str(self.top) + '\n'
        txt = txt + 'bottom = ' + str(self.bottom) + '\n'
        txt = txt + 'longitudeAxis = ' + str(self.longitudeAxis) + '\n'
        txt = txt + 'longitudeAxisCoord = ' + str(self.longitudeAxisCoord) + '\n'
        txt = txt + 'latitudeAxis = ' + str(self.latitudeAxis) + '\n'
        txt = txt + 'latitudeAxisCoord = ' + str(self.latitudeAxisCoord) + '\n'
        f.write(txt)
        f.close()

    def toMesh(self, z_coord=None):
        if z_coord is None:
            z_coord = 0
        out_mesh = aims.AimsTimeSurface_2()
        verts = aims.vector_POINT3DF()
        poly = aims.vector_AimsVector_U32_2()
        #-------- boundaries
        verts.append([self.left, self.top, z_coord])
        verts.append([self.right, self.top, z_coord])
        verts.append([self.left, self.bottom, z_coord])
        verts.append([self.right, self.bottom, z_coord])
        poly.append(np.array([0, 1], np.uint32))
        poly.append(np.array([0, 2], np.uint32))
        poly.append(np.array([1, 3], np.uint32))
        poly.append(np.array([2, 3], np.uint32))
        nb_verts = 4
        #-------- axes
        for ax in self.longitudeAxisCoord:
            if ax is not None:
                verts.append([ax, self.top, z_coord])
                verts.append([ax, self.bottom, z_coord])
                poly.append(np.array([nb_verts, nb_verts + 1], np.uint32))
                nb_verts = nb_verts + 2
        for ax in self.latitudeAxisCoord:
            if ax is not None:
                verts.append([self.left, ax, z_coord])
                verts.append([self.right, ax, z_coord])
                poly.append(np.array([nb_verts, nb_verts + 1], np.uint32))
                nb_verts = nb_verts + 2
        out_mesh.vertex().assign(verts)
        out_mesh.polygon().assign(poly)
        out_mesh.updateNormals()
        return out_mesh

    def plot(self, plt):
        for ax in self.longitudeAxisCoord:
            if ax is not None:
                plt.plot([ax, ax], [self.top, self.bottom])
        for ax in self.latitudeAxisCoord:
            if ax is not None:
                plt.plot([self.left, self.right], [ax, ax])
        plt.plot([self.left, self.right], [self.top, self.top], 'r')
        plt.plot([self.right, self.right], [self.top, self.bottom], 'g')
        plt.plot([self.left, self.right], [self.bottom, self.bottom], 'b')
        plt.plot([self.left, self.left], [self.top, self.bottom], 'm')

##########################################################################
# replaced by sulcus2Axis, should not be used anymore
##########################################################################
    def label2Axis(self, label):
        isLat = False
        isLon = False
        axisID = None
        #  longitude labels
        if label == 1 or label == 2:  # F.C.L.r.asc.
            axisID = 40
            isLon = True
        elif label == 9 or label == 10:  # F.Cal.ant.-Sc.Cal.
            axisID = 281
            isLon = True
        elif label == 19 or label == 20:  # F.P.O.
            axisID = 297
            isLon = True
        elif label == 25 or label == 26:  # S.C.
            axisID = 360
            isLon = True
        elif label == 41 or label == 42:  # S.F.orbitaire.
            axisID = 61
            isLon = True
        elif label == 43 or label == 44:  # S.F.marginal.
            axisID = 61
            isLon = True
        elif label == 15 or label == 16:  # F.I.P.Po.C.inf.
            axisID = 339
            isLon = True
        elif label == 63 or label == 64:  # S.Po.C.sup.
            axisID = 339
            isLon = True
        elif label == 57 or label == 58:  # S.Pe.C.inf.
            axisID = 16
            isLon = True
        elif label == 59 or label == 60:  # S.Pe.C.median.
            axisID = 16
            isLon = True
        elif label == 61 or label == 62:  # S.Pe.C.sup.
            axisID = 16
            isLon = True
        #  latitude labels
        elif label == 5 or label == 6:  # F.C.M.ant.
            axisID = 55
            isLat = True
        elif label == 11 or label == 12:  # F.Coll. 55
            axisID = 56
            isLat = True
        elif label == 31 or label == 32:  # S.Call.
            axisID = 1
            isLat = True
        elif label == 33 or label == 34:  # S.F.inf.
            axisID = 106
            isLat = True
        elif label == 39 or label == 40:  # S.F.inter.
            axisID = 92
            isLat = True
        elif label == 45 or label == 46:  # S.F.sup.
            axisID = 81
            isLat = True
        elif label == 53 or label == 54:  # S.O.T.lat.post.
            axisID = 81
            isLat = True
        elif label == 55 or label == 56:  # S.Olf.
            axisID = 81
            isLat = True
        elif label == 65 or label == 66:  # S.T.i.ant.
            axisID = 92
            isLat = True
        elif label == 67 or label == 68:  # S.T.i.post.
            axisID = 92
            isLat = True
        elif label == 71 or label == 72:  # S.T.s.
            axisID = 106
            isLat = True
        elif label == 73 or label == 74:  # S.T.s.ter.asc.ant.
            axisID = 106
            isLat = True
        elif label == 75 or label == 76:  # S.T.s.ter.asc.post.
            axisID = 92
            isLat = True
        else:
            print 'no axis defined for label ', label
        return(axisID, isLon, isLat)
    
####################################################################
#
# defines the association between each sulcus and its corresponding axis in the model
#
####################################################################
    def sulcus2Axis(self, sulc_name_in):
        isLat = False
        isLon = False
        axisID = None
        if sulc_name_in.find('left')>0:
            sulc_name = sulc_name_in[:-5]
        elif sulc_name_in.find('right')>0:
            sulc_name = sulc_name_in[:-6]
        else:
            sulc_name = sulc_name_in
        #  longitude labels
        if sulc_name == 'F.C.L.r.asc.':
            axisID = 40
            isLon = True
        elif sulc_name == 'F.Cal.ant.-Sc.Cal.':
            axisID = 281
            isLon = True
        elif sulc_name == 'F.P.O.':
            axisID = 297
            isLon = True
        elif sulc_name == 'S.C.':
            axisID = 360
            isLon = True
        elif sulc_name == 'S.F.orbitaire.':
            axisID = 61
            isLon = True
        elif sulc_name == 'S.F.marginal.':
            axisID = 61
            isLon = True
        elif sulc_name == 'F.I.P.Po.C.inf.':
            axisID = 339
            isLon = True
        elif sulc_name == 'S.Po.C.sup.':
            axisID = 339
            isLon = True
        elif sulc_name == 'S.Pe.C.inf.':
            axisID = 16
            isLon = True
        elif sulc_name == 'S.Pe.C.median.':
            axisID = 16
            isLon = True
        elif sulc_name == 'S.Pe.C.sup.':
            axisID = 16
            isLon = True
        #  latitude labels
        elif sulc_name == 'F.C.M.ant.':
            axisID = 55
            isLat = True
        elif sulc_name == 'F.Coll.':
            axisID = 56  #55
            isLat = True
        elif sulc_name == 'S.Call.':
            axisID = 1
            isLat = True
        elif sulc_name == 'S.F.inf.':
            axisID = 106
            isLat = True
        elif sulc_name == 'S.F.inter.':
            axisID = 92
            isLat = True
        elif sulc_name == 'S.F.sup.':
            axisID = 81
            isLat = True
        elif sulc_name == 'S.O.T.lat.post.':
            axisID = 81
            isLat = True
        elif sulc_name == 'S.Olf.':
            axisID = 81
            isLat = True
        elif sulc_name == 'S.T.i.ant.':
            axisID = 92
            isLat = True
        elif sulc_name == 'S.T.i.post.':
            axisID = 92
            isLat = True
        elif sulc_name == 'S.T.s.':
            axisID = 106
            isLat = True
        elif sulc_name == 'S.T.s.ter.asc.ant.':
            axisID = 106
            isLat = True
        elif sulc_name == 'S.T.s.ter.asc.post.':
            axisID = 92
            isLat = True
        else:
            print 'no axis defined for sulcus ', sulc_name
        return(axisID, isLon, isLat)

####################################################################
#
# set the coordinates of each axis from 'sulci' using 'method'
#
####################################################################
    def setAxisCoord(self, sulci=None, method=None):
        if method is None:  # setting axis coord as the weighted barycenter of corresponding sulci
            print 'barycenter'
            # longitudes
            self.longitudeAxisCoord = self.longitudeAxis[:]
            arrayLongitudeCstrAxis = np.array(sulci.longitudeCstrAxis)
            for ax_ind in range(len(self.longitudeAxis)):
                inds = np.where(arrayLongitudeCstrAxis == self.longitudeAxis[ax_ind])[0]
                if len(inds) == 0:
                    self.longitudeAxisCoord[ax_ind] = None
                else:
                    ax_verts = []
                    ax_barys = []
                    ax_weights = []
                    ax_verts_weights = []
                    for r_sc_ind in range(len(inds)):
                        ax_verts.append(sulci.sulcalLines[sulci.longitudeCstrIndex[inds[r_sc_ind]]].vertices)
                        ax_barys.append(sulci.sulcalLines[sulci.longitudeCstrIndex[inds[r_sc_ind]]].barycenter[0])
                        ax_weights.append(sulci.sulcalLines[sulci.longitudeCstrIndex[inds[r_sc_ind]]].weight)
                        ax_verts_weights.append(sulci.sulcalLines[sulci.longitudeCstrIndex[inds[r_sc_ind]]].weight * np.ones(sulci.sulcalLines[sulci.longitudeCstrIndex[inds[r_sc_ind]]].nbVertices))
                    self.longitudeAxisCoord[ax_ind] = (np.sum(np.array(ax_weights) * np.array(ax_barys))) / np.sum(ax_weights)
            # laitudes
            self.latitudeAxisCoord = self.latitudeAxis[:]
            arraylatitudeCstrAxis = np.array(sulci.latitudeCstrAxis)
            for ax_ind in range(len(self.latitudeAxis)):
                inds = np.where(arraylatitudeCstrAxis == self.latitudeAxis[ax_ind])[0]
                if len(inds) == 0:
                    self.latitudeAxisCoord[ax_ind] = None
                else:
                    ax_verts = []
                    ax_barys = []
                    ax_weights = []
                    ax_verts_weights = []
                    for r_sc_ind in range(len(inds)):
                        ax_verts.append(sulci.sulcalLines[sulci.latitudeCstrIndex[inds[r_sc_ind]]].vertices)
                        ax_barys.append(sulci.sulcalLines[sulci.latitudeCstrIndex[inds[r_sc_ind]]].barycenter[1])
                        ax_weights.append(sulci.sulcalLines[sulci.latitudeCstrIndex[inds[r_sc_ind]]].weight)
                        ax_verts_weights.append(sulci.sulcalLines[sulci.latitudeCstrIndex[inds[r_sc_ind]]].weight * np.ones(sulci.sulcalLines[sulci.latitudeCstrIndex[inds[r_sc_ind]]].nbVertices))
                    self.latitudeAxisCoord[ax_ind] = (np.sum(np.array(ax_weights) * np.array(ax_barys))) / np.sum(ax_weights)

        else:
            print 'method ',method,' is not implemented yet!'
            
# ####################################################################
# #
# # buildModel
# #
# ####################################################################
# def buildModel(list_mesh, list_texture_sulci):
# #    square_ratio = 4.5
#     length = 4.5
#     width = 1
#     neocortex_tex_value = 0
#     insula_tex_value = 180
#     cingular_tex_value = 1
#     SC_label = 25
# 
#     model = Model()
#     group_full_sulci = SulcalLinesSet()
#     nb_mesh = len(list_mesh)
#     for mesh_ind in range(nb_mesh):
#         mesh = list_mesh[mesh_ind]
#         texture_poles = list_texture_poles[mesh_ind]
#         texture_sulci = list_texture_sulci[mesh_ind]
#         
#         full_sulci = SulcalLinesSet()
#         full_sulci.extractFromTexture(texture_sulci, mesh)
# 
#         full_sulci.updateIndices(neocortex_indices)
#         vert = np.array(neoCortex_square.vertex())
#         full_sulci.updateVertices(vert)
# 
#         SC_ind = full_sulci.names.index(('S.C._'+side))   
#         SC_label = full_sulci.labels[SC_ind]
#         print 'SC_label: ', SC_label
#         full_sulci.sulcalLines[SC_ind].printArgs()
#         translation = -full_sulci.sulcalLines[SC_ind].barycenter[0]
#         vert = np.array(mesh.vertex())
#         vert[:, 0] = vert[:, 0] + translation # * np.ones(vert.shape[0])
#         full_sulci.updateVertices(vert)
#         full_sulci.sulcalLine2SulcalConstraint(model)
#         group_full_sulci.cat(full_sulci)
# 
#     model.setBoundary(vert[neoCortex_open_boundary[0][0], 0], vert[neoCortex_open_boundary[0][-1], 0], vert[neoCortex_open_boundary[2][0], 1], vert[neoCortex_open_boundary[0][0], 1])
#     model.setAxisCoord(group_full_sulci)
#     print 'model built from '+nb_mesh+' subjects'
#     model.printArgs()
# 
#     return model
