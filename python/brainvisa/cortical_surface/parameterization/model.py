'''
Created on Nov 16, 2012

@author: toz
'''
from __future__ import print_function

from __future__ import absolute_import
import numpy as np
from soma import aims
from six.moves import range
#from brainvisa.cortical_surface.parameterization import mapping as map
class Model(object):
    '''
    classdocs
    '''

    def __init__(self, modelVersion=None, left=None, right=None, top=None, bottom=None, longitudeAxisID=None, latitudeAxisID=None, longitudeAxisCoord=None, latitudeAxisCoord=None, longitudeAxisSulci=None, latitudeAxisSulci=None, insularPoleBoundaryCoord=None, cingularPoleBoundaryCoord=None):
        '''
        Constructor
        '''
        if modelVersion is None:
            self.modelVersion = 3
        else:
            self.modelVersion = modelVersion
        # AxisID
        if longitudeAxisID is None:
            self.longitudeAxisID = [16, 40, 61, 281, 297, 339, 360]
        else:
            self.longitudeAxisID = longitudeAxisID
        if latitudeAxisID is None:
            self.latitudeAxisID = [1, 55, 56, 81, 92, 106]
        else:
            self.latitudeAxisID = latitudeAxisID
        # AxisCoord
        if longitudeAxisCoord is None:
            self.longitudeAxisCoord = [None, None, None, None, None, None, None]
        else:
            self.longitudeAxisCoord = longitudeAxisCoord
        if latitudeAxisCoord is None:
            self.latitudeAxisCoord = [None, None, None, None, None, None]
        else:
            self.latitudeAxisCoord = latitudeAxisCoord
        # AxisSulci
        if longitudeAxisSulci is None:
            self.longitudeAxisSulci = {'F.C.L.r.asc.': (40,1.0), 'F.Cal.ant.-Sc.Cal.': (281,1.0), 'F.P.O.': (297,1.0), 'S.C.': (360,1.0), 'S.F.orbitaire.': (61,1.0), 'S.F.marginal.': (61,1.0), 'F.I.P.Po.C.inf.': (339,1.0), 'S.Po.C.sup.': (339,1.0), 'S.Pe.C.inf.': (16,1.0), 'S.Pe.C.median.': (16,1.0), 'S.Pe.C.sup.': (16,1.0)}
        else:
            self.longitudeAxisSulci = longitudeAxisSulci
        if latitudeAxisSulci is None:
            self.latitudeAxisSulci ={'F.C.M.ant.': (55,1.0) , 'F.Coll.': (56,1.0), 'S.Call.': (1,1.0), 'S.F.inf.': (106,1.0), 'S.F.inter.': (92,1.0), 'S.F.sup.': (81,1.0), 'S.O.T.lat.post.': (81,1.0), 'S.Olf.': (81,1.0), 'S.T.i.ant.': (92,1.0), 'S.T.i.post.': (92,1.0), 'S.T.s.': (106,1.0), 'S.T.s.ter.asc.ant.': (106,1.0), 'S.T.s.ter.asc.post.': (92,1.0)}
        else:
            self.latitudeAxisSulci = latitudeAxisSulci
        # boundary coordinates
        if left == None:
            self.left = 0.0
        else:
            self.left = left
        if right == None:
            self.right = 450
        else:
            self.right = right
        if top == None:
            self.top = 100
        else:
            self.top = top
        if bottom == None:
            self.bottom = 0.0
        else:
            self.bottom = bottom
        if insularPoleBoundaryCoord == None:
            self.insularPoleBoundaryCoord = 30.0
        else:
            self.insularPoleBoundaryCoord = insularPoleBoundaryCoord    
        if cingularPoleBoundaryCoord == None:
            self.cingularPoleBoundaryCoord = 30.0
        else:
            self.cingularPoleBoundaryCoord = cingularPoleBoundaryCoord    

            

    def printArgs(self):
        txt = 'modelVersion ' + str(self.modelVersion) + '\n'
#         for col in data[index,:]:
#             txt = txt+' '+str(col)
        txt = txt + 'left ' + str(self.left) + '\n'
        txt = txt + 'right ' + str(self.right) + '\n'
        txt = txt + 'top ' + str(self.top) + '\n'
        txt = txt + 'bottom ' + str(self.bottom) + '\n'
        txt = txt + 'insularPoleBoundaryCoord ' + str(self.insularPoleBoundaryCoord) + '\n'
        txt = txt + 'cingularPoleBoundaryCoord ' + str(self.cingularPoleBoundaryCoord) + '\n'
        txt_tmp = ','.join(str(i) for i in self.longitudeAxisID)
        txt = txt + 'longitudeAxisID '+ txt_tmp + '\n'
        txt_tmp = ','.join(str(i) for i in self.longitudeAxisCoord)
        txt = txt + 'longitudeAxisCoord ' + txt_tmp + '\n'
        if self.modelVersion == '0.1':
            pass
        else:
            txt_tmp = ','.join(i+':('+str(self.longitudeAxisSulci[i][0])+';'+str(self.longitudeAxisSulci[i][1])+')' for i in self.longitudeAxisSulci)
            txt = txt + 'longitudeAxisSulci ' + txt_tmp + '\n'
        txt_tmp = ','.join(str(i) for i in self.latitudeAxisID)
        txt = txt + 'latitudeAxisID ' + txt_tmp + '\n'
        txt_tmp = ','.join(str(i) for i in self.latitudeAxisCoord)
        txt = txt + 'latitudeAxisCoord ' + txt_tmp + '\n'
        if self.modelVersion == '0.1':
            pass
        else:
            txt_tmp = ','.join(i+':('+str(self.latitudeAxisSulci[i][0])+';'+str(self.latitudeAxisSulci[i][1])+')' for i in self.latitudeAxisSulci)
            txt = txt + 'latitudeAxisSulci ' + txt_tmp + '\n'

        # txt = 'modelVersion ' + str(self.modelVersion) + '\n'
        # txt = txt + 'left ' + str(self.left) + '\n'
        # txt = txt + 'right ' + str(self.right) + '\n'
        # txt = txt + 'top ' + str(self.top) + '\n'
        # txt = txt + 'bottom ' + str(self.bottom) + '\n'
        # txt = txt + 'insularPoleBoundaryCoord ' + str(self.insularPoleBoundaryCoord) + '\n'
        # txt = txt + 'cingularPoleBoundaryCoord ' + str(self.cingularPoleBoundaryCoord) + '\n'
        # txt = txt + 'longitudeAxisID '+ str(self.longitudeAxisID) + '\n'
        # txt = txt + 'longitudeAxisCoord ' + str(self.longitudeAxisCoord) + '\n'
        # txt = txt + 'longitudeAxisSulci ' + str(self.longitudeAxisSulci) + '\n'
        # txt = txt + 'latitudeAxisID ' + str(self.latitudeAxisID) + '\n'
        # txt = txt + 'latitudeAxisCoord ' + str(self.latitudeAxisCoord) + '\n'
        # txt = txt + 'latitudeAxisSulci ' + str(self.latitudeAxisSulci) + '\n'

        return txt

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
#         txt = 'modelVersion ' + str(self.modelVersion) + '\n'
# #         for col in data[index,:]:
# #             txt = txt+' '+str(col)
#         txt = txt + 'left ' + str(self.left) + '\n'
#         txt = txt + 'right ' + str(self.right) + '\n'
#         txt = txt + 'top ' + str(self.top) + '\n'
#         txt = txt + 'bottom ' + str(self.bottom) + '\n'
#         txt = txt + 'insularPoleBoundaryCoord ' + str(self.insularPoleBoundaryCoord) + '\n'
#         txt = txt + 'cingularPoleBoundaryCoord ' + str(self.cingularPoleBoundaryCoord) + '\n'
#         txt_tmp = ','.join(str(i) for i in self.longitudeAxisID)
#         txt = txt + 'longitudeAxisID '+ txt_tmp + '\n'
#         txt_tmp = ','.join(str(i) for i in self.longitudeAxisCoord)
#         txt = txt + 'longitudeAxisCoord ' + txt_tmp + '\n'
#         txt_tmp = ','.join(str(i) for i in self.latitudeAxisID)
#         txt = txt + 'latitudeAxisID ' + txt_tmp + '\n'
#         txt_tmp = ','.join(str(i) for i in self.latitudeAxisCoord)
#         txt = txt + 'latitudeAxisCoord ' + txt_tmp + '\n'
#        f.write(txt)
        f.write(self.printArgs())
        f.close()

####################################################################
#
# read model from file and return a dictionnary
#
####################################################################
    def read(self, input_file):
        try: 
            txt_list = []
            with open(input_file,'r') as inf:
                for line in inf:
                    txt_list.append(line.split())
            txt_dict = dict((key, value) for (key, value) in txt_list)

            txt_dict['right'] = float(txt_dict['right'])
            txt_dict['left'] = float(txt_dict['left'])
            txt_dict['top'] = float(txt_dict['top'])
            txt_dict['bottom'] = float(txt_dict['bottom'])
            txt_dict['insularPoleBoundaryCoord'] = float(txt_dict['insularPoleBoundaryCoord'])
            txt_dict['cingularPoleBoundaryCoord'] = float(txt_dict['cingularPoleBoundaryCoord'])

            data_txt = txt_dict['longitudeAxisID']
            data_num = []
            for x in data_txt.split(','):
                try:
                    data_num.append(int(x))
                except:
                    data_num.append(None)
            txt_dict['longitudeAxisID'] = data_num

            data_txt = txt_dict['latitudeAxisID']
            data_num = []
            for x in data_txt.split(','):
                try:
                    data_num.append(int(x))
                except:
                    data_num.append(None)
            txt_dict['latitudeAxisID'] = data_num

            data_txt = txt_dict['longitudeAxisCoord']
            data_num = []
            for x in data_txt.split(','):
                try:
                    data_num.append(float(x))
                except:
                    data_num.append(None)
            txt_dict['longitudeAxisCoord'] = data_num

            data_txt = txt_dict['latitudeAxisCoord']
            data_num = []
            for x in data_txt.split(','):
                try:
                    data_num.append(float(x))
                except:
                    data_num.append(None)
            txt_dict['latitudeAxisCoord'] = data_num
            if txt_dict['modelVersion'] == '0.1':
                output_model = Model(txt_dict['modelVersion'], txt_dict['left'], txt_dict['right'], txt_dict['top'], txt_dict['bottom'], txt_dict['longitudeAxisID'], txt_dict['latitudeAxisID'], txt_dict['longitudeAxisCoord'], txt_dict['latitudeAxisCoord'], None, None, txt_dict['insularPoleBoundaryCoord'], txt_dict['cingularPoleBoundaryCoord'])
            elif txt_dict['modelVersion'] == '2':
                print('read model version 2')
                data_txt = txt_dict['longitudeAxisSulci']
                data_num = []
                for x in data_txt.split(','):
                    data_num.append(x.split(':'))
                longitudeAxisSulci = dict((key, int(value)) for (key, value) in data_num)
                data_txt = txt_dict['latitudeAxisSulci']
                data_num = []
                for x in data_txt.split(','):
                    data_num.append(x.split(':'))
                latitudeAxisSulci = dict((key, int(value)) for (key, value) in data_num)
                output_model = Model(txt_dict['modelVersion'], txt_dict['left'], txt_dict['right'], txt_dict['top'], txt_dict['bottom'], txt_dict['longitudeAxisID'], txt_dict['latitudeAxisID'], txt_dict['longitudeAxisCoord'], txt_dict['latitudeAxisCoord'], longitudeAxisSulci, latitudeAxisSulci, txt_dict['insularPoleBoundaryCoord'], txt_dict['cingularPoleBoundaryCoord'])
            elif txt_dict['modelVersion'] == '3':
                print('read model version 3')
                data_txt = txt_dict['longitudeAxisSulci']
                data_num = []
                for x in data_txt.split(','):
                    y = x.split(':')
                    z = y[1].split(';')
                    data_num.append([y[0],z[0][1:],z[1][:-1]])
                longitudeAxisSulci = dict((key, (int(value),float(weight))) for (key, value, weight) in data_num)
                data_txt = txt_dict['latitudeAxisSulci']
                data_num = []
                for x in data_txt.split(','):
                    y = x.split(':')
                    z = y[1].split(';')
                    data_num.append([y[0],z[0][1:],z[1][:-1]])
                latitudeAxisSulci = dict((key, (int(value),float(weight))) for (key, value, weight) in data_num)
                output_model = Model(txt_dict['modelVersion'], txt_dict['left'], txt_dict['right'], txt_dict['top'], txt_dict['bottom'], txt_dict['longitudeAxisID'], txt_dict['latitudeAxisID'], txt_dict['longitudeAxisCoord'], txt_dict['latitudeAxisCoord'], longitudeAxisSulci, latitudeAxisSulci, txt_dict['insularPoleBoundaryCoord'], txt_dict['cingularPoleBoundaryCoord'])
            else:
                raise Exception('cannot read the model file :: wrong model version')
        except:
            raise Exception('cannot read the model file')

        return output_model

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
                plt.plot([ax, ax], [self.top, self.bottom], 'k')
        for ax in self.latitudeAxisCoord:
            if ax is not None:
                plt.plot([self.left, self.right], [ax, ax], 'k')
        plt.plot([self.left, self.right], [self.top, self.top], 'k', linewidth=3.0)
        plt.plot([self.right, self.right], [self.top, self.bottom], 'k', linewidth=3.0)
        plt.plot([self.left, self.right], [self.bottom, self.bottom], 'k', linewidth=3.0)
        plt.plot([self.left, self.left], [self.top, self.bottom], 'k', linewidth=3.0)

##########################################################################
# replaced by sulcus2Axis, should not be used anymore
##########################################################################
    # def label2Axis(self, label):
    #     isLat = False
    #     isLon = False
    #     axisID = None
    #     #  longitude labels
    #     if label == 1 or label == 2:  # F.C.L.r.asc.
    #         axisID = 40
    #         isLon = True
    #     elif label == 9 or label == 10:  # F.Cal.ant.-Sc.Cal.
    #         axisID = 281
    #         isLon = True
    #     elif label == 19 or label == 20:  # F.P.O.
    #         axisID = 297
    #         isLon = True
    #     elif label == 25 or label == 26:  # S.C.
    #         axisID = 360
    #         isLon = True
    #     elif label == 41 or label == 42:  # S.F.orbitaire.
    #         axisID = 61
    #         isLon = True
    #     elif label == 43 or label == 44:  # S.F.marginal.
    #         axisID = 61
    #         isLon = True
    #     elif label == 15 or label == 16:  # F.I.P.Po.C.inf.
    #         axisID = 339
    #         isLon = True
    #     elif label == 63 or label == 64:  # S.Po.C.sup.
    #         axisID = 339
    #         isLon = True
    #     elif label == 57 or label == 58:  # S.Pe.C.inf.
    #         axisID = 16
    #         isLon = True
    #     elif label == 59 or label == 60:  # S.Pe.C.median.
    #         axisID = 16
    #         isLon = True
    #     elif label == 61 or label == 62:  # S.Pe.C.sup.
    #         axisID = 16
    #         isLon = True
    #     #  latitude labels
    #     elif label == 5 or label == 6:  # F.C.M.ant.
    #         axisID = 55
    #         isLat = True
    #     elif label == 11 or label == 12:  # F.Coll. 55
    #         axisID = 56
    #         isLat = True
    #     elif label == 31 or label == 32:  # S.Call.
    #         axisID = 1
    #         isLat = True
    #     elif label == 33 or label == 34:  # S.F.inf.
    #         axisID = 106
    #         isLat = True
    #     elif label == 39 or label == 40:  # S.F.inter.
    #         axisID = 92
    #         isLat = True
    #     elif label == 45 or label == 46:  # S.F.sup.
    #         axisID = 81
    #         isLat = True
    #     elif label == 53 or label == 54:  # S.O.T.lat.post.
    #         axisID = 81
    #         isLat = True
    #     elif label == 55 or label == 56:  # S.Olf.
    #         axisID = 81
    #         isLat = True
    #     elif label == 65 or label == 66:  # S.T.i.ant.
    #         axisID = 92
    #         isLat = True
    #     elif label == 67 or label == 68:  # S.T.i.post.
    #         axisID = 92
    #         isLat = True
    #     elif label == 71 or label == 72:  # S.T.s.
    #         axisID = 106
    #         isLat = True
    #     elif label == 73 or label == 74:  # S.T.s.ter.asc.ant.
    #         axisID = 106
    #         isLat = True
    #     elif label == 75 or label == 76:  # S.T.s.ter.asc.post.
    #         axisID = 92
    #         isLat = True
    #     else:
    #         print('no axis defined for label ', label)
    #     return(axisID, isLon, isLat)
    
####################################################################
#
# defines the association between each sulcus and its corresponding axis in the model
#
####################################################################
    def sulcus2Axis(self, sulc_name_in):
        isLat = False
        isLon = False
        axisID = None
        slWeight = 1.0
        if sulc_name_in.find('left')>0:
            sulc_name = sulc_name_in[:-5]
        elif sulc_name_in.find('right')>0:
            sulc_name = sulc_name_in[:-6]
        else:
            sulc_name = sulc_name_in

        if sulc_name in self.latitudeAxisSulci:
            isLat = True
            axisID = self.latitudeAxisSulci[sulc_name][0]
            slWeight = self.latitudeAxisSulci[sulc_name][1]
        elif sulc_name in self.longitudeAxisSulci:
            isLon = True
            axisID = self.longitudeAxisSulci[sulc_name][0]
            slWeight = self.longitudeAxisSulci[sulc_name][1]
        else:
            print('no axis defined for sulcus ', sulc_name)
        return(axisID, isLon, isLat, slWeight)

####################################################################
#
# set the coordinates of each axis from 'sulci' using 'method'
#
####################################################################
    def setAxisCoord(self, sulci=None, method=None):
        if method is None:  # setting axis coord as the weighted barycenter of corresponding sulci
            print('barycenter')
            # longitudes
            self.longitudeAxisCoord = self.longitudeAxisID[:]
            arrayLongitudeCstrAxis = np.array(sulci.longitudeCstrAxis)
            for ax_ind in range(len(self.longitudeAxisID)):
                inds = np.where(arrayLongitudeCstrAxis == self.longitudeAxisID[ax_ind])[0]
                if len(inds) == 0:
                    self.longitudeAxisCoord[ax_ind] = None
                else:
#                    ax_verts = []
                    ax_barys = []
                    ax_weights = []
#                    ax_verts_weights = []
                    for r_sc_ind in range(len(inds)):
#                        ax_verts.append(sulci.sulcalLines[sulci.longitudeCstrIndex[inds[r_sc_ind]]].vertices)
                        ax_barys.append(sulci.sulcalLines[sulci.longitudeCstrIndex[inds[r_sc_ind]]].barycenter[0])
                        ax_weights.append(sulci.sulcalLines[sulci.longitudeCstrIndex[inds[r_sc_ind]]].weight)
#                        ax_verts_weights.append(sulci.sulcalLines[sulci.longitudeCstrIndex[inds[r_sc_ind]]].weight * np.ones(sulci.sulcalLines[sulci.longitudeCstrIndex[inds[r_sc_ind]]].nbVertices))
                    self.longitudeAxisCoord[ax_ind] = (np.sum(np.array(ax_weights) * np.array(ax_barys))) / np.sum(ax_weights)
            # laitudes
            self.latitudeAxisCoord = self.latitudeAxisID[:]
            arraylatitudeCstrAxis = np.array(sulci.latitudeCstrAxis)
            for ax_ind in range(len(self.latitudeAxisID)):
                inds = np.where(arraylatitudeCstrAxis == self.latitudeAxisID[ax_ind])[0]
                if len(inds) == 0:
                    self.latitudeAxisCoord[ax_ind] = None
                else:
#                    ax_verts = []
                    ax_barys = []
                    ax_weights = []
#                    ax_verts_weights = []
                    for r_sc_ind in range(len(inds)):
#                        ax_verts.append(sulci.sulcalLines[sulci.latitudeCstrIndex[inds[r_sc_ind]]].vertices)
                        ax_barys.append(sulci.sulcalLines[sulci.latitudeCstrIndex[inds[r_sc_ind]]].barycenter[1])
                        ax_weights.append(sulci.sulcalLines[sulci.latitudeCstrIndex[inds[r_sc_ind]]].weight)
#                        ax_verts_weights.append(sulci.sulcalLines[sulci.latitudeCstrIndex[inds[r_sc_ind]]].weight * np.ones(sulci.sulcalLines[sulci.latitudeCstrIndex[inds[r_sc_ind]]].nbVertices))
                    self.latitudeAxisCoord[ax_ind] = (np.sum(np.array(ax_weights) * np.array(ax_barys))) / np.sum(ax_weights)

        else:
            print('method ',method,' is not implemented yet!')

####################################################################
#
# converts the coords of axes in the rectangle into latitude and longitude degree in [0-180] or [0-360]
#
####################################################################
    def axisCoordToDegree(self):
        from brainvisa.cortical_surface.parameterization.mapping import coordinatesFromRect

        precision_threshold = 0.00000001

        tmp_vert_lon = []
        for l in self.longitudeAxisCoord:
            if l is not None:
                tmp_vert_lon.append([l,0])
        tmp_vert_lat = []
        for l in self.latitudeAxisCoord:
            if l is not None:
                tmp_vert_lat.append([0,l])
        nb_lon = len(tmp_vert_lon)
        nb_lat = len(tmp_vert_lat)
        tmp_vert_boundaries = np.array([[self.left, self.top], [self.right, self.bottom]])
        
#        tmp_vert = np.concatenate(( np.concatenate((np.array(tmp_vert_lon), np.array(tmp_vert_lat)), 0), tmp_vert_boundaries), 0)
        tmp_vert = np.concatenate(( np.concatenate((np.array(tmp_vert_lat), np.array(tmp_vert_lon)), 0), tmp_vert_boundaries), 0)

        degree_lon, degree_lat = coordinatesFromRect(tmp_vert, self.insularPoleBoundaryCoord, self.cingularPoleBoundaryCoord)
        latitude_axis_coords = degree_lat[0:nb_lat]
        longitude_axis_coords = degree_lon[nb_lat:nb_lat+nb_lon+1]
        longitude_axis_coords[np.absolute(longitude_axis_coords)<precision_threshold] = 0
        longitude_axis_coords[np.absolute(longitude_axis_coords-360)<precision_threshold] = 0
#        longitude_axis_coords = degree_lon[0:nb_lon]
#        latitude_axis_coords = degree_lat[nb_lon:nb_lon+nb_lat]
        return (longitude_axis_coords, latitude_axis_coords)
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
#         print('SC_label: ', SC_label)
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
#     print('model built from '+nb_mesh+' subjects')
#     model.printArgs()
# 
#     return model
