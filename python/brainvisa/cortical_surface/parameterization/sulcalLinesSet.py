'''
Created on 2 august 2012

@author: toz
'''
from soma import aims
import numpy as np
from brainvisa.cortical_surface.parameterization import sulcalLine as sln


class SulcalLinesSet(object):
    '''
    classdocs
    '''

    def __init__(self, sls=[sln.SulcalLine()]):
        '''
        Constructor
        '''
        self.sulcalLines = []
        self.labels = []
        self.names = []
        self.longitudeCstrIndex = []
        self.latitudeCstrIndex = []
        self.longitudeCstrAxis = []
        self.latitudeCstrAxis = []

        self.nbSulci = 0
        for sl in sls:
            if isinstance(sl, sln.SulcalConstraint):
                if sl.axisID == []:
                    print 'should not happen!!'
                    sl.label2Axis()
                if sl.isLon:
                    self.longitudeCstrIndex.append(self.nbSulci)
                    self.longitudeCstrAxis.append(sl.axisID)
                elif sl.isLat:
                    self.latitudeCstrIndex.append(self.nbSulci)
                    self.latitudeCstrAxis.append(sl.axisID)
            self.sulcalLines.append(sl)
            self.labels.append(sl.label)
            self.names.append(sl.name)
            self.nbSulci += 1

    def cat(self, sl_set):
        for sl in sl_set.sulcalLines:
            place = self.labels.index(sl.label)
            if place is None:
                if isinstance(sl, sln.SulcalConstraint):
                    if sl.axisID == []:
                        print 'should not happen!!'
                        sl.label2Axis()
                    if sl.isLon:
                        self.longitudeCstrIndex.append(self.nbSulci)
                    elif sl.isLat:
                        self.latitudeCstrIndex.append(self.nbSulci)
                self.sulcalLines.append(sl)
                self.labels.append(sl.label)
                self.nbSulci += 1
            else:
                self.sulcalLines[place].cat(sl)

    def updateVertices(self, vertices=None):
        for sl_ind in range(self.nbSulci):
            self.sulcalLines[sl_ind].updateVertices(vertices[self.sulcalLines[sl_ind].indices, :])
            self.sulcalLines[sl_ind].computeAttributes()

#     def updateIndices(self, neocortex_indices=None):
#         for sl_ind in range(self.nbSulci):
#             new_indices = []
#             for i in self.sulcalLines[sl_ind].indices:
#                 try:
#                     curr_ind = neocortex_indices.index(i)
#                     new_indices.append(curr_ind)
#                 except:
#                     # when a new indice is not in sulcalLines[sl_ind].indices,
#                     # add a -1 to allow sulcalLine.updateIndices() to delete 
#                     # the corresponding segment in self.segm
#                     new_indices.append(-1)
#                     print 'some vertices of sulcus nb ', sl_ind, ' are not in the new indices'
# #            print new_indices
# #            print 'max = ',np.max(new_indices)
#             self.sulcalLines[sl_ind].updateIndices(np.array(new_indices, np.int32))

#     def label2Axis(self):
#         cstr_indices = [self.longitudeCstrIndex, self.latitudeCstrIndex]
#         for sl_ind in cstr_indices:
#             self.sulcalLines[sl_ind].label2Axis()

    def printArgs(self):
        print 'SulcalLinesSet ::'
        print '    labels = ', self.labels
        print '    names = ', self.names
        print '    longitudeCstrIndex = ', self.longitudeCstrIndex
        print '    longitudeCstrAxis = ', self.longitudeCstrAxis
        print '    latitudeCstrIndex = ', self.latitudeCstrIndex
        print '    latitudeCstrAxis = ', self.latitudeCstrAxis
        print '    nbSulci = ', self.nbSulci
        print '    sulcalLines = ', self.sulcalLines

    def extractFromTexture(self, tex, mesh, sulc_labels_dict=None, labels=None, neigh=None):
        if neigh is None:
            neigh = aims.SurfaceManip.surfaceNeighbours(mesh)
        if labels is None:
            if tex.dtype != 'int16':
                print 'warning :: sulci texture type is not int16'
                atex = np.around(tex)
            else:
                atex = tex
            labels = np.unique(atex)
            labels = labels.tolist()
            labels.remove(0)
        sls = [sln.SulcalLine() for l in range(len(labels))]
        for l in range(len(labels)):
            sls[l].extractFromTexture(labels[l], tex, mesh, sulc_labels_dict, neigh)
        return self.__init__(sls)

    def sulcalLine2SulcalConstraint(self, modele=None, names=None):
        if modele is None:
            print 'no modele given, nothing to do!'
            print 'try sulcalLine2SulcalConstraint(modele)'
        else:
            if names is None:
                names = self.names
            for name in names:
                place = self.names.index(name)
                if place is None:
                    print 'sulcus ' + name + ' not present in this Sulci object!!'
                else:
                    (axisID, isLon, isLat) = modele.sulcus2Axis(name)
                    if isLon:
                        self.sulcalLines[place] = sln.SulcalConstraint(self.sulcalLines[place], axisID, isLon, isLat)
                        self.longitudeCstrIndex.append(place)
                        self.longitudeCstrAxis.append(axisID)
                    elif isLat:
                        self.sulcalLines[place] = sln.SulcalConstraint(self.sulcalLines[place], axisID, isLon, isLat)
                        self.latitudeCstrIndex.append(place)
                        self.latitudeCstrAxis.append(axisID)
                    else:
                        print 'sulcalLine with name ', name, 'is not a constraint'

    def toMesh(self):

        outplst = []
        outv_tmp=[]
        #outn_tmp=[]
        n = 0
        for sl_ind in range(self.nbSulci):
            sl_mesh = self.sulcalLines[sl_ind].toMesh()
            outv_tmp.append(np.array(sl_mesh.vertex()))
            outplst.append(np.array(sl_mesh.polygon()) + n)
            n += len(sl_mesh.vertex())
            #outn_tmp.append(sl_mesh.normals())
#             out_mesh.vertex(sl_ind).assign(sl_mesh.vertex())
#             out_mesh.polygon(sl_ind).assign(sl_mesh.polygon())
        out_mesh = aims.AimsTimeSurface_2()
        vv = aims.vector_POINT3DF()
        vp = aims.vector_AimsVector_U32_2()
        for x in np.vstack(outv_tmp):
            vv.append(x)    
        for x in np.vstack(outplst):
            vp.append(x)
        out_mesh.vertex().assign(vv) 
        out_mesh.polygon().assign(vp)
        out_mesh.updateNormals()
        return out_mesh

    def toTex(self):
        out_tex = aims.TimeTexture_S16()
#        out_tex = []
        for sl_ind in range(self.nbSulci):
            sl_tex = self.sulcalLines[sl_ind].toTex()
            for i in sl_tex[0]:
                out_tex[0].append(i)
#            out_tex[sl_ind] = sl_tex[0]
#            out_tex.extend(sl_tex[0])         
        return out_tex

    def plot(self, plt, modele=None):
        colors = ['b', 'g', 'r', 'm', 'y', 'k']
        i_col = 0
        for sl in self.sulcalLines:
            sl.plot(plt, colors[i_col], modele)
            i_col += 1
            if i_col > len(colors) - 1:
                i_col = 0

    def save(self, fileNameMesh, fileNameTex):
        ws = aims.Writer()
        ws.write(self.toMesh(), fileNameMesh)
        ws.write(self.toTex(), fileNameTex)