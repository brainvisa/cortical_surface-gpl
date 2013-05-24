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
            self.nbSulci += 1

    def cat(self, sl):
        for lab in sl.labels:
            place = self.labels.index(lab)
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
                self.sulcalLines[place].computeAttributes()

    def updateVertices(self, vertices=None):
        for sl_ind in range(self.nbSulci):
            self.sulcalLines[sl_ind].updateVertices(vertices[self.sulcalLines[sl_ind].indices, :])
            self.sulcalLines[sl_ind].computeAttributes()

    def updateIndices(self, neocortex_indices=None):
        for sl_ind in range(self.nbSulci):
            new_indices = []
            for i in self.sulcalLines[sl_ind].indices:
                try:
                    curr_ind = neocortex_indices.index(i)
                    new_indices.append(curr_ind)
                except:
                    # when a new indice is not in sulcalLines[sl_ind].indices,
                    # add a -1 to allow sulcalLine.updateIndices() to delete 
                    # the corresponding segment in self.segm
                    new_indices.append(-1)
                    print 'some vertices of sulcus nb ', sl_ind, ' are not in the new indices'
#            print new_indices
#            print 'max = ',np.max(new_indices)
            self.sulcalLines[sl_ind].updateIndices(np.array(new_indices, np.int32))

    def label2Axis(self):
        cstr_indices = [self.longitudeCstrIndex, self.latitudeCstrIndex]
        for sl_ind in cstr_indices:
            self.sulcalLines[sl_ind].label2Axis()

    def printArgs(self):
        print 'SulcalLinesSet ::'
        print '    labels = ', self.labels
        print '    longitudeCstrIndex = ', self.longitudeCstrIndex
        print '    longitudeCstrAxis = ', self.longitudeCstrAxis
        print '    latitudeCstrIndex = ', self.latitudeCstrIndex
        print '    latitudeCstrAxis = ', self.latitudeCstrAxis
        print '    nbSulci = ', self.nbSulci
        print '    sulcalLines = ', self.sulcalLines

    def extractFromTexture(self, tex, mesh, labels=None, neigh=None):
        if neigh is None:
            neigh = aims.SurfaceManip.surfaceNeighbours(mesh)
        if labels is None:
            atex = np.around(tex)
            labels = np.unique(atex)
            labels = labels.tolist()
            labels.remove(0)
        sls = [sln.SulcalLine() for l in range(len(labels))]
        for l in range(len(labels)):
            sls[l].extractFromTexture(labels[l], tex, mesh, neigh)
        return self.__init__(sls)

    def sulcalLine2SulcalConstraint(self, modele=None, labels=None):
        if modele is None:
            print 'no modele given, nothing to do!'
            print 'try sulcalLine2SulcalConstraint(modele)'
        else:
            if labels is None:
                labels = self.labels
            for lab in labels:
                place = self.labels.index(lab)
                if place is None:
                    print 'label ' + lab + ' not present in this Sulci object!!'
                else:
                    (axisID, isLon, isLat) = modele.label2Axis(lab)
                    if isLon:
                        self.sulcalLines[place] = sln.SulcalConstraint(self.sulcalLines[place], axisID, isLon, isLat)
                        self.longitudeCstrIndex.append(place)
                        self.longitudeCstrAxis.append(axisID)
                    elif isLat:
                        self.sulcalLines[place] = sln.SulcalConstraint(self.sulcalLines[place], axisID, isLon, isLat)
                        self.latitudeCstrIndex.append(place)
                        self.latitudeCstrAxis.append(axisID)
                    else:
                        print 'sulcalLine with label ', lab, 'is not a constraint'

    def toMesh(self):
        out_mesh = aims.AimsTimeSurface_2()
        for sl_ind in range(self.nbSulci):
            sl_mesh = self.sulcalLines[sl_ind].toMesh()
            out_mesh.vertex(sl_ind).assign(sl_mesh.vertex())
            out_mesh.polygon(sl_ind).assign(sl_mesh.polygon())
        out_mesh.updateNormals()
        return out_mesh

    def toTex(self):
        out_tex = aims.TimeTexture_S16()
        for sl_ind in range(self.nbSulci):
            sl_tex = self.sulcalLines[sl_ind].toTex()
            out_tex[sl_ind] = sl_tex[0]
        return out_tex

    def plot(self, plt, modele=None):
        colors = ['b', 'g', 'r', 'm', 'y', 'k']
        i_col = 0
        for sl in self.sulcalLines:
            sl.plot(plt, colors[i_col], modele)
            i_col += 1
            if i_col > len(colors) - 1:
                i_col = 0

    def save(self, fileName):
        ws = aims.Writer()
        ws.write(self.toMesh(), fileName + '.mesh')
        ws.write(self.toTex(), fileName + '.tex')