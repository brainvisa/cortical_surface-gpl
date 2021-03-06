#  This software and supporting documentation are distributed by
#      Institut Federatif de Recherche 49
#      CEA/NeuroSpin, Batiment 145,
#      91191 Gif-sur-Yvette cedex
#      France
#
# This software is governed by the CeCILL license version 2 under
# French law and abiding by the rules of distribution of free software.
# You can  use, modify and/or redistribute the software under the 
# terms of the CeCILL license version 2 as circulated by CEA, CNRS
# and INRIA at the following URL "http://www.cecill.info". 
#
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

from __future__ import absolute_import
from brainvisa.processes import *
from soma import aims
import numpy as np
import os
try:
    from brainvisa.cortical_surface.surface_tools import PDE_tools as pdeTls
except:
    pass

name = 'Parcels Area'
userLevel = 0

signature = Signature(
    'side', Choice('left', 'right'),
    'textures_parcels', ListOf(
        ReadDiskItem('hemisphere parcellation texture',
                     'aims Texture formats',
                     requiredAttributes={'regularized': 'false'}) ),
    'white_meshes',ListOf(ReadDiskItem( 'Hemisphere White Mesh',
                                       'aims mesh formats')),
    'parcels_to_ignore', ListOf( Integer() ),
    'normalization_by_total_area', Choice('no', 'yes'),
    'output_csv_file', WriteDiskItem( 'CSV file', 'CSV file' )
)

def linkParcels(proc, dummy):
    proc.signature['textures_parcels'] = ListOf(
        ReadDiskItem('hemisphere parcellation texture',
                     'aims Texture formats',
                     requiredAttributes={'regularized': 'false',
                                         'side': proc.side}))
    proc.changeSignature(proc.signature)

def initialization( self ):
    self.linkParameters(None, 'side', linkParcels)
    self.linkParameters('white_meshes','textures_parcels' )
    self.parcels_to_ignore = [0,255]

def execution( self, context ):
    def compute_parcels_area(mesh, atex_parcels, bad_labels):
        vert_voronoi = pdeTls.vertexVoronoi(mesh)
        parcels_labels = list(np.unique(atex_parcels))
        for l in bad_labels:
            parcels_labels.remove(l)
        parcels_area = np.zeros(len(parcels_labels))
        for ind_lab, lab in enumerate(parcels_labels):
            vert_inds = np.where(atex_parcels == lab)[0]
            parcels_area[ind_lab] = np.sum(vert_voronoi[vert_inds])
        total_area = np.sum(parcels_area)
        return parcels_area,  parcels_labels, total_area


    parcels_area = list()
    total_area = list()
    subjects = list()
    parcels_label = list()
    for ind_mesh, r_mesh in enumerate(self.white_meshes):
        context.write('working on mesh nb: ',ind_mesh+1)
        subjects.append(r_mesh.get('subject'))
        mesh = aims.read(r_mesh.fullPath())
        tex_parcels = aims.read(self.textures_parcels[ind_mesh].fullPath())
        atex_parcels = tex_parcels[0].arraydata()
        (subj_parcels_area, subj_parcels_labels, subj_total_area) = compute_parcels_area(mesh, atex_parcels, self.parcels_to_ignore)
        parcels_area.append(subj_parcels_area)
        total_area.append(subj_total_area)
        parcels_label.append(subj_parcels_labels)
    union_labels = np.array([], dtype=np.int32)
    for x in parcels_label:
        union_labels = np.append(union_labels, x)
    common_labels = np.unique(union_labels)
    f = open( self.output_csv_file.fullPath(), 'w' )
    f.write( 'subject,side,total_area' )
    for lab in common_labels:
        f.write( ','+str(lab) )
    f.write('\n')
    for ind_s, subj in enumerate(subjects):
        tot_a = total_area[ind_s]
        f.write( subj + ',' + self.side+ ',' + str(tot_a) )
        for lab in common_labels:
            ind_l = np.where(parcels_label[ind_s] == lab)[0]
            if ind_l.size > 0:
                if self.normalization_by_total_area == 'yes':
                    f.write( ','+str(parcels_area[ind_s][ind_l][0]/tot_a) )
                else:
                    f.write( ','+str(parcels_area[ind_s][ind_l][0]) )
            else:
                f.write( ','+str(0) )
        f.write('\n')
    f.close()

