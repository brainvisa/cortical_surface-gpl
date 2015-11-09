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

from brainvisa.processes import *
from soma import aims
import numpy as np
import os
try:
    from brainvisa.cortical_surface.surface_tools import PDE_tools as pdeTls
except:
    pass

name = 'Parcels Volume'
userLevel = 0

signature = Signature(
#    'side', Choice('left', 'right'),
    'volumes_parcels_left', ListOf(
        ReadDiskItem( 'Left Gyri Volume', 'Aims writable volume formats' ) ), #parcellation volume', 'Aims writable volume formats' ),
    'volumes_parcels_right', ListOf(
        ReadDiskItem( 'Right Gyri Volume', 'Aims writable volume formats' ) ), #parcellation volume', 'Aims writable volume formats' ),
#        ReadDiskItem('Parcellation volume',  'Aims writable volume formats'  ) ),
#    'normalization_by_total_GM_volume', Choice('no', 'yes'),
    'background_label', Integer(),
    'output_csv_file', WriteDiskItem( 'CSV file', 'CSV file' )
)

# def linkParcels(proc, dummy):
#     proc.signature['textures_parcels'] = ListOf(
#         ReadDiskItem('hemisphere marsAtlas parcellation texture',
#                      'aims Texture formats',
#                      requiredAttributes={'regularized': 'false',
#                                          'side': proc.side}))
#     proc.changeSignature(proc.signature)

def initialization( self ):
    self.linkParameters('volumes_parcels_right','volumes_parcels_left' )
    self.background_label = 0
#    self.linkParameters(None, 'side', linkParcels)
#    self.linkParameters('white_meshes','textures_parcels' )


def execution( self, context ):
    def compute_parcels_volume(avol, vox_volume):
        #context.write(avol.shape)
        #context.write(avol)
        bins=np.bincount(avol.ravel())
        parcels_labels = np.nonzero(bins)[0]
        counts = bins[parcels_labels]
        #context.write(parcels_labels)
        parcels_volume = counts*vox_volume
        return parcels_volume,  parcels_labels


    parcels_volume = list()
    parcels_label = list()
    subjects = list()
    for ind_vol, vol_file in enumerate(self.volumes_parcels_left):
        subj = vol_file.get('subject')
        subjects.append(subj)
        context.write('working on subject nb: '+str(ind_vol+1)+', '+subj)
        vol = aims.read(vol_file.fullPath())
        vox_volume = np.prod(vol.header()[ 'voxel_size' ])
        context.write('voxel volume = '+str(vox_volume))
        avol = np.array( vol, copy=False )
        vol_r = aims.read(self.volumes_parcels_right[ind_vol].fullPath())
        avol_r = np.array( vol_r, copy=False )
        avol = avol + avol_r
        (subj_parcels_volume, subj_parcels_labels) = compute_parcels_volume(avol, vox_volume)
        parcels_volume.append(subj_parcels_volume)
        parcels_label.append(subj_parcels_labels)
    union_labels = np.array([], dtype=np.int32)
    for x in parcels_label:
        union_labels = np.append(union_labels, x)
    common_labels = np.unique(union_labels)
    ind_bck_label = np.where(common_labels == self.background_label)[0]
    if ind_bck_label.size == 0:
        context.write('background_label not found in the volumes')
    else:
        common_labels = np.delete(common_labels, ind_bck_label)
    context.write('common_labels : ')
    context.write(common_labels)
    f = open( self.output_csv_file.fullPath(), 'w' )
    f.write( 'subject' )
    for lab in common_labels:
        f.write( ','+str(lab) )
    f.write('\n')
    for ind_s, subj in enumerate(subjects):
        f.write( subj )
        for lab in common_labels:
            ind_l = np.where(parcels_label[ind_s] == lab)[0]
            if ind_l.size > 0:
                f.write( ','+str(parcels_volume[ind_s][ind_l][0]) )
            else:
                f.write( ','+str(0) )
        f.write('\n')
    f.close()

