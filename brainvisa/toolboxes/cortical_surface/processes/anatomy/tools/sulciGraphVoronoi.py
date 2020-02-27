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

from __future__ import print_function

from __future__ import absolute_import
from brainvisa.processes import *
from soma import aims
import numpy as np

name = 'Sulci Voronoi'

userLevel = 2


signature = Signature(
    'graph', ReadDiskItem( 'Labelled Cortical folds graph', 'Graph and data' ),
    'side', Choice('left', 'right'),
    'white_mesh',ReadDiskItem( 'Hemisphere White Mesh', 'aims mesh formats' ),
    'mri', ReadDiskItem( 'T1 MRI Bias Corrected', 'Aims readable volume formats' ),
    'sulcus_identification',Choice('name','label'),
    'labels_translation_map',ReadDiskItem( 'Label Translation' ,'Label Translation'),
    'bucket_label_type', Choice('All', 'aims_junction', 'aims_bottom', 'aims_ss', 'aims_other'),
    'sulci_voronoi',WriteDiskItem( 'sulci voronoi texture', 'aims Texture formats'),
    'input_int_to_label_translation', ReadDiskItem('log file', 'text file'),
    'int_to_label_translation', WriteDiskItem('log file', 'text file'),


)

def initialization( self ):
    def linkSide( proc, dummy ):
        if proc.graph is not None:
            side = proc.graph.get( 'side' )
            return side

    self.linkParameters( 'mri', 'white_mesh' )
    self.linkParameters( 'white_mesh', 'graph' )
    self.linkParameters( 'sulci_voronoi', 'white_mesh' )
    self.linkParameters( 'side', 'graph', linkSide )
    self.findValue( 'labels_translation_map', {'filename_variable' : 'sulci_model_2008'} )
    self.linkParameters('input_int_to_label_translation', 'graph')
    self.linkParameters('int_to_label_translation', 'graph')
    self.bucket_label_type = 'All'
    self.sulcus_identification = 'label'
    self.setOptional('input_int_to_label_translation')
    self.setOptional('int_to_label_translation')


def execution( self, context ):
    def nearest_neighbor(vert_template,vert_pits):
        vertex_number1=vert_template.shape[0]
        print('vert_template.shape', vert_template.shape[0])
        v_number=vert_pits.shape[0]
        print('vert_pits.shape', vert_pits.shape[0])
        nn=[]
        for v in vert_pits:
                nn_tmp = np.argmin(np.sum(np.square(np.tile(v,(vertex_number1,1))-vert_template),1))
                nn.append(nn_tmp)
        return nn

    volumeGraphLabelBasins=context.temporary('NIFTI-1 image')

    context.write('computing Graph Label Volume')

    if (self.bucket_label_type=='aims_junction'):
        graphBucketLabel = [ 'siGraph2Label','-g', self.graph.fullPath(),'-a', self.sulcus_identification,'-tv', self.mri.fullPath(),'-tr', self.labels_translation_map.fullPath(),'-o', volumeGraphLabelBasins.fullPath(),
      '-b', 'aims_junction']
    elif (self.bucket_label_type=='aims_bottom'):
        graphBucketLabel = [ 'siGraph2Label','-g', self.graph.fullPath(),'-a', self.sulcus_identification,'-tv', self.mri.fullPath(),'-tr', self.labels_translation_map.fullPath(),'-o', volumeGraphLabelBasins.fullPath(),
      '-b', 'aims_bottom']
    elif (self.bucket_label_type=='aims_ss'):
        graphBucketLabel = [ 'siGraph2Label','-g', self.graph.fullPath(),'-a', self.sulcus_identification,'-tv', self.mri.fullPath(),'-tr', self.labels_translation_map.fullPath(),'-o', volumeGraphLabelBasins.fullPath(),
      '-b', 'aims_ss']
    elif (self.bucket_label_type=='aims_other'):
        graphBucketLabel = [ 'siGraph2Label','-g', self.graph.fullPath(),'-a', self.sulcus_identification,'-tv', self.mri.fullPath(),'-tr', self.labels_translation_map.fullPath(),'-o', volumeGraphLabelBasins.fullPath(),
      '-b', 'aims_other']
    else :
        graphBucketLabel = [ 'siGraph2Label','-g', self.graph.fullPath(),'-a', self.sulcus_identification,'-tv', self.mri.fullPath(),'-tr', self.labels_translation_map.fullPath(),'-o', volumeGraphLabelBasins.fullPath(),
      '-b', 'aims_junction','-b', 'aims_bottom','-b', 'aims_ss', '-b', 'aims_other']

    if self.input_int_to_label_translation is not None:
        graphBucketLabel += ['-it', self.input_int_to_label_translation]
    if self.int_to_label_translation:
        graphBucketLabel += ['-ot', self.int_to_label_translation]

    context.system(*graphBucketLabel)

    context.write('nearest-neighbor interpolation between mesh and Graph Label Volume')
    vol_sulci = aims.read(volumeGraphLabelBasins.fullPath())
    avol_sulci = np.array( vol_sulci, copy=False )
    vs = vol_sulci.header()[ 'voxel_size' ]
    (x, y, z) = np.where( avol_sulci[:, :, :, 0] > 0 )

    labels = avol_sulci[x, y, z,0]

    nb_verts = len(labels)
    vert_vol_sulci = np.zeros((nb_verts,3))
    # get coords in millimeters
    vert_vol_sulci[:, 0] = x * vs[0]
    vert_vol_sulci[:, 1] = y * vs[1]
    vert_vol_sulci[:, 2] = z * vs[2 ]

    mesh = aims.read(self.white_mesh.fullPath())
    vert_mesh = np.array(mesh.vertex())
    nn = nearest_neighbor(vert_vol_sulci, vert_mesh)

    a_tex_out = labels[nn]
    tex_out = aims.TimeTexture_S16()
    tex_out[0].assign(a_tex_out)
    aims.write(tex_out, self.sulci_voronoi.fullPath())
    context.write('Done')










