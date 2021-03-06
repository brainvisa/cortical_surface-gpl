
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
from brainvisa.configuration import neuroConfig

name = 'Create Cortical Constraints Texture'

userLevel = 0

signature = Signature(
    'Side', Choice("Both","Left","Right"),
    'Lgraph', ReadDiskItem( 'Cortical folds graph', 'Graph',requiredAttributes={ 'side': 'left' } ),
    'Rgraph', ReadDiskItem( 'Cortical folds graph', 'Graph',requiredAttributes={ 'side': 'right' } ),
    'sulcus_identification',Choice('name','label'),
    'mri_corrected', ReadDiskItem( 'T1 MRI Bias Corrected', 'Aims readable volume formats' ),
    'left_white_mesh',ReadDiskItem( 'Left Hemisphere White Mesh' , 'aims Mesh Formats'),
    'right_white_mesh',ReadDiskItem( 'Right Hemisphere White Mesh' , 'aims Mesh Formats'),
    'left_white_sulci_mer',WriteDiskItem( 'Left hemisphere longitude constraints texture' ,'Texture',requiredAttributes={ 'side': 'left' } ),
    'left_white_sulci_par',WriteDiskItem( 'Left hemisphere latitude constraints texture' ,'Texture',requiredAttributes={ 'side': 'left' } ),
    'right_white_sulci_mer',WriteDiskItem( 'Right hemisphere longitude constraints texture' ,'Texture',requiredAttributes={ 'side': 'right' } ),
    'right_white_sulci_par',WriteDiskItem( 'Right hemisphere latitude constraints texture' ,'Texture',requiredAttributes={ 'side': 'right' } ),
    'translation',ReadDiskItem('Label Translation','Label Translation' ),
    'left_sulci_label_to_sulci_name',WriteDiskItem( 'Sulci To White Texture Translation', 'Text File',requiredAttributes={ 'side': 'left' }),
    'right_sulci_label_to_sulci_name',WriteDiskItem( 'Sulci To White Texture Translation', 'Text File',requiredAttributes={ 'side': 'right' }),
    'coord',Choice("Both","Longitude","Latitude"),
    'Metric',Choice('Euclidean','Internode')
)

def initialization( self ):
    self.linkParameters( 'left_white_mesh', 'Lgraph' )
    self.linkParameters( 'Rgraph', 'Lgraph' )
    self.linkParameters( 'right_white_mesh', 'Rgraph' )
    self.linkParameters( 'mri_corrected', 'left_white_mesh' )
    self.linkParameters( 'left_white_sulci_mer', 'left_white_mesh' )
    self.linkParameters( 'left_white_sulci_par', 'left_white_mesh' )
    self.linkParameters( 'right_white_sulci_mer', 'right_white_mesh' )
    self.linkParameters( 'right_white_sulci_par', 'right_white_mesh' )
    self.linkParameters( 'left_sulci_label_to_sulci_name', 'Lgraph' )
    self.linkParameters( 'right_sulci_label_to_sulci_name', 'Rgraph' )
    self.setOptional('Rgraph', 'Lgraph', 'left_white_mesh', 'right_white_mesh', 'left_white_sulci_mer','left_white_sulci_par', 'right_white_sulci_mer', 'right_white_sulci_par'  )
    self.sulcus_identification = 'label'
    self.translation = neuroConfig.getSharePath() + '/shfj-' + neuroConfig.shortVersion + '/nomenclature/translation/surfaceReferential.trl'
     
def execution( self, context ):
    Affine_estimation_coef = 0.9
    if self.Side in ('Left','Both'):
        context.write('Left hemisphere')
        if self.coord in ('Longitude','Both'):
            context.write('processing longitude constraints...')
            self.gyri_model = neuroConfig.getSharePath() + '/shfj-' + neuroConfig.shortVersion + '/models/3.0/gyrus/surfaceRef_mer.gyr'
            context.system('siMeshSulciProjection', '-i', self.left_white_mesh.fullPath() ,  '-g' , self.Lgraph.fullPath(), '-l' , self.translation.fullPath() , '-m', self.gyri_model , '-s' , self.sulcus_identification, '-v' , self.mri_corrected.fullPath(), '-o' ,self.left_white_sulci_mer.fullPath() ,'-V' , 1,'-M',2,'-n',5,'-a', Affine_estimation_coef, '-e', 20,'-t', self.left_sulci_label_to_sulci_name.fullPath() ,'-p',1 )
        if self.coord in ('Latitude','Both'):
            context.write('processing latitude constraints...')
            self.gyri_model = neuroConfig.getSharePath() + '/shfj-' + neuroConfig.shortVersion + '/models/3.0/gyrus/surfaceRef_par.gyr'
            context.system('siMeshSulciProjection', '-i', self.left_white_mesh.fullPath() ,  '-g' , self.Lgraph.fullPath(), '-l' , self.translation.fullPath() , '-m', self.gyri_model , '-s' , self.sulcus_identification, '-v' , self.mri_corrected.fullPath(), '-o' ,self.left_white_sulci_par.fullPath() ,'-V' , 1,'-M',2,'-n',5,'-a', Affine_estimation_coef, '-e', 20,'-t', self.left_sulci_label_to_sulci_name.fullPath() ,'-p',1 )
        context.write('Done')

    if self.Side in ('Right','Both'):
        context.write('Right hemisphere')
        if self.coord in ('Longitude','Both'):
            context.write('processing longitude constraints...')
            self.gyri_model = neuroConfig.getSharePath() + '/shfj-' + neuroConfig.shortVersion + '/models/3.0/gyrus/surfaceRef_mer.gyr'
            context.system('siMeshSulciProjection', '-i', self.right_white_mesh.fullPath() ,  '-g' , self.Rgraph.fullPath(), '-l' , self.translation.fullPath() , '-m', self.gyri_model , '-s' , self.sulcus_identification, '-v' , self.mri_corrected.fullPath(), '-o' ,self.right_white_sulci_mer.fullPath() ,'-V', 1,'-M',2,'-n',5,'-a', Affine_estimation_coef, '-e', 20,'-t',self.right_sulci_label_to_sulci_name.fullPath() ,'-p',1 )
        if self.coord in ('Latitude','Both'):
            context.write('processing latitude constraints...')
            self.gyri_model = neuroConfig.getSharePath() + '/shfj-' + neuroConfig.shortVersion + '/models/3.0/gyrus/surfaceRef_par.gyr'
            context.system('siMeshSulciProjection', '-i', self.right_white_mesh.fullPath() ,  '-g' , self.Rgraph.fullPath(), '-l' , self.translation.fullPath() , '-m', self.gyri_model , '-s' , self.sulcus_identification, '-v' , self.mri_corrected.fullPath(), '-o' ,self.right_white_sulci_par.fullPath() ,'-V', 1,'-M',2,'-n',5,'-a', Affine_estimation_coef, '-e', 20,'-t',self.right_sulci_label_to_sulci_name.fullPath() ,'-p',1 )
        context.write('Done')
