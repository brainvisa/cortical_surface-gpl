
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
from neuroProcesses import *
import shfjGlobals     

name = 'Create Right Hemisphere Cortical Constraints Texture'

userLevel = 2

signature = Signature(
    'Side', Choice("Right"),
    'Rgraph', ReadDiskItem( 'Cortical folds graph', 'Graph',requiredAttributes={ 'side': 'right' } ),
    'sulcus_identification',Choice('name','label'),
    'mri_corrected', ReadDiskItem( 'T1 MRI Bias Corrected', 'Aims readable volume formats' ),
    'right_white_mesh',ReadDiskItem( 'Right Hemisphere White Mesh' , shfjGlobals.aimsMeshFormats),
    'right_white_sulci_mer',WriteDiskItem( 'Right hemisphere longitude constraints texture' ,'Texture',requiredAttributes={ 'side': 'right' } ),
    'right_white_sulci_par',WriteDiskItem( 'Right hemisphere latitude constraints texture' ,'Texture',requiredAttributes={ 'side': 'right' } ),
    'translation',ReadDiskItem('Surface Label Translation','Label Translation' ),
    'ParModel',ReadDiskItem('Latitude Constraint Gyri Model','Gyri Model'),
    'MerModel',ReadDiskItem('Longitude Constraint Gyri Model','Gyri Model'),
    'right_sulci_label_to_sulci_name',WriteDiskItem( 'Sulci To White Texture Translation', 'Text File',requiredAttributes={ 'side': 'right' }),
    'coord',Choice("Both","Longitude","Latitude"),
    'Metric',Choice('Euclidean','Internode')
)

def initialization( self ):
    self.linkParameters( 'right_white_mesh', 'Rgraph' )
    self.linkParameters( 'mri_corrected', 'right_white_mesh' )
    self.linkParameters( 'right_white_sulci_mer', 'right_white_mesh' )
    self.linkParameters( 'right_white_sulci_par', 'right_white_mesh' )
    self.linkParameters( 'right_sulci_label_to_sulci_name', 'Rgraph' )
    self.setOptional('Rgraph', 'right_white_mesh', 'right_white_sulci_mer', 'right_white_sulci_par'  )
    self.sulcus_identification = 'label'
    self.findValue( 'translation', {} )
    self.setOptional('translation')
    self.findValue( 'ParModel', {} )
    self.setOptional('ParModel')
    self.findValue( 'MerModel', {} )
    self.setOptional('MerModel')
     
def execution( self, context ):
    Affine_estimation_coef = 0.9
    context.write('Right hemisphere')
    if self.coord in ('Longitude','Both'):
        context.write('processing longitude constraints...')
        context.system('siMeshSulciProjection', '-i', self.right_white_mesh.fullPath() ,  '-g' , self.Rgraph.fullPath(), '-l' , self.translation.fullPath() , '-m', self.MerModel.fullPath() , '-s' , self.sulcus_identification, '-v' , self.mri_corrected.fullPath(), '-o' ,self.right_white_sulci_mer.fullPath() ,'-V', 1,'-M',2,'-n',5,'-a', Affine_estimation_coef, '-e', 20,'-t',self.right_sulci_label_to_sulci_name.fullPath() ,'-p',1 )
    if self.coord in ('Latitude','Both'):
        context.write('processing latitude constraints...')
        context.system('siMeshSulciProjection', '-i', self.right_white_mesh.fullPath() ,  '-g' , self.Rgraph.fullPath(), '-l' , self.translation.fullPath() , '-m', self.ParModel.fullPath() , '-s' , self.sulcus_identification, '-v' , self.mri_corrected.fullPath(), '-o' ,self.right_white_sulci_par.fullPath() ,'-V', 1,'-M',2,'-n',5,'-a', Affine_estimation_coef, '-e', 20,'-t',self.right_sulci_label_to_sulci_name.fullPath() ,'-p',1 )
    context.write('Done')
