# Copyright CEA and IFR 49 (2000-2005)
#
#  This software and supporting documentation were developed by
#      CEA/DSV/SHFJ and IFR 49
#      4 place du General Leclerc
#      91401 Orsay cedex
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

name = 'Create Label Texture'

userLevel = 1

signature = Signature(
    'Side', Choice("Both","Left","Right"),
    'Lgraph', ReadDiskItem( 'Cortical folds graph', 'Graph',requiredAttributes={ 'side': 'left' } ),
    'Rgraph', ReadDiskItem( 'Cortical folds graph', 'Graph',requiredAttributes={ 'side': 'right' } ),
    'sulcus_identification',Choice('name','label'),
    'mri_corrected', ReadDiskItem( 'T1 MRI Bias Corrected', 'Aims readable volume formats' ),
    'left_white_mesh',ReadDiskItem( 'Left Hemisphere White Mesh' , shfjGlobals.aimsMeshFormats),
    'right_white_mesh',ReadDiskItem( 'Right Hemisphere White Mesh' , shfjGlobals.aimsMeshFormats),
    'MNI', Choice("No","Yes"),
    'mni_mesh', ReadDiskItem('MNI Cortex Mesh', 'Aims mesh formats'),
    'left_white_sulci',WriteDiskItem( 'Sulci White Texture' ,'Texture',requiredAttributes={ 'side': 'left' } ),
    'right_white_sulci',WriteDiskItem( 'Sulci White Texture' ,'Texture',requiredAttributes={ 'side': 'right' } ),
    'translation',ReadDiskItem('Label Translation','Label Translation' ),
    'left_sulci_label_to_sulci_name',WriteDiskItem( 'Sulci To White Texture Translation', 'Text File',requiredAttributes={ 'side': 'left' }),
    'right_sulci_label_to_sulci_name',WriteDiskItem( 'Sulci To White Texture Translation', 'Text File',requiredAttributes={ 'side': 'right' }),
    'gyri_model',ReadDiskItem('Gyri Model','Gyri Model' ),
    'Metric',Choice('Euclidean','Internode')
)

def initialization( self ):
     self.linkParameters( 'left_white_mesh', 'Lgraph' )
     self.linkParameters( 'Rgraph', 'Lgraph' )
     self.linkParameters( 'mri_corrected', 'Lgraph' )
     self.linkParameters( 'right_white_mesh', 'Rgraph' )
     self.linkParameters( 'left_white_sulci', 'left_white_mesh' )
     self.linkParameters( 'right_white_sulci', 'right_white_mesh' )
     self.linkParameters( 'mni_mesh', 'mri_corrected' )
     self.linkParameters( 'left_sulci_label_to_sulci_name', 'mri_corrected' )
     self.linkParameters( 'right_sulci_label_to_sulci_name', 'mri_corrected' )
     self.setOptional('mni_mesh','Rgraph', 'Lgraph', 'left_white_mesh', 'right_white_mesh', 'left_white_sulci', 'right_white_sulci'  )
     self.sulcus_identification = 'label'
     self.translation = '/home/Panabase/data_icbm/morphometry_gyri/Traduction_Gyri_icbm_16-05-03.trl'
     self.gyri_model = '/home/Panabase/data_icbm/morphometry_gyri/Gyri_icbm_16-05-03.gyr'
     
def execution( self, context ):
     Min_points_in_connected_component = 5
     Affine_estimation_coef = 0.9
     Distance_Euclidienne = 20
     Distance_au_plan     = 1
     Volume_Closing = 1
     Mesh_Closing = 2
     if self.MNI == 'Yes':
          self.left_white_mesh = self.mni_mesh
          self.right_white_mesh = self.mni_mesh

     if self.Metric == 'Euclidean':
          if self.Side in ('Left','Both'):
               context.system('siMeshSulciProjection', '-i', self.left_white_mesh.fullPath() ,  '-g' , self.Lgraph.fullPath(), '-l' , self.translation.fullPath() , '-m', self.gyri_model.fullPath() , '-s' , self.sulcus_identification, '-v' , self.mri_corrected.fullPath(), '-o' ,self.left_white_sulci.fullPath() ,'-V' , Volume_Closing,'-M',Mesh_Closing,'-n',Min_points_in_connected_component,'-a', Affine_estimation_coef, '-e', Distance_Euclidienne,'-t', self.left_sulci_label_to_sulci_name.fullPath() ,'-p',Distance_au_plan )
          if self.Side in ('Right','Both'):
               context.system('siMeshSulciProjection', '-i', self.right_white_mesh.fullPath() ,  '-g' , self.Rgraph.fullPath(), '-l' , self.translation.fullPath() , '-m', self.gyri_model.fullPath() , '-s' , self.sulcus_identification, '-v' , self.mri_corrected.fullPath(), '-o' ,self.right_white_sulci.fullPath() ,'-V', Volume_Closing,'-M',Mesh_Closing,'-n',Min_points_in_connected_component,'-a', Affine_estimation_coef, '-e', Distance_Euclidienne,'-t',self.right_sulci_label_to_sulci_name.fullPath()  ,'-p',Distance_au_plan )
     else:
          if self.Side in ('Left','Both'):
               context.system('siMeshSulciProjection', '-i', self.left_white_mesh.fullPath() ,  '-g' , self.Lgraph.fullPath(), '-l' , self.translation.fullPath() , '-m', self.gyri_model.fullPath() , '-s' , self.sulcus_identification, '-v' , self.mri_corrected.fullPath(), '-o' ,self.left_white_sulci.fullPath() ,'-V' , Volume_Closing,'-M',Mesh_Closing,'-n',Min_points_in_connected_component,'-a', Affine_estimation_coef, '-e', Distance_Euclidienne,'-t',self.left_sulci_label_to_sulci_name.fullPath() ,'--connexity','-p',Distance_au_plan )
          if self.Side in ('Right','Both'):
               context.system('siMeshSulciProjection', '-i', self.right_white_mesh.fullPath() ,  '-g' , self.Rgraph.fullPath(), '-l' , self.translation.fullPath() , '-m', self.gyri_model.fullPath() , '-s' , self.sulcus_identification, '-v' , self.mri_corrected.fullPath(), '-o' ,self.right_white_sulci.fullPath() ,'-V', Volume_Closing,'-M',Mesh_Closing,'-n',Min_points_in_connected_component,'-a', Affine_estimation_coef, '-e', Distance_Euclidienne,'-t', self.right_sulci_label_to_sulci_name.fullPath() ,'--connexity','-p',Distance_au_plan )
