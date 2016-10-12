# -*- coding: utf-8 -*-
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
from brainvisa.data.neuroHierarchy import databases

name = 'Create Label Texture'

userLevel = 1

signature = Signature(
    'Side', Choice("Both","Left","Right"),
    'Lgraph', ReadDiskItem( 'Cortical folds graph', 'Graph',requiredAttributes={ 'side': 'left' } ),
    'Rgraph', ReadDiskItem( 'Cortical folds graph', 'Graph',requiredAttributes={ 'side': 'right' } ),
    'sulcus_identification',Choice('name','label'),
    'mri_corrected', ReadDiskItem( 'T1 MRI Bias Corrected', 'Aims readable volume formats' ),
    'left_white_mesh',ReadDiskItem( 'Left Hemisphere White Mesh' , 'aims Mesh Formats'),
    'right_white_mesh',ReadDiskItem( 'Right Hemisphere White Mesh' , 'aims Mesh Formats'),
    'MNI', Choice("No","Yes"),
    'mni_mesh', ReadDiskItem('MNI Cortex Mesh', 'Aims mesh formats'),
    'left_white_sulci',WriteDiskItem( 'Sulci White Texture' ,'aims texture formats',requiredAttributes={ 'side': 'left' } ),
    'right_white_sulci',WriteDiskItem( 'Sulci White Texture' ,'aims texture formats',requiredAttributes={ 'side': 'right' } ),
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
     self.linkParameters( 'left_white_sulci', 'Lgraph' )
     self.linkParameters( 'right_white_sulci', 'Rgraph' )
     self.linkParameters( 'mni_mesh', 'mri_corrected' )
     self.linkParameters( 'left_sulci_label_to_sulci_name', 'Lgraph' )
     self.linkParameters( 'right_sulci_label_to_sulci_name', 'Rgraph' )
     self.setOptional('mni_mesh','Rgraph', 'Lgraph', 'left_white_mesh', 'right_white_mesh', 'left_white_sulci', 'right_white_sulci'  )
     self.sulcus_identification = 'label'
     self.translation = ReadDiskItem('Label Translation','Label Translation' ).findValue( { 'filename_variable' : 'gyri' } )
     self.gyri_model = databases.getDiskItemFromUuid( '172c4168-a9d3-dc41-464c-1226ad07c19c' )

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
          lmesh = 'mni_mesh'
          rmesh = 'mni_mesh'
     else:
          lmesh = 'left_white_mesh'
          rmesh = 'right_white_mesh'

     if self.Metric == 'Euclidean':
          if self.Side in ('Left','Both'):
               if not self.Lgraph or not self.left_white_mesh or not self.left_white_sulci:
                 raise ValidationError( 'Lgraph, ' + lmesh + ' and left_white_sulci are required to process the right hemisphere' )
               context.system('siMeshSulciProjection', '-i', self.left_white_mesh.fullPath() ,  '-g' , self.Lgraph.fullPath(), '-l' , self.translation.fullPath() , '-m', self.gyri_model.fullPath() , '-s' , self.sulcus_identification, '-v' , self.mri_corrected.fullPath(), '-o' ,self.left_white_sulci.fullPath() ,'-V' , Volume_Closing,'-M',Mesh_Closing,'-n',Min_points_in_connected_component,'-a', Affine_estimation_coef, '-e', Distance_Euclidienne,'-t', self.left_sulci_label_to_sulci_name.fullPath() ,'-p',Distance_au_plan )
          if self.Side in ('Right','Both'):
               if not self.Rgraph or not self.right_white_mesh or not self.right_white_sulci:
                 raise ValidationError( 'Rgraph, ' + rmesh + ' and right_white_sulci are required to process the right hemisphere' )
               context.system('siMeshSulciProjection', '-i', self.right_white_mesh.fullPath() ,  '-g' , self.Rgraph.fullPath(), '-l' , self.translation.fullPath() , '-m', self.gyri_model.fullPath() , '-s' , self.sulcus_identification, '-v' , self.mri_corrected.fullPath(), '-o' ,self.right_white_sulci.fullPath() ,'-V', Volume_Closing,'-M',Mesh_Closing,'-n',Min_points_in_connected_component,'-a', Affine_estimation_coef, '-e', Distance_Euclidienne,'-t',self.right_sulci_label_to_sulci_name.fullPath()  ,'-p',Distance_au_plan )
     else:
          if self.Side in ('Left','Both'):
               if not self.Lgraph or not self.left_white_mesh or not self.left_white_sulci:
                 raise ValidationError( 'Lgraph, ' + lmesh + ' and left_white_sulci are required to process the right hemisphere' )
               context.system('siMeshSulciProjection', '-i', self.left_white_mesh.fullPath() ,  '-g' , self.Lgraph.fullPath(), '-l' , self.translation.fullPath() , '-m', self.gyri_model.fullPath() , '-s' , self.sulcus_identification, '-v' , self.mri_corrected.fullPath(), '-o' ,self.left_white_sulci.fullPath() ,'-V' , Volume_Closing,'-M',Mesh_Closing,'-n',Min_points_in_connected_component,'-a', Affine_estimation_coef, '-e', Distance_Euclidienne,'-t',self.left_sulci_label_to_sulci_name.fullPath() ,'--connexity','-p',Distance_au_plan )
          if self.Side in ('Right','Both'):
               if not self.Rgraph or not self.right_white_mesh or not self.right_white_sulci:
                 raise ValidationError( 'Rgraph, ' + rmesh + ' and right_white_sulci are required to process the right hemisphere' )
               context.system('siMeshSulciProjection', '-i', self.right_white_mesh.fullPath() ,  '-g' , self.Rgraph.fullPath(), '-l' , self.translation.fullPath() , '-m', self.gyri_model.fullPath() , '-s' , self.sulcus_identification, '-v' , self.mri_corrected.fullPath(), '-o' ,self.right_white_sulci.fullPath() ,'-V', Volume_Closing,'-M',Mesh_Closing,'-n',Min_points_in_connected_component,'-a', Affine_estimation_coef, '-e', Distance_Euclidienne,'-t', self.right_sulci_label_to_sulci_name.fullPath() ,'--connexity','-p',Distance_au_plan )
