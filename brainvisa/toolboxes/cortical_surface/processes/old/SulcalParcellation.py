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

name = 'Sulcal Parcellation'

userLevel = 2

signature = Signature(
    'Side', Choice("Both","Left","Right"),
    'Projection',Choice("Yes","No"),
    'Parcellation',Choice("White Matter Surface (2D)","Cortical Ribbon (3D)"),
    'Lgraph', ReadDiskItem( 'Cortical folds graph', 'Graph',requiredAttributes={ 'side': 'left' } ),
    'Rgraph', ReadDiskItem( 'Cortical folds graph', 'Graph',requiredAttributes={ 'side': 'right' } ),
    'mri_corrected', ReadDiskItem( 'T1 MRI Bias Corrected', 'GIS image' ),
    'sulcus_identification',Choice('label','name'),
    'translation',ReadDiskItem('Label Translation','Label Translation' ),
    'gyri_model',ReadDiskItem('Gyri Model','Gyri Model' ),
    'left_sulci_label_to_sulci_name',WriteDiskItem( 'Sulci To White Texture Translation', 'Text File',requiredAttributes={ 'side': 'left' }),
    'right_sulci_label_to_sulci_name',WriteDiskItem( 'Sulci To White Texture Translation', 'Text File',requiredAttributes={ 'side': 'right' }),
    'left_white_mesh',ReadDiskItem( 'Left Hemisphere White Mesh' ,
                                    shfjGlobals.aimsMeshFormats),
    'right_white_mesh',ReadDiskItem( 'Right Hemisphere White Mesh' ,
                                     shfjGlobals.aimsMeshFormats),
    'left_hemi_mesh', WriteDiskItem( 'Left Hemisphere Mesh', 'MESH mesh' ),
    'right_hemi_mesh', WriteDiskItem( 'Right Hemisphere Mesh', 'MESH mesh' ),
    'left_white_sulci',WriteDiskItem( 'Sulci White Texture' ,'Texture',
                                      requiredAttributes={ 'side': 'left' } ),
    'right_white_sulci',WriteDiskItem( 'Sulci White Texture' ,'Texture',
                                       requiredAttributes= \
                                       { 'side': 'right' } ),
    'left_grey_white', ReadDiskItem( 'Left Grey White Mask', 'GIS Image' ),
    'right_grey_white', ReadDiskItem( 'Right Grey White Mask', 'GIS Image' ),
    'left_sulcal_patch_texture',WriteDiskItem( 'Sulci White Texture Patch' ,'Texture', requiredAttributes={ 'side': 'left' } ),
    'right_sulcal_patch_texture',WriteDiskItem( 'Sulci White Texture Patch' ,'Texture', requiredAttributes={ 'side': 'right' } ),
    'left_sulcal_patch_graph',WriteDiskItem( 'Sulcal Patch Graph' ,'Graph',requiredAttributes= { 'side': 'left' } ),
    'right_sulcal_patch_graph',WriteDiskItem( 'Sulcal Patch Graph' ,'Graph',requiredAttributes= { 'side': 'right' } ),
    'left_sulcal_patch_volume',WriteDiskItem( 'Sulci White Volume Patch' ,'GIS Image',requiredAttributes={ 'side': 'left' } ),
    'right_sulcal_patch_volume',WriteDiskItem( 'Sulci White Volume Patch' ,'GIS Image',  requiredAttributes= { 'side': 'right' } )
)

def initialization( self ):
     self.linkParameters( 'left_white_mesh', 'Lgraph' )
     self.linkParameters( 'left_hemi_mesh', 'Lgraph' )
     self.linkParameters( 'Rgraph', 'Lgraph' )
     self.linkParameters( 'mri_corrected', 'Lgraph' )
     self.linkParameters( 'right_white_mesh', 'Rgraph' )
     self.linkParameters( 'right_hemi_mesh', 'Rgraph' )
     self.linkParameters( 'left_sulcal_patch_texture', 'Lgraph' )
     self.linkParameters( 'right_sulcal_patch_texture', 'Rgraph' )
     self.linkParameters( 'left_sulcal_patch_graph', 'Lgraph' )
     self.linkParameters( 'right_sulcal_patch_graph', 'Rgraph' )
     self.linkParameters( 'left_white_sulci', 'Lgraph' )
     self.linkParameters( 'right_white_sulci', 'Rgraph' )
     self.linkParameters( 'left_sulcal_patch_volume', 'Lgraph' )
     self.linkParameters( 'right_sulcal_patch_volume', 'Rgraph' )
     self.linkParameters( 'left_grey_white', 'Lgraph' )
     self.linkParameters( 'right_grey_white', 'Rgraph' )
     self.linkParameters( 'left_sulci_label_to_sulci_name', 'mri_corrected' )
     self.linkParameters( 'right_sulci_label_to_sulci_name', 'mri_corrected' )
     self.sulcus_identification = 'label'
     self.setOptional('right_grey_white','left_grey_white', 
                      'left_sulcal_patch_volume','right_sulcal_patch_volume', 
                      'Rgraph', 'Lgraph',
                      'left_white_mesh','right_white_mesh',
                      'left_sulcal_patch_volume','right_sulcal_patch_volume',
                      'left_hemi_mesh','right_hemi_mesh', 
                      'left_sulcal_patch_graph','right_sulcal_patch_graph',
                      'left_sulcal_patch_texture','right_sulcal_patch_texture',
                      'left_white_sulci','right_white_sulci')
     self.translation = '/home/appli/shared-main/nomenclature/translation/gyri.trl'
     self.gyri_model = '/home/appli/shared-main/models/gyrus/gyri.gyr'
     self.Projection = 'Yes'

def execution( self, context ): 

     call_list = ['siSulcalParcellation',
                  '-m', self.gyri_model.fullPath() ]
          
     if self.Parcellation == 'Cortical Ribbon (3D)':
          call_list += ['--3D']
         
     if self.Side in ('Left','Both'):
          if ( self.Projection == 'Yes' ):
               context.runProcess( 'CreateLabelTexture', 
                                   Side = 'Left',
                                   Lgraph = self.Lgraph,
                                   Rgraph = self.Rgraph,
                                   left_white_mesh = self.left_white_mesh,
                                   right_white_mesh = self.right_white_mesh,
                                   left_white_sulci = self.left_white_sulci,
                                   right_white_sulci = self.right_white_sulci,
                                   translation = self.translation,
                                   left_sulci_label_to_sulci_name = self.left_sulci_label_to_sulci_name,
                                   right_sulci_label_to_sulci_name = self.right_sulci_label_to_sulci_name,
                                   gyri_model = self.gyri_model,
                                   mri_corrected = self.mri_corrected,
                                   sulcus_identification = \
                                   self.sulcus_identification)
               
          io = ['-i', self.left_white_mesh.fullPath() ,
                '-s' ,self.left_white_sulci.fullPath() ,
                '-o' ,self.left_sulcal_patch_texture.fullPath() ,
                '-g',self.left_sulcal_patch_graph.fullPath(),
                '-p',self.left_sulcal_patch_volume.fullPath(),
                '-b', self.left_hemi_mesh.fullPath(),
                '--sulcitraduction', self.left_sulci_label_to_sulci_name.fullPath() ]
     
          
          if self.Parcellation == 'Cortical Ribbon (3D)':
               io += ['-V',self.left_grey_white.fullPath()]
          
          apply( context.system, call_list+io )
          

     if self.Side in ('Right','Both'):
          if ( self.Projection == 'Yes' ):
               context.runProcess( 'CreateLabelTexture', 
                                   Side = 'Right',
                                   Lgraph = self.Lgraph,
                                   Rgraph = self.Rgraph,
                                   left_white_mesh = self.left_white_mesh,
                                   right_white_mesh = self.right_white_mesh,
                                   left_white_sulci = self.left_white_sulci,
                                   right_white_sulci = self.right_white_sulci,
                                   translation = self.translation,
                                   left_sulci_label_to_sulci_name = self.left_sulci_label_to_sulci_name,
                                   right_sulci_label_to_sulci_name = self.right_sulci_label_to_sulci_name,
                                   gyri_model = self.gyri_model,
                                   mri_corrected = self.mri_corrected,
                                   sulcus_identification = \
                                   self.sulcus_identification)
          io = ['-i', self.right_white_mesh.fullPath() ,
                "-s" ,self.right_white_sulci.fullPath() ,
                '-o' ,self.right_sulcal_patch_texture.fullPath() ,
                '-g',self.right_sulcal_patch_graph.fullPath(),
                '-p',self.right_sulcal_patch_volume.fullPath(),
                '-b', self.right_hemi_mesh.fullPath(),
                '--sulcitraduction', self.right_sulci_label_to_sulci_name.fullPath()  ]

          if self.Parcellation == 'Cortical Ribbon (3D)':
               io += ['-V',self.right_grey_white.fullPath()]
               
          apply( context.system, call_list+io )
    
