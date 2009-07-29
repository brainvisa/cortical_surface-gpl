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

name = '2D Geodesic Primal Sketch'

userLevel = 2

signature = Signature(
     'Side', Choice("Both","Left","Right"),
     'Graph', Choice("No","Yes"),
     'mri_corrected', ReadDiskItem( 'T1 MRI Bias Corrected', 'GIS image' ),
     'left_white_mesh',ReadDiskItem( 'Left Hemisphere White Mesh' , shfjGlobals.aimsMeshFormats),
     'right_white_mesh',ReadDiskItem( 'Right Hemisphere White Mesh' , shfjGlobals.aimsMeshFormats),
     'left_white_mesh_inflated',ReadDiskItem( 'Inflated Hemisphere White Mesh' , shfjGlobals.aimsMeshFormats, requiredAttributes={ 'side': 'left' }),
     'right_white_mesh_inflated',ReadDiskItem( 'Inflated Hemisphere White Mesh' , shfjGlobals.aimsMeshFormats, requiredAttributes={ 'side': 'right' }),
     'left_white_curvature', ReadDiskItem( 'White Curvature Texture', 'Texture', requiredAttributes={ 'side': 'left' } ), 
     'right_white_curvature', ReadDiskItem( 'White Curvature Texture', 'Texture', requiredAttributes={ 'side': 'right' } ),
     'left_white_curvature_blobs', WriteDiskItem( 'Blob White Curvature Texture', 'Texture', requiredAttributes={ 'side': 'left' } ), 
     'right_white_curvature_blobs', WriteDiskItem( 'Blob White Curvature Texture', 'Texture', requiredAttributes={ 'side': 'right' } ),
     'left_white_curvature_ss', WriteDiskItem( 'Scale Space White Curvature Texture', 'Texture', requiredAttributes={ 'side': 'left' } ), 
     'right_white_curvature_ss', WriteDiskItem( 'Scale Space White Curvature Texture', 'Texture', requiredAttributes={ 'side': 'right' } ),
     'left_white_primal', WriteDiskItem( 'Primal Sketch', 'Graph', requiredAttributes={ 'side': 'left' } ), 
     'right_white_primal', WriteDiskItem( 'Primal Sketch', 'Graph', requiredAttributes={ 'side': 'right' } ),
     'left_white_glb', WriteDiskItem( 'Grey Level Blob Graph', 'Graph', requiredAttributes={ 'side': 'left' } ), 
     'right_white_glb', WriteDiskItem( 'Grey Level Blob Graph', 'Graph', requiredAttributes={ 'side': 'right' } ),
     'Begin_scale', Float(),
     'End_scale', Float(),
     'Scale_space_sample_step', Float(),
     'Recursivity', Integer(),
     'Laplacian_threshold', Float(),
     'Texture_threshold', Float(),
     'Smooth_step', Float(),
     'Grey_level_blob', Choice("Minima","Maxima"),
     'Scale', Choice("Logarithmic","Normal"),
     'Surface_type', Choice("tore","surface"),
     'Grow_mode', Choice("scale","translate","pushnormal"),
     'Coef_translation', Float()
)

def initialization( self ):
     self.linkParameters( 'left_white_mesh', 'mri_corrected' )
     self.linkParameters( 'right_white_mesh', 'mri_corrected' )
     self.linkParameters( 'left_white_mesh_inflated', 'left_white_mesh' )
     self.linkParameters( 'right_white_mesh_inflated', 'right_white_mesh' )
     self.linkParameters( 'left_white_curvature', 'left_white_mesh' )
     self.linkParameters( 'right_white_curvature', 'right_white_mesh' )
     self.linkParameters( 'right_white_curvature_blobs', 'right_white_mesh' )
     self.linkParameters( 'left_white_curvature_blobs', 'left_white_mesh' )
     self.linkParameters( 'right_white_curvature_ss', 'right_white_mesh' )
     self.linkParameters( 'left_white_curvature_ss', 'left_white_mesh' )
     self.linkParameters( 'right_white_primal', 'right_white_mesh' )
     self.linkParameters( 'left_white_primal', 'left_white_mesh' )
     self.linkParameters( 'right_white_glb', 'right_white_mesh' )
     self.linkParameters( 'left_white_glb', 'left_white_mesh' )
     self.Recursivity = 1
     self.Texture_threshold = 1e10 
     self.Laplacian_threshold = 0.95
     self.Smooth_step = 0.02
     self.Coef_translation = 0.5
     self.Begin_scale = 10
     self.End_scale = 800
     self.Scale_space_sample_step = 1.19
     if self.Grey_level_blob == "Minima":
          self.convexity = -1
     else:
           self.convexity = 1
     

def execution( self, context ): 

     if self.Graph == "Yes": 
          if self.Scale == "Normal":
               if self.Side in ('Left','Both'):
                    context.system('AimsTexture2Primal', '-i', self.left_white_curvature.fullPath() , '-v', self.mri_corrected.fullPath(), '-M', self.left_white_mesh.fullPath(), '-b',self.Begin_scale,'-e',self.End_scale,'-d', self.Scale_space_sample_step, '-S', self.Smooth_step,'-m', self.Recursivity,'-H',self.Texture_threshold,'-W', self.Laplacian_threshold, '-s',self.Surface_type, '-G', self.Grow_mode, '-t',self.Coef_translation, '-c', self.convexity, '-o', self.left_white_curvature_ss.fullPath(), '-B',self.left_white_curvature_blobs.fullPath(),'-p', self.left_white_primal.fullPath(),'-g',self.left_white_glb.fullPath(),'-F',self.left_white_mesh_inflated.fullPath() )
               if self.Side in ('Right','Both'):
                    context.system('AimsTexture2Primal', '-i', self.right_white_curvature.fullPath() , '-v', self.mri_corrected.fullPath(), '-M', self.right_white_mesh.fullPath(), '-b',self.Begin_scale,'-e',self.End_scale,'-d', self.Scale_space_sample_step, '-S', self.Smooth_step,'-m', self.Recursivity,'-H',self.Texture_threshold,'-W', self.Laplacian_threshold, '-s',self.Surface_type, '-G', self.Grow_mode, '-t',self.Coef_translation, '-c', self.convexity, '-o', self.right_white_curvature_ss.fullPath(), '-B',self.right_white_curvature_blobs.fullPath(),'-p', self.right_white_primal.fullPath(),'-g',self.right_white_glb.fullPath(),'-F',self.right_white_mesh_inflated.fullPath())
          else:
               if self.Side in ('Left','Both'):
                    context.system('AimsTexture2Primal', '-i', self.left_white_curvature.fullPath() , '-v', self.mri_corrected.fullPath(), '-M', self.left_white_mesh.fullPath(), '-b',self.Begin_scale,'-e',self.End_scale,'-d', self.Scale_space_sample_step, '-S', self.Smooth_step,'-m', self.Recursivity,'-H',self.Texture_threshold,'-W', self.Laplacian_threshold, '-s',self.Surface_type, '-G', self.Grow_mode, '-t',self.Coef_translation, '-c', self.convexity, '-o', self.left_white_curvature_ss.fullPath(), '-B',self.left_white_curvature_blobs.fullPath(),'--log','-p', self.left_white_primal.fullPath(),'-g',self.left_white_glb.fullPath(),'-F',self.left_white_mesh_inflated.fullPath() )
               if self.Side in ('Right','Both'):
                    context.system('AimsTexture2Primal', '-i', self.right_white_curvature.fullPath() , '-v', self.mri_corrected.fullPath(), '-M', self.right_white_mesh.fullPath(), '-b',self.Begin_scale,'-e',self.End_scale,'-d', self.Scale_space_sample_step, '-S', self.Smooth_step,'-m', self.Recursivity,'-H',self.Texture_threshold,'-W', self.Laplacian_threshold, '-s',self.Surface_type, '-G', self.Grow_mode, '-t',self.Coef_translation, '-c', self.convexity, '-o', self.right_white_curvature_ss.fullPath(), '-B',self.right_white_curvature_blobs.fullPath(), '--log','-p', self.right_white_primal.fullPath(),'-g',self.right_white_glb.fullPath(),'-F',self.right_white_mesh_inflated.fullPath())
     else:
        if self.Scale == "Normal":
               if self.Side in ('Left','Both'):
                    context.system('AimsTexture2Primal', '-i', self.left_white_curvature.fullPath() , '-v', self.mri_corrected.fullPath(), '-M', self.left_white_mesh.fullPath(), '-b',self.Begin_scale,'-e',self.End_scale,'-d', self.Scale_space_sample_step, '-S', self.Smooth_step,'-m', self.Recursivity,'-H',self.Texture_threshold,'-W', self.Laplacian_threshold, '-s',self.Surface_type, '-G', self.Grow_mode, '-t',self.Coef_translation, '-c', self.convexity, '-o', self.left_white_curvature_ss.fullPath(), '-B',self.left_white_curvature_blobs.fullPath(),'-F',self.left_white_mesh_inflated.fullPath() )
               if self.Side in ('Right','Both'):
                    context.system('AimsTexture2Primal', '-i', self.right_white_curvature.fullPath() , '-v', self.mri_corrected.fullPath(), '-M', self.right_white_mesh.fullPath(), '-b',self.Begin_scale,'-e',self.End_scale,'-d', self.Scale_space_sample_step, '-S', self.Smooth_step,'-m', self.Recursivity,'-H',self.Texture_threshold,'-W', self.Laplacian_threshold, '-s',self.Surface_type, '-G', self.Grow_mode, '-t',self.Coef_translation, '-c', self.convexity, '-o', self.right_white_curvature_ss.fullPath(), '-B',self.right_white_curvature_blobs.fullPath(),'-F',self.right_white_mesh_inflated.fullPath())
        else:
             if self.Side in ('Left','Both'):
                  context.system('AimsTexture2Primal', '-i', self.left_white_curvature.fullPath() , '-v', self.mri_corrected.fullPath(), '-M', self.left_white_mesh.fullPath(), '-b',self.Begin_scale,'-e',self.End_scale,'-d', self.Scale_space_sample_step, '-S', self.Smooth_step,'-m', self.Recursivity,'-H',self.Texture_threshold,'-W', self.Laplacian_threshold, '-s',self.Surface_type, '-G', self.Grow_mode, '-t',self.Coef_translation, '-c', self.convexity, '-o', self.left_white_curvature_ss.fullPath(), '-B',self.left_white_curvature_blobs.fullPath(),'--log','-F',self.left_white_mesh_inflated.fullPath() )
             if self.Side in ('Right','Both'):
                  context.system('AimsTexture2Primal', '-i', self.right_white_curvature.fullPath() , '-v', self.mri_corrected.fullPath(), '-M', self.right_white_mesh.fullPath(), '-b',self.Begin_scale,'-e',self.End_scale,'-d', self.Scale_space_sample_step, '-S', self.Smooth_step,'-m', self.Recursivity,'-H',self.Texture_threshold,'-W', self.Laplacian_threshold, '-s',self.Surface_type, '-G', self.Grow_mode, '-t',self.Coef_translation, '-c', self.convexity, '-o', self.right_white_curvature_ss.fullPath(), '-B',self.right_white_curvature_blobs.fullPath(), '--log','-F',self.right_white_mesh_inflated.fullPath())  
