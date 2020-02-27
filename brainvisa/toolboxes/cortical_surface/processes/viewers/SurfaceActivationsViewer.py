# -*- coding: iso-8859-1 -*-
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

# By A.M. - 02-02-2011
# Inpired from inflateCorticalSurface.py, AnatomistShowInflatedWhiteMesh.py

from __future__ import absolute_import
from brainvisa.processes import *
from brainvisa import anatomist

name = 'Surface MultiActivations Viewer'
roles = ('viewer',)
userLevel = 0

################################################################################

# Argument declaration
signature = Signature(
  'activations_texture_1_red',ReadDiskItem( 'Functional Time Texture',
    'Aims texture formats' ),
  'activations_texture_2_green',ReadDiskItem( 'Functional Time Texture',
    'Aims texture formats' ),
  'activations_texture_3_blue',ReadDiskItem( 'Functional Time Texture',
    'Aims texture formats' ),
  'input_mesh',ReadDiskItem( 'Hemisphere White Mesh' , 'Aims mesh formats' ),
  'compute_inflated_mesh', Boolean(),
  'inflated_mesh',WriteDiskItem( 'Inflated Hemisphere White Mesh',
    'Aims mesh formats' ),
  'curvature_texture',WriteDiskItem( 'White Curvature Texture',
    'Aims texture formats' ),
  'iterations', Integer(),
  'save_inflate_sequence', Boolean(),
)

# Default values
def initialization( self ):
  self.linkParameters( 'input_mesh', 'activations_texture_1_red' )
  self.linkParameters( 'input_mesh', 'activations_texture_2_green' )
  self.linkParameters( 'input_mesh', 'activations_texture_3_blue' )
  self.linkParameters( 'inflated_mesh', 'input_mesh' )
  self.linkParameters( 'curvature_texture', 'input_mesh' )
  self.compute_inflated_mesh = 1
  self.iterations = 500
  self.save_inflate_sequence = 0
  self.setOptional( 'activations_texture_2_green','activations_texture_3_blue',
                    'inflated_mesh', 'curvature_texture' )


# Surface activations viewer process
#

def execution( self, context ):

  # Quite standard values (not available for manipulation on the interface for easier use)
  normal_force = 0.01
  spring_force = 0.01
  smoothing_force = 0.5

  if os.path.exists(self.inflated_mesh.fullName() + '.loc'):
    context.write( "Inflated cortical surface locked")
  elif self.compute_inflated_mesh and self.inflated_mesh is not None and self.curvature_texture is not None:
  ##############################################################################
  # If the user asks to compute the mesh or if it does not exist ###############
  ##############################################################################
    if self.save_inflate_sequence:
      context.system( 'AimsInflate', '-i', self.input_mesh.fullPath(), 
                                     '-o', self.inflated_mesh.fullPath(), 
                                     '-t', self.iterations, 
                                     '-Kn', normal_force, 
                                     '-Ksp', spring_force, 
                                     '-Ksm', smoothing_force, 
                                     '-c', self.curvature_texture.fullPath(), 
                                     '-S')
    else:
      context.system( 'AimsInflate', '-i', self.input_mesh.fullPath(), 
                                     '-o', self.inflated_mesh.fullPath(), 
                                     '-t', self.iterations, 
                                     '-Kn', normal_force, 
                                     '-Ksp', spring_force, 
                                     '-Ksm', smoothing_force,
                                     '-c', self.curvature_texture.fullPath())

  if self.inflated_mesh is not None:
      mesh = self.inflated_mesh
  else:
      mesh = self.input_mesh

  # Final visualization
  #
  
  a = anatomist.Anatomist()
  if (self.curvature_texture is not None and self.activations_texture_1_red is not None):

    # We get the different palettes...
    palette_curvature = a.getPalette("B-W LINEAR") #or "actif_ret"
    palette_activations_1 = a.getPalette("RED TEMPERATURE") 
    palette_activations_2 = a.getPalette("Green-White-linear") # or "Green-White-linear-fusion"
    palette_activations_3 = a.getPalette("Blue-White") # or "Blue-White-fusion"

    mymesh = a.loadObject( mesh )

    if palette_curvature is not None:
      duplicate = True
    tex_curvature = a.loadObject( self.curvature_texture, duplicate=duplicate ) 
    if palette_activations_1 is not None:
      duplicate = True
    tex_activations_1 = a.loadObject( self.activations_texture_1_red, duplicate=duplicate ) 

    # Setting textures palettes
    tex_curvature.setPalette( palette_curvature )
    tex_activations_1.setPalette( palette_activations_1 )
	# It is possible to impose the min and max values:
    #minValcurvature = 0.4
    #maxValcurvature = 0.6
    #minValactivations = 0.5
    #maxValactivations = 1
    #tex_curvature.setPalette( palette_curvature, minVal = minValcurvature, maxVal = maxValcurvature )
    #tex_activations_1.setPalette( palette_activations_1, minVal = minValactivations, maxVal = maxValactivations )
    # or
    #a.setObjectPalette([tex_curvature], palette_curvature, minVal=minValcurvature, maxVal=maxValcurvature) 
    #a.setObjectPalette([tex_activations_1], palette_activations_1, minVal=minValactivations, maxVal=maxValactivations) 

    # Setting activations 2 and 3 if present
    if self.activations_texture_2_green is not None:
        if palette_activations_2 is not None:
          duplicate = True
        tex_activations_2 = a.loadObject( self.activations_texture_2_green, duplicate=duplicate ) 
        tex_activations_2.setPalette( palette_activations_2 )

    if self.activations_texture_3_blue is not None:
        if palette_activations_3 is not None:
          duplicate = True
        tex_activations_3 = a.loadObject( self.activations_texture_3_blue, duplicate=duplicate ) 
        tex_activations_3.setPalette( palette_activations_3 )

    # Multitexture fusion of the curvature texture with the 1, 2 or 3 activation textures
    # Help found here:
    # http://brainvisa.info/doc/anatomist/html/fr/programmation/commands.html#TexturingParams
    # and here:
    # http://www.brainvisa.info/doc/brainvisa-4.0/epydoc/brainvisa.anatomist-pysrc.html
    if (self.activations_texture_2_green is not None) and (self.activations_texture_3_blue is not None):
        fusionMultiTextureActivations = a.fusionObjects( [tex_curvature, tex_activations_1, tex_activations_2, tex_activations_3], method='FusionMultiTextureMethod' )
#        a.execute( "FusionMultiTextureMethod", object=fusionMultiTextureActivations, texture=1, mode="add", rate = 1.0 )
        a.execute("TexturingParams", objects=[fusionMultiTextureActivations], texture_index = 3, mode = "add", rate = 1) 
        a.execute("TexturingParams", objects=[fusionMultiTextureActivations], texture_index = 2, mode = "add", rate = 1) 
        a.execute("TexturingParams", objects=[fusionMultiTextureActivations], texture_index = 1, mode = "add", rate = 1) 
        a.execute("TexturingParams", objects=[fusionMultiTextureActivations], texture_index = 0, mode = "geometric", rate = 1) 
    elif self.activations_texture_2_green is not None:
        fusionMultiTextureActivations = a.fusionObjects( [tex_curvature, tex_activations_1, tex_activations_2], method='FusionMultiTextureMethod' )
#        a.execute( "FusionMultiTextureMethod", object=fusionMultiTextureActivations, texture=1, mode="add", rate = 1.0 )
        a.execute("TexturingParams", objects=[fusionMultiTextureActivations], texture_index = 2, mode = "add", rate = 1) 
        a.execute("TexturingParams", objects=[fusionMultiTextureActivations], texture_index = 1, mode = "add", rate = 1) 
        a.execute("TexturingParams", objects=[fusionMultiTextureActivations], texture_index = 0, mode = "geometric", rate = 1) 
    elif self.activations_texture_3_blue is not None:
        fusionMultiTextureActivations = a.fusionObjects( [tex_curvature, tex_activations_1, tex_activations_3], method='FusionMultiTextureMethod' )
#        a.execute( "FusionMultiTextureMethod", object=fusionMultiTextureActivations, texture=1, mode="add", rate = 1.0 )
        a.execute("TexturingParams", objects=[fusionMultiTextureActivations], texture_index = 2, mode = "add", rate = 1) 
        a.execute("TexturingParams", objects=[fusionMultiTextureActivations], texture_index = 1, mode = "add", rate = 1) 
        a.execute("TexturingParams", objects=[fusionMultiTextureActivations], texture_index = 0, mode = "geometric", rate = 1) 
    else:
        fusionMultiTextureActivations = a.fusionObjects( [tex_curvature, tex_activations_1], method='FusionMultiTextureMethod' )
#        a.execute( "FusionMultiTextureMethod", object=fusionMultiTextureActivations, mode="add", rate = 1.0 )
        a.execute("TexturingParams", objects=[fusionMultiTextureActivations], texture_index = 1, mode = "add", rate = 1) 
        a.execute("TexturingParams", objects=[fusionMultiTextureActivations], texture_index = 0, mode = "geometric", rate = 1) 

    # TextSurf fusion of the multitexture with the mesh
    fusionTexSurfFinal = a.fusionObjects( [mymesh, fusionMultiTextureActivations], method='FusionTexSurfMethod' )

    # Creation of a 3D window
    mywindow = a.createWindow( '3D' )
    a.execute("WindowConfig", windows=[mywindow], cursor_visibility=0 ) # Removes cursor

    # The fusion is added into the window
    mywindow.addObjects([fusionTexSurfFinal])

    # Display palette control window
    a.execute("PopupPalette", objects=[tex_curvature] )
    a.execute("PopupPalette", objects=[tex_activations_1] )
    if self.activations_texture_2_green is not None:
        a.execute("PopupPalette", objects=[tex_activations_2] )
    if self.activations_texture_3_blue is not None:
        a.execute("PopupPalette", objects=[tex_activations_3] )

    if (self.activations_texture_2_green is not None) and (self.activations_texture_3_blue is not None):
        self._dontdestroy = [ mymesh, tex_curvature, tex_activations_1, tex_activations_2, tex_activations_3, fusionMultiTextureActivations, fusionTexSurfFinal, mywindow ]
    elif self.activations_texture_2_green is not None:
        self._dontdestroy = [ mymesh, tex_curvature, tex_activations_1, tex_activations_2, fusionMultiTextureActivations, fusionTexSurfFinal, mywindow ]
    elif self.activations_texture_3_blue is not None:
        self._dontdestroy = [ mymesh, tex_curvature, tex_activations_1, tex_activations_3, fusionMultiTextureActivations, fusionTexSurfFinal, mywindow ]
    else: # if both activations 2 and 3 are not None
        self._dontdestroy = [ mymesh, tex_curvature, tex_activations_1, fusionMultiTextureActivations, fusionTexSurfFinal, mywindow ]

  else:
    if self.curvature_texture is not None:
      return a.viewTextureOnMesh( mesh, self.curvature_texture,
                                  a.getPalette('B-W LINEAR'))
    elif self.activations_texture_1_red is not None:
      return a.viewTextureOnMesh( mesh, self.activations_texture_1_red,
                                  a.getPalette('Rainbow1-fusion'))
    else:
      return a.viewMesh( mesh )

################################################################################
