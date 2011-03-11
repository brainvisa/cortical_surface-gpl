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

from neuroProcesses import *
from brainvisa import anatomist

name = 'Surface Activations Viewer'
userLevel = 0

################################################################################

# Argument declaration
signature = Signature(
  'activations_texture',WriteDiskItem( 'Functional Time Texture',
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
  self.linkParameters( 'input_mesh', 'activations_texture' )
  self.linkParameters( 'inflated_mesh', 'input_mesh' )
  self.linkParameters( 'curvature_texture', 'input_mesh' )
  self.compute_inflated_mesh = 1
  self.iterations = 500
  self.save_inflate_sequence = 0


# Surface activations viewer process
#

def execution( self, context ):

  # Quite standard values (not available for manipulation on the interface for easier use)
  normal_force = 0.01
  spring_force = 0.01
  smoothing_force = 0.5

  if os.path.exists(self.inflated_mesh.fullName() + '.loc'):
    context.write( "Inflated cortical surface locked")
  elif self.compute_inflated_mesh or not os.path.exists(self.inflated_mesh.fullName() + '.mesh'): 
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

  # Final visualization
  #
  
  a = anatomist.Anatomist()
  if (self.curvature_texture is not None and self.activations_texture is not None):

    # We get the different palettes...
    palette_curvature = a.getPalette("B-W LINEAR") 
    palette_activations = a.getPalette("Rainbow1-fusion") 

    mymesh = a.loadObject( self.inflated_mesh ) 

    if palette_curvature is not None:
      duplicate = True
    tex_curvature = a.loadObject( self.curvature_texture, duplicate=duplicate ) 
    if palette_activations is not None:
      duplicate = True
    tex_activations = a.loadObject( self.activations_texture, duplicate=duplicate ) 

    # Setting textures palettes
    tex_curvature.setPalette( palette_curvature )
    tex_activations.setPalette( palette_activations )
	# It is possible to impose the min and max values:
    #minValcurvature = 0.4
    #maxValcurvature = 0.6
    #minValactivations = 0.5
    #maxValactivations = 1
    #tex_curvature.setPalette( palette_curvature, minVal = minValcurvature, maxVal = maxValcurvature )
    #tex_activations.setPalette( palette_activations, minVal = minValactivations, maxVal = maxValactivations )
    # or
    #a.setObjectPalette([tex_curvature], palette_curvature, minVal=minValcurvature, maxVal=maxValcurvature) 
    #a.setObjectPalette([tex_activations], palette_activations, minVal=minValactivations, maxVal=maxValactivations) 

    # Multitexture fusion of the two textures
    fusionMultiTexture2textures = a.fusionObjects( [tex_curvature, tex_activations], method='FusionMultiTextureMethod' )

    # TextSurf fusion of the multetexture with the mesh
    fusionTexSurfFinal = a.fusionObjects( [mymesh, fusionMultiTexture2textures], method='FusionTexSurfMethod' )

    # Creation of a 3D window
    mywindow = a.createWindow( '3D' )
    a.execute("WindowConfig", windows=[mywindow], cursor_visibility=0 ) # Removes cursor

    # The fusion is added into the window
    mywindow.addObjects([fusionTexSurfFinal])

    # Display palette control window
    a.execute("PopupPalette", objects=tex_curvature )
    a.execute("PopupPalette", objects=tex_activations )

    self._dontdestroy = [ mymesh, tex_curvature, tex_activations, fusionMultiTexture2textures, fusionTexSurfFinal, mywindow ]

  else:
    if self.curvature_texture is not None:
      return a.viewTextureOnMesh( self.inflated_mesh, self.curvature_texture, 
                                  a.getPalette('B-W LINEAR'))
    elif self.activations_texture is not None:
      return a.viewTextureOnMesh( self.inflated_mesh, self.activations_texture, 
                                  a.getPalette('Rainbow1-fusion'))
    else:
      return a.viewMesh( self.inflated_mesh )

################################################################################
