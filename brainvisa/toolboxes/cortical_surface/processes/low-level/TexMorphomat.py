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

name = 'Texture Mathematic Morphology'

userLevel = 1

signature = Signature(
     'Texture', ReadDiskItem( 'Texture', 'Texture' ), 
     'Mesh', ReadDiskItem( 'Mesh', shfjGlobals.aimsMeshFormats ), 
     'output_texture',WriteDiskItem( 'Texture', 'Texture' ),
     'Mode', Choice("Dilation",
                    "Erosion",
                    "Opening",
                    "Closing",
                    "Voronoi_Diagram"),
     'Radius', Float(),
     'background_label', Integer(),
     'forbidden_label', Integer(),
     'Metric', Choice("Euclidean","Mesh_Connectivity"),
)

def initialization( self ):
     self.linkParameters( 'Mesh', 'Texture' )
     self.background_label = 0
     self.forbidden_label = -1
     self.setOptional('Radius')
def execution( self, context ):
    if (self.Metric == 'Euclidean'): 
         if (self.Mode == 'Dilation'):
              context.system('AimsTextureDilation','-i',self.Mesh.fullPath(),'-t',self.Texture.fullPath(),
                             '-o',self.output_texture.fullPath(),'-s',self.Radius,'-b',self.background_label,
                             '-f',self.forbidden_label)
         if (self.Mode == 'Erosion'):
              context.system('AimsTextureErosion','-i',self.Mesh.fullPath(),'-t',self.Texture.fullPath(),
                             '-o',self.output_texture.fullPath(),'-s',self.Radius,'-b',self.background_label,
                             '-f',self.forbidden_label)
         if (self.Mode == 'Voronoi_Diagram'):
              context.system('AimsTextureVoronoi','-m',self.Mesh.fullPath(),'-t',self.Texture.fullPath(),
                             '-o',self.output_texture.fullPath(),'-b',self.background_label,
                             '-f',self.forbidden_label)
         if (self.Mode == 'Closing'):
              tex_temp = context.temporary( 'Texture')
              context.system('AimsTextureDilation','-i',self.Mesh.fullPath(),'-t',self.Texture.fullPath(),
                             '-o',tex_temp.fullPath(),'-s',self.Radius,'-b',self.background_label,
                             '-f',self.forbidden_label)
              context.system('AimsTextureErosion','-i',self.Mesh.fullPath(),'-t',tex_temp.fullPath(),
                             '-o',self.output_texture.fullPath(),'-s',self.Radius,'-b',self.background_label,
                             '-f',self.forbidden_label)
         if (self.Mode == 'Opening'):
              tex_temp = context.temporary( 'Texture')
              context.system('AimsTextureErosion','-i',self.Mesh.fullPath(),'-t',self.Texture.fullPath(),
                             '-o',tex_temp.fullPath(),'-s',self.Radius,'-b',self.background_label,
                             '-f',self.forbidden_label)
              context.system('AimsTextureDilation','-i',self.Mesh.fullPath(),'-t',tex_temp.fullPath(),
                             '-o',self.output_texture.fullPath(),'-s',self.Radius,'-b',self.background_label,
                             '-f',self.forbidden_label)
    else:
         if (self.Mode == 'Dilation'):
              context.system('AimsTextureDilation','-i',self.Mesh.fullPath(),'-t',self.Texture.fullPath(),
                             '-o',self.output_texture.fullPath(),'-s',self.Radius,'-b',self.background_label,
                             '-f',self.forbidden_label,'--connexity')
         if (self.Mode == 'Erosion'):
              context.system('AimsTextureErosion','-i',self.Mesh.fullPath(),'-t',self.Texture.fullPath(),
                             '-o',self.output_texture.fullPath(),'-s',self.Radius,'-b',self.background_label,
                             '-f',self.forbidden_label,'--connexity')
         if (self.Mode == 'Voronoi_Diagram'):
              context.system('AimsTextureVoronoi','-m',self.Mesh.fullPath(),'-t',self.Texture.fullPath(),
                             '-o',self.output_texture.fullPath(),'-b',self.background_label,
                             '-f',self.forbidden_label,'--connexity')
         if (self.Mode == 'Closing'):
              tex_temp = context.temporary( 'Texture')
              context.system('AimsTextureDilation','-i',self.Mesh.fullPath(),'-t',self.Texture.fullPath(),
                             '-o',tex_temp.fullPath(),'-s',self.Radius,'-b',self.background_label,
                             '-f',self.forbidden_label,'--connexity')
              context.system('AimsTextureErosion','-i',self.Mesh.fullPath(),'-t',tex_temp.fullPath(),
                             '-o',self.output_texture.fullPath(),'-s',self.Radius,'-b',self.background_label,
                             '-f',self.forbidden_label,'--connexity')
         if (self.Mode == 'Opening'):
              tex_temp = context.temporary( 'Texture')
              context.system('AimsTextureErosion','-i',self.Mesh.fullPath(),'-t',self.Texture.fullPath(),
                             '-o',tex_temp.fullPath(),'-s',self.Radius,'-b',self.background_label,
                             '-f',self.forbidden_label,'--connexity')
              context.system('AimsTextureDilation','-i',self.Mesh.fullPath(),'-t',tex_temp.fullPath(),
                             '-o',self.output_texture.fullPath(),'-s',self.Radius,'-b',self.background_label,
                             '-f',self.forbidden_label,'--connexity')
