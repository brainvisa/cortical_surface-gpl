# -*- coding: utf-8 -*-
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


def validation():
  try:
    from brainvisa.cortical_surface.surface_tools import texture_tools
  except:
    raise ValidationError( 'brainvisa.cortical_surface.parameterization.surface_tools module can not be imported.' )
  
from brainvisa.cortical_surface.surface_tools import texture_tools as textureTls
from brainvisa.processes import *
import shfjGlobals  
from soma import aims
import numpy as np


#from brainvisa import anatomist

name = 'Insular Pole From Manual Drawing'

userLevel = 0

# def validation():
#     anatomist.validation()
    
signature = Signature(
    'input_texture',ReadDiskItem( 'Texture', 'Aims Texture formats' ),
    'white_mesh',ReadDiskItem( 'Hemisphere White Mesh' , shfjGlobals.aimsMeshFormats),
    'pole',WriteDiskItem( 'insula pole texture','Aims Texture formats' )
)

def initialization( self ):
    self.linkParameters( 'white_mesh', 'input_texture')
    self.linkParameters( 'pole', 'input_texture')
 

    
def execution( self, context ):
    re = aims.Reader()
    ws = aims.Writer()
    context.write('Reading textures and mesh')
    input_tex = re.read(self.input_texture.fullPath())
    mesh = re.read(self.white_mesh.fullPath())

    atex = np.zeros(input_tex[0].arraydata().shape)
    atex[input_tex[0].arraydata() > 0] = 180
    context.write('Topological correction...')  
    insular_tex_clean, insula_tex_boundary = textureTls.textureTopologicalCorrection(mesh, atex, 180)

    tex_out = aims.TimeTexture_S16()
    tex_out[0].assign(insular_tex_clean)
    ws.write(tex_out, self.pole.fullPath())
    context.write('... Done')


            
      
