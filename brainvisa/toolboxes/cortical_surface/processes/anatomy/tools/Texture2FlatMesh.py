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
    import brainvisa.cortical_surface.parameterization.mapping
  except:
    raise ValidationError( 'brainvisa.cortical_surface.parameterization.mapping module can not be imported.' )
  
from brainvisa.processes import *
import shfjGlobals  
from soma import aims
import numpy as np

try:
  from brainvisa.cortical_surface.parameterization.mapping import texture2ROI
except:
  pass

#from brainvisa import anatomist

name = 'Texture to Flat Mesh'

userLevel = 2

# def validation():
#     anatomist.validation()
    
signature = Signature(
    'input_texture',ReadDiskItem( 'hemisphere Sulcal Lines texture', 'Texture' ),                     
#    'input_texture',ReadDiskItem('Texture', 'Texture'),
    'corresp_indices_texture',ReadDiskItem( 'Rectangular flat indices texture', 'Texture'),
    'boundary_texture',ReadDiskItem( 'Rectangular boundary texture', 'Texture'),
    'output_texture',WriteDiskItem( 'hemisphere Sulcal Lines Rectangular Flat texture', 'Texture' )
    #'output_texture',WriteDiskItem( 'Rectangular flat texture', 'Texture')
)

def initialization( self ):
    self.linkParameters( 'boundary_texture','input_texture')
    self.linkParameters( 'corresp_indices_texture','input_texture')
    self.linkParameters( 'boundary_texture','input_texture')
    self.linkParameters( 'output_texture','input_texture')

    
def execution( self, context ):
  
    re = aims.Reader()
    ws = aims.Writer()
    context.write('Reading textures and mesh')
    input_tex = re.read(self.input_texture.fullPath())
    tex_corresp_indices = re.read(self.corresp_indices_texture.fullPath())
    boundary_tex = re.read(self.boundary_texture.fullPath())
    context.write('texture to flat mesh')
    
    '''
    boundaries (see mapping.path2Boundary for details:"
    boundary[0] == insula_boundary
    boundary[1] == neocortex_poles_path always from insula to cingular pole
    boundary[2] == cingular_boundary
    boundary[3] == new vertices always from insula to cingular pole
    '''
    new_verts_indices = np.where(boundary_tex[3].arraydata()>0)[0]
    inv_new_verts_indices = np.where(boundary_tex[1].arraydata()>0)[0]
    boundary0_indices = np.where(boundary_tex[0].arraydata()>0)[0]
    boundary1_indices = np.where(boundary_tex[2].arraydata()>0)[0]
    '''
    tex_corresp_indices contains the indices of the vertices in white_mesh for:
        neoCortex_square in time 0
        insula_indices in time 1
        cingular_indices in time 2
    '''
    rectangular_mesh_indices = np.where( tex_corresp_indices[0].arraydata() )[0]
    print rectangular_mesh_indices
    print len(rectangular_mesh_indices)
    print len(new_verts_indices)
    nb_vert_square = len(rectangular_mesh_indices) + len(new_verts_indices)
    print nb_vert_square
    print new_verts_indices
    print inv_new_verts_indices
    output_tex = aims.TimeTexture_S16()
    for t in  range( input_tex.size() ):
        output_tex_tmp = input_tex[t].arraydata()[rectangular_mesh_indices]
        tmp_tex = np.zeros(nb_vert_square)
        tmp_tex[range( len(rectangular_mesh_indices) )] = output_tex_tmp
        tmp_tex[range( len(rectangular_mesh_indices) )] = output_tex_tmp
 #       tmp_tex = np.hstack((output_tex_tmp, output_tex_tmp[inv_new_verts_indices]))
        output_tex[t].assign(tmp_tex)
    ws.write(output_tex, self.output_texture.fullPath())
    
    context.write('Done')
            
      
