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


from __future__ import absolute_import
from six.moves import range
def validation():
  try:
    import brainvisa.cortical_surface.surface_tools.basic_tools
  except:
    raise ValidationError( 'brainvisa.cortical_surface.surface_tools.basic_tools module can not be imported.' )
  
from brainvisa.processes import *
from soma import aims
import numpy as np

try:
  from brainvisa.cortical_surface.surface_tools import basic_tools as basicTls
  from brainvisa.cortical_surface.surface_tools import texture_tools as texTls
except:
  pass

#from brainvisa import anatomist

name = 'Poles Textures Sanity Check'

userLevel = 0

# def validation():
#     anatomist.validation()
    
signature = Signature(

    'white_mesh',ReadDiskItem( 'Hemisphere White Mesh', 'aims mesh formats' ),
    'side', Choice('left', 'right'),
    'cingular_pole_texture_in',ReadDiskItem( 'Cingular pole texture', 'aims Texture formats'),
    'insula_pole_texture_in',ReadDiskItem( 'Insula pole texture', 'aims Texture formats'),
    'dilation_size', Integer(),
    'cingular_pole_texture_out',WriteDiskItem( 'Cingular pole texture', 'aims Texture formats'),
    'insula_pole_texture_out',WriteDiskItem( 'Insula pole texture', 'aims Texture formats')
)

def initialization( self ):
    def linkSide( proc, dummy ):
        if proc.white_mesh is not None:
            return proc.white_mesh.get( 'side' )
    self.linkParameters( 'side', 'white_mesh', linkSide )
    self.linkParameters( 'cingular_pole_texture_in', 'white_mesh')
    self.linkParameters( 'insula_pole_texture_in', 'white_mesh')
    self.linkParameters( 'cingular_pole_texture_out', 'white_mesh')
    self.linkParameters( 'insula_pole_texture_out', 'white_mesh')
    self.dilation_size = 2
    
def execution( self, context ):
  
    re = aims.Reader()
    ws = aims.Writer()
    context.write('Reading textures and mesh')
    cingular_tex_clean_in = re.read(self.cingular_pole_texture_in.fullPath())
    insula_tex_clean_in = re.read(self.insula_pole_texture_in.fullPath())
    mesh = re.read(self.white_mesh.fullPath())
    acingular_tex_clean_in = cingular_tex_clean_in[0].arraydata().copy()
    ainsula_tex_clean_in = insula_tex_clean_in[0].arraydata().copy()

    cingular_tex_value = acingular_tex_clean_in.max()
    insula_tex_value = ainsula_tex_clean_in.max()
    neigh = aims.SurfaceManip.surfaceNeighbours(mesh)

    save_cingular = False
    cingular_boundary = basicTls.textureBoundary(mesh, acingular_tex_clean_in, cingular_tex_value, neigh)
    if len(cingular_boundary) > 1:
        save_cingular = True
        context.write('Topological correction needed for cingular pole...')
        acingular_tex_clean_in, cingular_tex_boundary = texTls.textureTopologicalCorrection(mesh, acingular_tex_clean_in, cingular_tex_value, 0, neigh)
    save_insula = False
    insula_boundary = basicTls.textureBoundary(mesh, ainsula_tex_clean_in, insula_tex_value, neigh)
    if len(insula_boundary) > 1:
        save_insula = True
        context.write('Topological correction needed for insular pole...')
        ainsula_tex_clean_in, insula_tex_boundary = texTls.textureTopologicalCorrection(mesh, ainsula_tex_clean_in, insula_tex_value, 0, neigh)

    ## dilation of the two poles will ensure that the gap between them is least of 2*dilation_size vertices
    for it in range(self.dilation_size):
        simple_boundary = basicTls.textureSimpleBoundary(mesh, acingular_tex_clean_in, cingular_tex_value, neigh)
        acingular_tex_clean_in = texTls.dilateTexture(acingular_tex_clean_in, neigh, simple_boundary)
        simple_boundary = basicTls.textureSimpleBoundary(mesh, ainsula_tex_clean_in, insula_tex_value, neigh)
        ainsula_tex_clean_in = texTls.dilateTexture(ainsula_tex_clean_in, neigh, simple_boundary)

    cingular_inds_in = np.where(acingular_tex_clean_in > 0)[0]
    insula_inds_in = np.where(ainsula_tex_clean_in >0 )[0]
    # test if there is an intersection between the two poles
    inter = list(set(cingular_inds_in).intersection(insula_inds_in))
    if inter:
        context.write('Cingular and Insula poles are connected by '+str(len(inter))+' vertices')
        context.write('dilation of the intersection')
        inter_tex = np.zeros(mesh.vertex().size())
        inter_tex[inter] = 1
        for it in range(self.dilation_size):
            simple_boundary = basicTls.textureSimpleBoundary(mesh, inter_tex, 1, neigh)
            inter_tex = texTls.dilateTexture(inter_tex, neigh, simple_boundary)
        context.write('deleting intersection from insular and cingular poles')
        inds_inter = np.where(inter_tex >0 )[0]
        acingular_tex_clean_out = cingular_tex_clean_in[0].arraydata()
        ainsula_tex_clean_out = insula_tex_clean_in[0].arraydata()
        acingular_tex_clean_out[inds_inter] = 0
        ainsula_tex_clean_out[inds_inter] = 0
        ## topological correction of the new poles textures
        context.write('topological correction of the new cingular pole texture')
        acingular_tex_clean_out, boundary = texTls.textureTopologicalCorrection(mesh, acingular_tex_clean_out, cingular_tex_value)
        context.write('topological correction of the new insula pole texture')
        ainsula_tex_clean_out, boundary = texTls.textureTopologicalCorrection(mesh, ainsula_tex_clean_out,  insula_tex_value)
        tex_out = aims.TimeTexture_S16()
        tex_out[0].assign(acingular_tex_clean_out.astype(np.int16))
        ws.write(tex_out, self.cingular_pole_texture_out.fullPath())
        tex_out = aims.TimeTexture_S16()
        tex_out[0].assign(ainsula_tex_clean_out.astype(np.int16))
        ws.write(tex_out, self.insula_pole_texture_out.fullPath())

    else:
        context.write('Cingular and Insular poles are separated by more than '+str(2*self.dilation_size)+' vertices, OK')
        if save_cingular:
            acingular_tex_clean_in[acingular_tex_clean_in > 0] = 1
            tex_out = aims.TimeTexture_S16()
            tex_out[0].assign(acingular_tex_clean_in.astype(np.int16))
            ws.write(tex_out, self.cingular_pole_texture_out.fullPath())
        if save_insula:
            ainsula_tex_clean_in[ainsula_tex_clean_in > 0] = 180
            tex_out = aims.TimeTexture_S16()
            tex_out[0].assign(ainsula_tex_clean_in.astype(np.int16))
            ws.write(tex_out, self.insula_pole_texture_out.fullPath())
    context.write('Done')

