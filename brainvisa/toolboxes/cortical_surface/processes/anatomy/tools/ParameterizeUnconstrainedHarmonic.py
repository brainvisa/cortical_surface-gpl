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
  from brainvisa.cortical_surface.parameterization.mapping import hip#hipHop
#  from brainvisa.cortical_surface.surface_tools import surface_tools as surfTls
except:
  pass

#from brainvisa import anatomist

name = 'Harmonic Intrinsic Parameterization (HIP)'

userLevel = 2

# def validation():
#     anatomist.validation()
    
signature = Signature(
                      
    'side', Choice('right', 'left'),    
    'white_mesh',ReadDiskItem( 'Hemisphere White Mesh', shfjGlobals.aimsMeshFormats),    
    'cingular_pole_texture',ReadDiskItem( 'Hippocampus pole texture', 'Texture'),
    'insular_pole_texture',ReadDiskItem( 'Insula pole texture', 'Texture'),
    'rectangular_mesh',WriteDiskItem( 'Rectangular flat mesh', shfjGlobals.aimsMeshFormats),
    'boundary_texture',WriteDiskItem( 'Rectangular boundary texture', 'Texture'),
    'corresp_indices_texture',WriteDiskItem( 'Rectangular flat indices texture', 'Texture')
)

def initialization( self ):
    self.linkParameters( 'cingular_pole_texture', 'white_mesh')
    self.linkParameters( 'insular_pole_texture', 'white_mesh')
    self.linkParameters( 'rectangular_mesh','white_mesh')
    self.linkParameters( 'boundary_texture','white_mesh')
    self.linkParameters( 'corresp_indices_texture','white_mesh')

    
def execution( self, context ):
  
    re = aims.Reader()
    ws = aims.Writer()
    context.write('Reading textures and mesh')
    cing_pole = re.read(self.cingular_pole_texture.fullPath())
    insula_pole = re.read(self.insular_pole_texture.fullPath())
    mesh = re.read(self.white_mesh.fullPath())
    context.write('HIP')
 
     
    (neoCortex_square, neoCortex_open_boundary, neocortex_indices, insula_indices, cingular_indices, insula_mesh, cingular_mesh, neoCortex_mesh) = hip(mesh, insula_pole[0].arraydata(), cing_pole[0].arraydata())
    context.write('Writing meshes and textures')
    
    '''boundaries (see mapping.path2Boundary for details:"
    boundary[0] == insula_boundary
    boundary[1] == neocortex_poles_path always from insula to cingular pole
    boundary[2] == cingular_boundary
    boundary[3] == new vertices always from insula to cingular pole
    '''
    tex_boundary = aims.TimeTexture_S16()
    for ind,bound in enumerate(neoCortex_open_boundary):
        tmp_tex = np.zeros(len(neoCortex_square.vertex()))
        tmp_tex[bound] = 1
        tex_boundary[ind].assign(tmp_tex)
    ws.write(tex_boundary, self.boundary_texture.fullPath())
    ws.write( neoCortex_square, self.rectangular_mesh.fullPath() )
    '''
    tex_corresp_indices contains the indices of the vertices in white_mesh for:
        neoCortex_square in time 0
        insula_indices in time 1
        cingular_indices in time 2
    '''
    tex_corresp_indices = aims.TimeTexture_S16()
    tmp_tex = np.zeros(len(mesh.vertex()))
    tmp_tex[neocortex_indices] = 1
    tex_corresp_indices[0].assign(tmp_tex)
    tmp_tex = np.zeros(len(mesh.vertex()))
    tmp_tex[insula_indices] = 1
    tex_corresp_indices[1].assign(tmp_tex)
    tmp_tex = np.zeros(len(mesh.vertex()))
    tmp_tex[cingular_indices] = 1
    tex_corresp_indices[2].assign(tmp_tex)
    ws.write(tex_corresp_indices, self.corresp_indices_texture.fullPath())

#     re = aims.Reader()
#     ws = aims.Writer()
#     context.write('Reading textures and mesh')
#     cing_pole = re.read(self.cingular_pole_texture.fullPath())
#     insula_pole = re.read(self.insular_pole_texture.fullPath())
#     texture_sulci = re.read(self.white_sulcalines.fullPath())
#     mesh = re.read(self.white_mesh.fullPath())
#     context.write('HIP-HOP')
#     context.write(np.unique(texture_sulci[0].arraydata()))
# #     from brainvisa.cortical_surface.parameterization import sulcalLinesSet as slSet
# #     full_sulci = slSet.SulcalLinesSet()
# #     full_sulci.extractFromTexture(texture_sulci[0].arraydata(), mesh)
# #     full_sulci.printArgs()
# 
#     
#     lon, lat = hipHop(mesh, insula_pole[0].arraydata(), cing_pole[0].arraydata(), texture_sulci[0].arraydata(), self.side)
#     context.write('Writing textures')
#     tex_lon = aims.TimeTexture_FLOAT()
#     tex_lon[0].assign(lon)
#     ws.write(tex_lon, self.longitude.fullPath())
#     tex_lat = aims.TimeTexture_FLOAT()
#     tex_lat[0].assign(lat)
#     ws.write(tex_lat, self.latitude.fullPath())

#    spherical_verts = sphericalMeshFromCoords(lat, lon, 50):
    
    context.write('Done')
            
      
