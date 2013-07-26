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
  from brainvisa.cortical_surface.parameterization import mapping as map
  from brainvisa.cortical_surface.surface_tools import surface_tools as surfTls
except:
  pass

#from brainvisa import anatomist

name = 'CoordinatesFromHiphopMapping'

userLevel = 2

# def validation():
#     anatomist.validation()
    
signature = Signature(
    'cstr_rectangular_mesh',ReadDiskItem( 'Rectangular flat cstr mesh', shfjGlobals.aimsMeshFormats),        
    'boundary_texture',ReadDiskItem( 'Rectangular boundary texture', 'Texture'),
    'corresp_indices_texture',ReadDiskItem( 'Rectangular flat indices texture', 'Texture'),
    'white_mesh_parts',ReadDiskItem( 'White Mesh Parts', shfjGlobals.aimsMeshFormats),
    'latitude_insula_boundary', Float(),
    'latitude_cingular_pole_boundary', Float(),
    'latitude',WriteDiskItem( 'Latitude coordinate texture','Texture'),
    'longitude',WriteDiskItem( 'Longitude coordinate texture','Texture')

)

def initialization( self ):
    self.linkParameters('boundary_texture', 'cstr_rectangular_mesh')
    self.linkParameters('corresp_indices_texture', 'cstr_rectangular_mesh')
    self.linkParameters('white_mesh_parts', 'cstr_rectangular_mesh')
    self.linkParameters('latitude', 'cstr_rectangular_mesh')
    self.linkParameters('longitude', 'cstr_rectangular_mesh')
    self.latitude_insula_boundary = 30
    self.latitude_cingular_pole_boundary = 30

    
def execution( self, context ):
    re = aims.Reader()
    ws = aims.Writer()
    context.write('Reading textures and mesh')
    neoCortex_square_cstr = re.read(self.cstr_rectangular_mesh.fullPath())
    mesh_parts = re.read(self.white_mesh_parts.fullPath())
    tex_corresp_indices = re.read(self.corresp_indices_texture.fullPath())
    boundary_tex = re.read(self.boundary_texture.fullPath())
    boundary = []
    for t in  range( boundary_tex.size() ):
        boundary.append(np.where(boundary_tex[t].arraydata()>0)[0])
    
    '''
    mesh_parts[0] = neoCortex
    mesh_parts[1] = insula
    mesh_parts[2] = cingular pole
    '''
    insula_mesh = aims.AimsTimeSurface_3()
    insula_mesh.vertex().assign( mesh_parts.vertex(1) )
    insula_mesh.normal().assign( mesh_parts.normal(1) )
    insula_mesh.polygon().assign( mesh_parts.polygon(1) )
    cingular_mesh = aims.AimsTimeSurface_3()
    cingular_mesh.vertex().assign( mesh_parts.vertex(2) )
    cingular_mesh.normal().assign( mesh_parts.normal(2) )
    cingular_mesh.polygon().assign( mesh_parts.polygon(2) )

    '''
    tex_corresp_indices contains the indices of the vertices in white_mesh for:
        neoCortex_square in time 0
        insula_indices in time 1
        cingular_indices in time 2
    '''
    neocortex_indices = np.where( tex_corresp_indices[0].arraydata() )[0]
    insula_indices= np.where( tex_corresp_indices[1].arraydata() )[0]
    cingular_indices = np.where( tex_corresp_indices[2].arraydata() )[0]
    
    nb_vert_full_mesh = len(neocortex_indices) + len(insula_indices) + len(cingular_indices) 
    lon, lat = map.computeCoordinates(nb_vert_full_mesh, neocortex_indices, neoCortex_square_cstr, boundary, self.latitude_insula_boundary, self.latitude_cingular_pole_boundary)

    '''
    if side is left
    invert the bpundary
    '''
    context.write('mapping the insula to a disk')
    insula_lon = map.texture2ROI(lon, insula_indices)
    insula_boundary = surfTls.meshBoundary(insula_mesh)
    context.write('insula_lon[insula_boundary]')
    context.write(insula_lon[insula_boundary])
    (insula_lon, insula_lat, insula_disk) = map.mesh2Disk(insula_mesh, insula_boundary, insula_lon)
    context.write('insula_lon = [',np.min(insula_lon),', ',np.max(insula_lon),']')
    context.write('insula_lat = [',np.min(insula_lat),', ',np.max(insula_lat),']')

    lon[insula_indices] = insula_lon
    lat[insula_indices] = insula_lat * self.latitude_insula_boundary
    context.write('mapping the cingular pole to a disk')    
    cingular_boundary = surfTls.meshBoundary(cingular_mesh)
    cingular_lon = map.texture2ROI(lon, cingular_indices)
    (cingular_lon, cingular_lat, cingular_disk) = map.mesh2Disk(cingular_mesh, cingular_boundary, cingular_lon)
    context.write('cingular_lon = [', np.min(cingular_lon),', ',np.max(cingular_lon),']')
    context.write('cingular_lat = [', np.min(cingular_lat),', ',np.max(cingular_lat),']')

    lon[cingular_indices] = cingular_lon
    lat[cingular_indices] = 180 - cingular_lat * self.latitude_cingular_pole_boundary
    context.write('Writing textures')
    tex_lon = aims.TimeTexture_FLOAT()
    tex_lon[0].assign(lon)
    ws.write(tex_lon, self.longitude.fullPath())
    tex_lat = aims.TimeTexture_FLOAT()
    tex_lat[0].assign(lat)
    ws.write(tex_lat, self.latitude.fullPath())
    ws.write(insula_disk, 'insula_disk.mesh')
    ws.write(cingular_disk, 'cingular_disk.mesh')
    