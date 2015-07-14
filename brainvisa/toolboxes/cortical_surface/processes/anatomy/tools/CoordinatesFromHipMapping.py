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
from soma import aims
import numpy as np

try:
  from brainvisa.cortical_surface.parameterization import mapping as map
  from brainvisa.cortical_surface.surface_tools import basic_tools as basicTls
  from brainvisa.cortical_surface.parameterization import model as md
except:
  pass

#from brainvisa import anatomist

name = 'Coordinates From Hip Mapping'

userLevel = 0

# def validation():
#     anatomist.validation()
    
signature = Signature(
    'rectangular_mesh',ReadDiskItem( 'Rectangular flat mesh', 'aims mesh formats' ),
    'boundary_texture',ReadDiskItem( 'Rectangular boundary texture', 'aims Texture formats'),
    'corresp_indices_texture',ReadDiskItem( 'Rectangular flat indices texture', 'aims Texture formats'),
    'white_mesh_parts',ReadDiskItem( 'White Mesh Parts', 'aims mesh formats' ),
    'model_file',ReadDiskItem( 'HipHop Model', 'Text File'),
#    'latitude_insula_boundary', Float(),
#    'latitude_cingular_pole_boundary', Float(),
    'latitude',WriteDiskItem( 'Latitude coordinate HIP texture', 'aims Texture formats' ),
    'longitude',WriteDiskItem( 'Longitude coordinate HIP texture', 'aims Texture formats' )

)

def initialization( self ):
    self.linkParameters('boundary_texture', 'rectangular_mesh')
    self.linkParameters('corresp_indices_texture', 'rectangular_mesh')
    self.linkParameters('white_mesh_parts', 'rectangular_mesh')
    self.linkParameters('latitude', 'rectangular_mesh')
    self.linkParameters('longitude', 'rectangular_mesh')
#    self.latitude_insula_boundary = 30
#    self.latitude_cingular_pole_boundary = 30

    
def execution( self, context ):
    re = aims.Reader()
    ws = aims.Writer()
    context.write('Reading textures, meshes and model')
    model = md.Model().read(self.model_file.fullPath())
    for line in model.printArgs().splitlines():
        context.write(line)

    neoCortex_square_cstr = re.read(self.rectangular_mesh.fullPath())
    mesh_parts = re.read(self.white_mesh_parts.fullPath())
    tex_corresp_indices = re.read(self.corresp_indices_texture.fullPath())
    boundary_tex = re.read(self.boundary_texture.fullPath())
    '''
    boundaries (see mapping.path2Boundary for details:"
    boundary[0] == insula_boundary
    boundary[1] == neocortex_poles_path always from insula to cingular pole
    boundary[2] == cingular_boundary
    boundary[3] == new vertices always from insula to cingular pole
    '''
    boundary = []
    for t in  range( boundary_tex.size() ):
        boundary.append(np.where(boundary_tex[t].arraydata()>0)[0])
    
    '''
    mesh_parts[0] = neoCortex
    mesh_parts[1] = insula
    mesh_parts[2] = cingular pole
    '''
    insula_mesh = aims.AimsTimeSurface_3_VOID()
    insula_mesh.vertex().assign( mesh_parts.vertex(1) )
    insula_mesh.normal().assign( mesh_parts.normal(1) )
    insula_mesh.polygon().assign( mesh_parts.polygon(1) )
    cingular_mesh = aims.AimsTimeSurface_3_VOID()
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
#    print np.concatenate((neocortex_indices, insula_indices),0)
    nb_vert_full_mesh = len(np.unique(np.concatenate((np.concatenate((neocortex_indices, insula_indices), 0), cingular_indices), 0)))
#    context.write('nb_vert_full_mesh'+str(nb_vert_full_mesh))
#     nb_vert_full_mesh = len(neocortex_indices) + len(insula_indices) + len(cingular_indices) 
#     context.write('nb_vert_full_mesh'+str(nb_vert_full_mesh))
    lon, lat = map.computeCoordinates(nb_vert_full_mesh, neocortex_indices, neoCortex_square_cstr, boundary, model.insularPoleBoundaryCoord, model.cingularPoleBoundaryCoord)

#     context.write('lon.shape'+str(lon.shape))
#     context.write('lat.shape'+str(lon.shape))
#     context.write('neoCortex_square_cstr.vertex().size '+str(neoCortex_square_cstr.vertex().size()))
    
    '''
    if side is left
    invert the boundary
    '''

    context.write('mapping the insula to a disk')
    insula_lon = map.texture2ROI(lon, insula_indices)
    insula_boundary = basicTls.meshBoundary(insula_mesh)[0]
#     context.write(insula_boundary)
#     context.write(len(insula_boundary))
# 
# 
#     context.write('insula_lon[insula_boundary]')
#     context.write(len(insula_lon[insula_boundary]))
#     context.write('lon[boundary[0]]')
#     context.write(len(lon[boundary[0]]))
#     context.write('lon[boundary[0]] = [',np.min(lon[boundary[0]]),', ',np.max(lon[boundary[0]]),']')
#     context.write('insula_lon[insula_boundary] = [',np.min(insula_lon[insula_boundary]),', ',np.max(insula_lon[insula_boundary]),']')
#     #(a,b) = np.sort(insula_lon[insula_boundary])
#     context.write(insula_boundary)

#     i_tmp = insula_lon[insula_boundary]
#     for i,b in enumerate(insula_boundary):
#         context.write(i_tmp[i])
    insula_bound_rad = np.pi * (insula_lon[insula_boundary] - 180) / 180 
    circle = np.array([np.cos(insula_bound_rad), np.sin(insula_bound_rad)])
#     context.write(circle)
    bound_mesh = aims.AimsTimeSurface_2_VOID()
    vv = aims.vector_POINT3DF()
    ee = aims.vector_AimsVector_U32_2()
    for i in range(len(circle)):
        vv.append([circle[i, 0], circle[i, 1], 0])
    for i in range(len(insula_boundary) - 1):
        ee.append(np.array([i, i + 1], np.uint32))
    bound_mesh.vertex().assign(vv)
    bound_mesh.polygon().assign(ee)
    bound_mesh.updateNormals()
    #ws.write(bound_mesh, '/home/toz/ammon_Lwhite_insula_bound_cirlce.mesh' )

    
    
    #ws.write(basicTls.meshBoundaryMesh(insula_mesh, [insula_boundary]), '/home/toz/ammon_Lwhite_insula_bound.mesh' )
    tex_insula_boundary_lon = aims.TimeTexture_FLOAT()
    tex_insula_boundary_lon[0].assign(insula_lon[insula_boundary])
    #ws.write(tex_insula_boundary_lon,  '/home/toz/ammon_Lwhite_insula_boundary_lon.tex')

    tex_insula_lon = aims.TimeTexture_FLOAT()
    tex_insula_lon[0].assign(insula_lon)
    #ws.write(tex_insula_lon,  '/home/toz/ammon_Lwhite_insula_lon.tex')

    (insula_lon, insula_lat, insula_disk) = map.mesh2Disk(insula_mesh, insula_boundary, insula_lon)
#     context.write('insula_lon = [',np.min(insula_lon),', ',np.max(insula_lon),']')
#     context.write('insula_lat = [',np.min(insula_lat),', ',np.max(insula_lat),']')

#     import matplotlib.pyplot as plt
#     plt.figure(1)
#     plt.hold(True)
#     vert = np.array(insula_disk.vertex())
#     bord = insula_boundary
# #    vert = array(mesh_out.vertex())
#     #a = plt.gca()
# #     for p in poly:
# #         plt.plot(vert[[p[0],p[1]],0],vert[[p[0],p[1]],1],'b')
# #         plt.plot(vert[[p[1],p[2]],0],vert[[p[1],p[2]],1],'b')
# #         plt.plot(vert[[p[2],p[0]],0],vert[[p[2],p[0]],1],'b')
# #        a.plot(vert[[p[0],p[1]],0],vert[[p[0],p[1]],1])
#     plt.plot(vert[:,0],vert[:,1],'ro')
#     plt.plot(vert[bord,0],vert[bord,1],'r')
#     plt.show()

    lon[insula_indices] = insula_lon
    lat[insula_indices] = insula_lat * model.insularPoleBoundaryCoord
    context.write('mapping the cingular pole to a disk')    
    cingular_boundary = basicTls.meshBoundary(cingular_mesh)[0]
    cingular_lon = map.texture2ROI(lon, cingular_indices)
    (cingular_lon, cingular_lat, cingular_disk) = map.mesh2Disk(cingular_mesh, cingular_boundary, cingular_lon)
#     context.write('cingular_lon = [', np.min(cingular_lon),', ',np.max(cingular_lon),']')
#     context.write('cingular_lat = [', np.min(cingular_lat),', ',np.max(cingular_lat),']')

    lon[cingular_indices] = cingular_lon
    lat[cingular_indices] = 180 - cingular_lat * model.cingularPoleBoundaryCoord
    context.write('Writing textures')
    tex_lon = aims.TimeTexture_FLOAT()
    tex_lon[0].assign(lon)
    ws.write(tex_lon, self.longitude.fullPath())
    tex_lat = aims.TimeTexture_FLOAT()
    tex_lat[0].assign(lat)
    ws.write(tex_lat, self.latitude.fullPath())
#    ws.write(insula_disk, 'insula_disk.mesh')
#    ws.write(cingular_disk, 'cingular_disk.mesh')
    