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

try :
  from brainvisa.cortical_surface.parameterization.mapping import sphericalMeshFromCoords
except :
  pass


#from brainvisa import anatomist

name = 'Spherical mesh from HIP-HOP parameterization coordinates'

userLevel = 0
    
signature = Signature(
                      
    'white_mesh',ReadDiskItem( 'Hemisphere White Mesh' , shfjGlobals.aimsMeshFormats ),
    'side', Choice('left', 'right'),
    'latitude',ReadDiskItem( 'Latitude coordinate texture','Texture' ),
    'longitude',ReadDiskItem( 'Longitude coordinate texture','Texture' ),
    'sphere_ray', Float(),
    'spherical_mesh', WriteDiskItem( 'spherical mesh', 'Aims mesh formats' )
)

def initialization( self ):
    def linkSide( proc, dummy ):
        if proc.white_mesh is not None:
            return proc.white_mesh.get( 'side' )
    self.linkParameters( 'side', 'white_mesh', linkSide )
    self.linkParameters( 'longitude','white_mesh')
    self.linkParameters( 'latitude','white_mesh')
    self.linkParameters( 'spherical_mesh','white_mesh')
    self.sphere_ray = 100 
    
def execution( self, context ):
       
    re = aims.Reader()
    ws = aims.Writer()
    context.write('Reading textures and mesh')
    tex_lon = re.read(self.longitude.fullPath())
    tex_lat = re.read(self.latitude.fullPath())
    spherical_verts = sphericalMeshFromCoords(tex_lat[0].arraydata(), tex_lon[0].arraydata(), self.sphere_ray)#, self.side)
    vv = aims.vector_POINT3DF()
    for x in spherical_verts:
        vv.append(x)
#     for i in range(spherical_verts.shape[0]):
#         vv.append([spherical_verts[i, 0], spherical_verts[i, 1], spherical_verts[i, 2]])
    mesh = re.read(self.white_mesh.fullPath())
    new_mesh = aims.AimsTimeSurface_3()
    new_mesh.vertex().assign(vv)
    new_mesh.polygon().assign(mesh.polygon())
    if self.side == 'right':
        aims.SurfaceManip.invertSurfacePolygons(new_mesh)
#         poly = np.array(mesh.polygon())
#         poly_tmp = poly.copy()
# #        context.write(poly_tmp[0,:])
#         poly_tmp[:,0] = poly[:,1]
#         poly_tmp[:,1] = poly[:,0]
#         pp = aims.vector_AimsVector_U32_3()
#         for i in poly_tmp:
#             pp.append(i)
# #        context.write(np.array(pp)[0,:])
#         new_mesh.polygon().assign(pp)
    new_mesh.updateNormals()
    ws.write( new_mesh, self.spherical_mesh.fullPath() )

    context.write('Done')
            
      