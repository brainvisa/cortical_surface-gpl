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
from brainvisa.processes import *
import shfjGlobals 
from brainvisa.cortical_surface.parameterization.mapping import sphericalMeshFromCoords
from soma import aims

#from brainvisa import anatomist

name = 'Spherical mesh from HIP-HOP parameterization coordinates'

userLevel = 3

# def validation():
#     anatomist.validation()
    
signature = Signature(
                      
    'side', Choice('left', 'right'),
    'white_mesh',ReadDiskItem( 'Hemisphere White Mesh' , shfjGlobals.aimsMeshFormats ),
    'latitude',ReadDiskItem( 'Latitude coordinate texture','Texture' ),
    'longitude',ReadDiskItem( 'Longitude coordinate texture','Texture' ),
    'sphere_ray', Float(),
    'spherical_mesh', WriteDiskItem( 'spherical mesh', 'Aims mesh formats' )
)

def initialization( self ):
    self.side='left'
    self.linkParameters( 'latitude','white_mesh')
    self.linkParameters( 'longitude','white_mesh')
    self.linkParameters( 'spherical_mesh','white_mesh')
    self.sphere_ray = 100 
    
def execution( self, context ):
       
    re = aims.Reader()
    ws = aims.Writer()
    context.write('Reading textures and mesh')
    tex_lon = re.read(self.longitude.fullPath())
    tex_lat = re.read(self.latitude.fullPath())
    spherical_verts = sphericalMeshFromCoords(tex_lat[0].arraydata(), tex_lon[0].arraydata(), self.sphere_ray)

    vv = aims.vector_POINT3DF()
    for x in spherical_verts:
        vv.append(x)
    mesh = re.read(self.white_mesh.fullPath())
    mesh.vertex().assign(vv)
    mesh.updateNormals()
    ws.write( mesh, self.spherical_mesh.fullPath() )

    context.write('Done')
            
      