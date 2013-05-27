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
from soma import aims
import numpy as np
 
#from brainvisa import anatomist

name = 'Harmonic Intrinsic Parameterization (HIP)'

userLevel = 3

#def validation():
#    anatomist.validation()
    
signature = Signature(
                      
    'side', Choice('right', 'left'),    
    'white_mesh',ReadDiskItem( 'Hemisphere White Mesh' , shfjGlobals.aimsMeshFormats),    
    'cingular_pole_texture',ReadDiskItem( 'Hippocampus pole texture','Texture'),
    'insular_pole_texture',ReadDiskItem( 'Insula pole texture','Texture'),
    'white_sulcalines',ReadDiskItem( 'hemisphere Sulcal Lines texture', 'Texture' ),
    'sulcus_labels',ReadDiskItem( 'Graph Label Translation', 'Text File'),
    'latitude',WriteDiskItem( 'Latitude coordinate texture','Texture'),
    'longitude',WriteDiskItem( 'Longitude coordinate texture','Texture')
)

def initialization( self ):
    self.linkParameters( 'cingular_pole_texture','white_mesh')
    self.linkParameters( 'insular_pole_texture','white_mesh')
    self.linkParameters( 'white_sulcalines','white_mesh')
    self.linkParameters( 'sulcus_labels','white_mesh')
    self.linkParameters( 'latitude','white_mesh')
    self.linkParameters( 'longitude','white_mesh')

    
def execution( self, context ):
#    sys.path.append('/home/toz/workspace/MyTestProject/cortical_surface')
    for p in sys.path:
        print p
    print np.__version__ 
    from brainvisa.cortical_surface.parameterization.mapping import hipHop
    print 'mapping imported'
   
  
    re = aims.Reader()
    ws = aims.Writer()
    context.write('Reading textures and mesh')
    cing_pole = re.read(self.cingular_pole_texture.fullPath())
    insula_pole = re.read(self.insular_pole_texture.fullPath())
    texture_sulci = re.read(self.white_sulcalines.fullPath())
    mesh = re.read(self.white_mesh.fullPath())
    context.write('HIP-HOP')
    lon, lat = hipHop(mesh, insula_pole[0].arraydata(), cing_pole[0].arraydata(), texture_sulci[0].arraydata())
    context.write('Writing textures')
    tex_lon = aims.TimeTexture_FLOAT()
    tex_lon[0].assign(lon)
    ws.write(tex_lon, self.longitude.fullPath())
    tex_lat = aims.TimeTexture_FLOAT()
    tex_lat[0].assign(lat)
    ws.write(tex_lat, self.latitude.fullPath())

#    spherical_verts = sphericalMeshFromCoords(lat, lon, 50):
    
    context.write('Done')
            
      
