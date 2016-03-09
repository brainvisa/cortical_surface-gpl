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
#from brainvisa import anatomist
import sigraph
import sys
from soma import aims, aimsalgo
import numpy as np

name = 'Project Texture Onto Atlas From Sphere'

userLevel = 0

# def validation():
#     anatomist.validation()
    
signature = Signature(
    'subject_texture', ReadDiskItem('Texture','Aims Texture formats'),
    'subject_spherical_mesh', ReadDiskItem( 'Hemisphere spherical mesh', 'Aims mesh formats' ),
    'atlas_spherical_mesh', ReadDiskItem( 'Hemisphere spherical mesh', 'Aims mesh formats' ),
    'texture_on_atlas', WriteDiskItem('Texture','Aims Texture formats')

)

def initialization( self ):
    self.linkParameters( 'subject_spherical_mesh','subject_texture' )
#    self.findValue( 'spherical_template', {'filename_variable' : 'ico100_7'} )
    
def execution( self, context ):
       
    re = aims.Reader()
    ws = aims.Writer()
    
    subject_spherical_mesh = re.read(self.subject_spherical_mesh.fullPath())
    atlas_spherical_mesh = re.read(self.atlas_spherical_mesh.fullPath())
    mi = aims.MeshInterpoler(subject_spherical_mesh, atlas_spherical_mesh)
    mi.project() # calcule les correspondances et coord barycentriques
    subject_texture = re.read(self.subject_texture.fullPath())
    a = subject_texture[0].arraydata()
    isint = a.dtype.type in ( np.int, np.int8, np.int16, np.int32, np.int64 )
    if isint:
        context.write('texture type is interger, using nearest neighbour interpolation')
        texture_on_atlas = mi.resampleTexture(subject_texture, mi.NearestNeighbour)
    else:
        context.write('texture type is float, using linear interpolation')
        texture_on_atlas = mi.resampleTexture(subject_texture)
    ws.write( texture_on_atlas, self.texture_on_atlas.fullPath() )

    context.write('Done')
            
      
