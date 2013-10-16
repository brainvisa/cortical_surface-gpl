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

name = 'Remesh From Sphere'

userLevel = 0

# def validation():
#     anatomist.validation()
    
signature = Signature(
                      
    'white_mesh',ReadDiskItem( 'Hemisphere White Mesh' , shfjGlobals.aimsMeshFormats ),    
    'spherical_mesh', ReadDiskItem( 'Spherical mesh', 'Aims mesh formats' ),
    'spherical_template', ReadDiskItem( 'Spherical mesh', 'Aims mesh formats' ),
    'remeshed_mesh', WriteDiskItem( 'Remeshed mesh', 'Aims mesh formats' )
)

def initialization( self ):
    self.linkParameters( 'white_mesh','spherical_mesh')
    self.linkParameters( 'spherical_mesh', 'white_mesh')
    self.linkParameters( 'remeshed_mesh','spherical_mesh')
    self.linkParameters( 'remeshed_mesh','white_mesh')
#    self.findValue( 'spherical_template', {'filename_variable' : 'ico100_7'} )
    self.sphere_ray = 100 
    
def execution( self, context ):
       
    re = aims.Reader()
    ws = aims.Writer()
    
    spherical_mesh = re.read(self.spherical_mesh.fullPath())
    spherical_template_mesh = re.read(self.spherical_template.fullPath())
    mi = aims.MeshInterpoler(spherical_mesh, spherical_template_mesh)
    mi.project() # calcule les correspondances et coord barycentriques
    white_mesh = re.read(self.white_mesh.fullPath())
    outmesh = mi.resampleMesh(white_mesh)
    ws.write( outmesh, self.remeshed_mesh.fullPath() )

    context.write('Done')
            
      
