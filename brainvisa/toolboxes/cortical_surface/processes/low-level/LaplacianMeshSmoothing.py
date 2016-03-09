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
    from brainvisa.cortical_surface.surface_tools import PDE_tools as pdeTls
  except:
    raise ValidationError( 'brainvisa.cortical_surface.surface_tools.PDE_tools can not be imported.' )

from brainvisa.processes import *
try:
    from brainvisa.cortical_surface.surface_tools import PDE_tools as pdeTls
except:
    pass

name = 'Laplacian Mesh Smoothing'
userLevel = 0
signature = Signature(
    'mesh',ReadDiskItem( 'Mesh', 'aims mesh formats' ),
    'dt', Float(),
    'nb_iterations', Integer(),
    'smoothed_mesh', WriteDiskItem('Mesh', 'aims mesh formats')
)

def initialization( self ):
    self.dt = 0.7
    self.nb_iterations = 60
    
def execution( self, context ):
       
    re = aims.Reader()
    ws = aims.Writer()
    
    mesh = re.read(self.mesh.fullPath())
    mesh_smooth = pdeTls.laplacianMeshSmoothing(mesh ,self.nb_iterations, self.dt)
    ws.write(mesh_smooth, self.smoothed_mesh.fullPath())
    context.write('Done')
            
      
