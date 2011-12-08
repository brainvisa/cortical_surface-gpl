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
from neuroProcesses import *
import shfjGlobals   
from brainvisa import anatomist
import sigraph

name = 'Harmonic mapping under orthogonal constraints left hemisphere'

userLevel = 2

def validation():
    anatomist.validation()
    
signature = Signature(
                      
    'Side', Choice('left'),
    'Lwhite_mesh',ReadDiskItem( 'Hemisphere White Mesh' , shfjGlobals.aimsMeshFormats),
    
)

def initialization( self ):

    self.side = 'left'
    
def execution( self, context ):
        
    context.write('Done')
            
      