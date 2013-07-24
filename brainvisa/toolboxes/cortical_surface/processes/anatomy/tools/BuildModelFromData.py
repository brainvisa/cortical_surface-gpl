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
    import brainvisa.cortical_surface.parameterization.model
  except:
    raise ValidationError( 'brainvisa.cortical_surface.parameterization.model module can not be imported.' )
  
from brainvisa.processes import *
import shfjGlobals  
from soma import aims
import numpy as np

try:
  from brainvisa.cortical_surface.parameterization.mapping import hop
  from brainvisa.cortical_surface.parameterization import sulcalLinesSet as slSet
  from brainvisa.cortical_surface.parameterization import model as md
#  from brainvisa.cortical_surface.surface_tools import surface_tools as surfTls
except:
  pass
    
name = 'Build a Model from a set of subjects'

userLevel = 2

signature = Signature(
                      
    'side', Choice('right', 'left'),    
    'rectangular_mesh',ListOf( ReadDiskItem( 'Rectangular flat mesh', shfjGlobals.aimsMeshFormats) ),
    'boundary_texture',ListOf( ReadDiskItem( 'Rectangular boundary texture', 'Texture') ),
#    'corresp_indices_texture',ReadDiskItem( 'Rectangular flat indices texture', 'Texture'),
    'white_sulcalines',ListOf( ReadDiskItem( 'hemisphere Sulcal Lines Rectangular Flat texture', 'Texture' ) ),
#    'white_sulcalines',ReadDiskItem( 'hemisphere Sulcal Lines texture', 'Texture' ),
    'sulcus_labels',ListOf( ReadDiskItem( 'Graph Label Translation', 'Text File') ),
    'model_file',WriteDiskItem( 'Graph Label Translation', 'Text File')
)

def initialization( self ):
    self.linkParameters( 'boundary_texture','rectangular_mesh')
#    self.linkParameters( 'corresp_indices_texture','rectangular_mesh')
    self.linkParameters( 'white_sulcalines', 'rectangular_mesh')
    self.linkParameters( 'sulcus_labels', 'rectangular_mesh')
#    self.linkParameters( 'cstr_rectangular_mesh','rectangular_mesh')

    
def execution( self, context ):

    re = aims.Reader()


    model = md.Model()
    for ind_mesh,r_mesh in enumerate(self.rectangular_mesh):
        mesh = re.read(r_mesh)
        tex_square_sulci = re.read(self.white_sulcalines.fullPath())
        sulci_dict = surfTls.readSulcusLabelTranslationFile(self.sulcus_labels.fullPath())
        full_sulci = slSet.SulcalLinesSet()
        full_sulci.extractFromTexture(tex_square_sulci[0].arraydata(), mesh, sulci_dict)
#     lon, lat = hipHop(mesh, insula_pole[0].arraydata(), cing_pole[0].arraydata(), texture_sulci[0].arraydata(), self.side)
