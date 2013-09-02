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
  from brainvisa.cortical_surface.surface_tools import surface_tools as surfTls
except:
  pass
    
name = 'Build a Model from a set of subjects'

userLevel = 2

signature = Signature(
                      
    'side', Choice('right', 'left'),    
    'rectangular_mesh',ListOf( ReadDiskItem( 'Rectangular flat mesh', shfjGlobals.aimsMeshFormats) ),
    'boundary_texture',ListOf( ReadDiskItem( 'Rectangular boundary texture', 'Texture') ),
#    'corresp_indices_texture',ReadDiskItem( 'Rectangular flat indices texture', 'Texture'),
    'flat_white_sulcalines',ListOf( ReadDiskItem( 'hemisphere Sulcal Lines Rectangular Flat texture', 'Texture' ) ),
#    'white_sulcalines',ReadDiskItem( 'hemisphere Sulcal Lines texture', 'Texture' ),
    'sulcus_labels',ListOf( ReadDiskItem( 'Graph Label Translation', 'Text File') ),
    'model_file',WriteDiskItem( 'Graph Label Translation', 'Text File')
)

def initialization( self ):
    self.linkParameters( 'boundary_texture','rectangular_mesh')
#    self.linkParameters( 'corresp_indices_texture','rectangular_mesh')
    self.linkParameters( 'flat_white_sulcalines', 'rectangular_mesh')
    self.linkParameters( 'sulcus_labels', 'rectangular_mesh')
#    self.linkParameters( 'cstr_rectangular_mesh','rectangular_mesh')

    
def execution( self, context ):

    re = aims.Reader()
    nb_mesh = len(self.rectangular_mesh)

    model = md.Model()
    for ind_mesh,r_mesh in enumerate(self.rectangular_mesh):
        context.write('working on mesh nb: ',ind_mesh+1)
        mesh = re.read(r_mesh.fullPath())
        tex_square_sulci = re.read(self.flat_white_sulcalines[ind_mesh].fullPath())
        sulci_dict = surfTls.readSulcusLabelTranslationFile(self.sulcus_labels[ind_mesh].fullPath())
        full_sulci = slSet.SulcalLinesSet()
        full_sulci.extractFromTexture(tex_square_sulci[0].arraydata(), mesh, sulci_dict)
        SC_ind = full_sulci.names.index(('S.C._'+self.side))   
        SC_label = full_sulci.labels[SC_ind]
        print 'SC_label: ', SC_label
#        full_sulci.sulcalLines[SC_ind].printArgs()
        translation = -full_sulci.sulcalLines[SC_ind].barycenter[0]
        vert = np.array(mesh.vertex())
        vert[:, 0] = vert[:, 0] + translation # * np.ones(vert.shape[0])
        full_sulci.updateVertices(vert)
        full_sulci.sulcalLine2SulcalConstraint(model)
        if ind_mesh==0:
            group_full_sulci = full_sulci
        else:
            group_full_sulci.cat(full_sulci)

    model.setBoundary(np.min(vert[:, 0]), np.max(vert[:, 0]), np.min(vert[:, 1]), np.max(vert[:, 1]))
    model.setAxisCoord(group_full_sulci)
    context.write('model built from '+str(nb_mesh)+' subjects')
    model.printArgs()
    model.saveToFile(self.model_file.fullPath())
    model2 = md.Model()
    test = model2.read(self.model_file.fullPath())
    context.write(test.printArgs())
