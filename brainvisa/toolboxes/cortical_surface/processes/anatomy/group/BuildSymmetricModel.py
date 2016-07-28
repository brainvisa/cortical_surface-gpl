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
  from brainvisa.cortical_surface.surface_tools import readSulcusLabelTranslationFile as rSLT
except:
  pass
    
name = 'Build a Symmetric Model'

userLevel = 0

signature = Signature(

    'model_left_in',ReadDiskItem( 'HipHop Model', 'Text File'),
    'model_right_in',ReadDiskItem( 'HipHop Model', 'Text File'),
    'model_symmetric',WriteDiskItem( 'HipHop Model', 'Text File'),
    'model_symmetric_mesh',WriteDiskItem( 'Mesh','Mesh mesh'),
#    'union_sulcal_lines_mesh',WriteDiskItem( 'Mesh', shfjGlobals.aimsMeshFormats),
#    'union_sulcal_lines_texture',WriteDiskItem( 'Texture', 'Texture')
)

def initialization( self ):

    self.setOptional('model_symmetric_mesh')
    
def execution( self, context ):

    #model_left = md.Model()
    model_left= md.Model().read(self.model_left_in.fullPath())
    #model_right = md.Model()
    model_right = md.Model().read(self.model_right_in.fullPath())

    ## check that the models from the two sides can be combined
    ok = 1
    if model_left.modelVersion == model_right.modelVersion:
        context.write('OK : the models from left and right side have the same modelVersion')
        if float(model_left.modelVersion) < 2:
            context.warning('models version is old, the axis - sulci correspondence are set to default')
    else:
        context.error('the models from left and right side are not compatible !')
        ok = 0
    if cmp(model_left.longitudeAxisID, model_right.longitudeAxisID) == 0 and cmp(model_left.latitudeAxisID, model_right.latitudeAxisID) == 0:
        context.write('OK : the models from left and right side have the same longitudeAxisID and latitudeAxisID')
    else:
        context.error('the models from left and right side are not compatible !')
        ok = 0
    if cmp(model_left.longitudeAxisSulci, model_right.longitudeAxisSulci) == 0 and cmp(model_left.latitudeAxisSulci, model_right.latitudeAxisSulci) == 0:
        context.write('OK : the models from left and right side have the same longitudeAxisSulci and latitudeAxisSulci')
    else:
        context.error('the models from left and right side are not compatible !')
        ok = 0


####################################################################
# the symmetric model is the average of left and right models
# it works directly because the rectangular representation of right hemisphere is automatically flipped !
####################################################################
    if ok:
        def mean_list(list1,list2):
            output = list()
            for (x, y) in zip(list1,list2):
                if x is None or y is None:
                    output.append(None)
                else:
                    output.append((x + y) / 2.0)
            return output


        # intialize the output model from _left one as _left and _right are compatible
        model_sym = md.Model()
        model_sym.longitudeAxisID = model_left.longitudeAxisID
        model_sym.latitudeAxisID = model_left.latitudeAxisID


        left = (model_left.left + model_right.left) / 2
        right = (model_left.right + model_right.right) / 2
        top = (model_left.top + model_right.top) / 2
        bottom = (model_left.bottom + model_right.bottom) / 2
        model_sym.setBoundary(left, right, top, bottom)
        model_sym.longitudeAxisCoord =mean_list(model_left.longitudeAxisCoord, model_right.longitudeAxisCoord)
        model_sym.latitudeAxisCoord =  mean_list(model_left.latitudeAxisCoord, model_right.latitudeAxisCoord)


        model_sym.saveToFile(self.model_symmetric.fullPath())
        context.write('------------------- output model -------------------')
        for line in model_sym.printArgs().splitlines():
            context.write(line)

        if self.model_symmetric_mesh is not None:
            context.write('saving model mesh')
            context.write(model_sym.toMesh())
            context.write(self.model_symmetric_mesh.fullPath())
            aims.write(model_sym.toMesh(),self.model_symmetric_mesh.fullPath())

