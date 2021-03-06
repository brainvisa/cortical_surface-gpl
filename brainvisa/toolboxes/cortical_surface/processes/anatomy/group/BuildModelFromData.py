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


from __future__ import absolute_import
def validation():
  try:
    import brainvisa.cortical_surface.parameterization.model
  except:
    raise ValidationError( 'brainvisa.cortical_surface.parameterization.model module can not be imported.' )
  
from brainvisa.processes import *
from soma import aims
import numpy as np

try:
  from brainvisa.cortical_surface.parameterization.mapping import hop
  from brainvisa.cortical_surface.parameterization import sulcalLinesSet as slSet
  from brainvisa.cortical_surface.parameterization import model as md
  from brainvisa.cortical_surface.surface_tools import readSulcusLabelTranslationFile as rSLT
except:
  pass
    
name = 'Build a Model from a set of subjects'

userLevel = 0

signature = Signature(
                      
    'rectangular_mesh',ListOf( ReadDiskItem( 'Rectangular flat mesh', 'aims Mesh Formats') ),
    'side', Choice('left', 'right'),
    'boundary_texture',ListOf( ReadDiskItem( 'Rectangular boundary texture', 'Texture') ),
#    'corresp_indices_texture',ReadDiskItem( 'Rectangular flat indices texture', 'Texture'),
    'flat_white_sulcalines',ListOf( ReadDiskItem( 'hemisphere Sulcal Lines Rectangular Flat texture', 'Texture' ) ),
#    'white_sulcalines',ReadDiskItem( 'hemisphere Sulcal Lines texture', 'Texture' ),
    'sulcus_labels',ListOf( ReadDiskItem( 'Graph Label Translation', 'Text File') ),
    'model_file_in',ReadDiskItem( 'HipHop Model', 'Text File'),
    'model_file',WriteDiskItem( 'HipHop Model', 'Text File'),
    'model_file_mesh',WriteDiskItem( 'mesh', 'Mesh mesh'),
    'union_sulcal_lines_mesh',WriteDiskItem( 'mesh', 'Mesh mesh'),
    'union_sulcal_lines_texture',WriteDiskItem( 'Texture', 'Texture') 
)

def initialization( self ):
    def linkSide( proc, dummy ):
        if proc.rectangular_mesh is not None \
                and len( proc.rectangular_mesh ) != 0:
            return proc.rectangular_mesh[0].get( 'side' )
    self.linkParameters( 'side', 'rectangular_mesh', linkSide )
    self.linkParameters( 'boundary_texture','rectangular_mesh')
#    self.linkParameters( 'corresp_indices_texture','rectangular_mesh')
    self.linkParameters( 'flat_white_sulcalines', 'rectangular_mesh')
    self.linkParameters( 'sulcus_labels', 'rectangular_mesh')
    self.setOptional('model_file_mesh', 'union_sulcal_lines_mesh', 'union_sulcal_lines_texture', 'model_file_in')
    
def execution( self, context ):

    re = aims.Reader()
    nb_mesh = len(self.rectangular_mesh)


    model = md.Model()
    if self.model_file_in is not None:
        model = model.read(self.model_file_in.fullPath())

    left = 0
    right = 0
    top = 0
    bottom = 0
    for ind_mesh,r_mesh in enumerate(self.rectangular_mesh):
        context.write('working on mesh nb: ',ind_mesh+1)
        mesh = re.read(r_mesh.fullPath())
        tex_square_sulci = re.read(self.flat_white_sulcalines[ind_mesh].fullPath())
        sulci_dict = rSLT.readSulcusLabelTranslationFile(self.sulcus_labels[ind_mesh].fullPath(), invert=True)
        # extract only the sulci that are used in the model
        model_sulci_dict = {}
        for s in model.longitudeAxisSulci.keys():
            sulci_dict_key = s+'_'+self.side
            if sulci_dict_key in sulci_dict:
                model_sulci_dict[sulci_dict[sulci_dict_key]] = sulci_dict_key
        for s in model.latitudeAxisSulci.keys():
            sulci_dict_key = s+'_'+self.side
            if sulci_dict_key in sulci_dict:
                model_sulci_dict[sulci_dict[sulci_dict_key]] = sulci_dict_key
        full_sulci = slSet.SulcalLinesSet()
        full_sulci.extractFromTexture(tex_square_sulci[0].arraydata(), mesh, model_sulci_dict)
        context.write('Translating the barycenter of S.C. to 0')
        SC_ind = full_sulci.names.index(('S.C._'+self.side))   
        SC_label = full_sulci.labels[SC_ind]
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
        left +=np.min(vert[:, 0])
        right += np.max(vert[:, 0])
        top += np.max(vert[:, 1])
        bottom += np.min(vert[:, 1])

####################################################################
# the coordinates of the boundaries of the model corresponds to the barycenter 
# of the boundaries of the rectangular meshes given in self.rectangular_mesh
####################################################################    
    left = left / nb_mesh
    right = right / nb_mesh
    top = top / nb_mesh
    bottom = bottom / nb_mesh
    context.write('right-left '+str(right-left))
    context.write('top-bottom '+str(top-bottom))
    model.setBoundary(left, right, top, bottom)
    model.setAxisCoord(group_full_sulci)
    context.write('model built from '+str(nb_mesh)+' subjects')
    model.saveToFile(self.model_file.fullPath())
    context.write('------------------- output model -------------------')
    for line in model.printArgs().splitlines():
        context.write(line)

    if self.model_file_mesh is not None:
        context.write('saving model mesh')
        model.saveToMesh(self.model_file_mesh.fullPath())
    if self.union_sulcal_lines_mesh is not None:
        context.write('saving the union of sulcal lines mesh')
        ws = aims.Writer()
        ws.write(group_full_sulci.toMesh(), self.union_sulcal_lines_mesh.fullPath())
    if self.union_sulcal_lines_texture is not None:
        context.write('saving the union of sulcal lines texture')
        ws = aims.Writer()
        ws.write(group_full_sulci.toTex(), self.union_sulcal_lines_texture.fullPath())