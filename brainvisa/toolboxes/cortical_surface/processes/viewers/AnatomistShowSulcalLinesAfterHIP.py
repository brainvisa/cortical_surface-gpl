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
from brainvisa.processes import *
from brainvisa import anatomist
from soma import aims
import numpy as np
from brainvisa.cortical_surface.parameterization import sulcalLinesSet as slSet
from brainvisa.cortical_surface.parameterization import model as md
from brainvisa.cortical_surface.surface_tools import readSulcusLabelTranslationFile as rSLT

name = 'Anatomist Show Sulcal Lines After HIP'
roles = ('viewer',)
userLevel = 0

def validation():
  anatomist.validation()

signature = Signature(
    'rectangular_mesh',ListOf( ReadDiskItem( 'Rectangular flat mesh', 'aims Mesh Formats') ),
    'side', Choice('left', 'right'),
#    'corresp_indices_texture',ReadDiskItem( 'Rectangular flat indices texture', 'Texture'),
    'flat_white_sulcalines',ListOf( ReadDiskItem( 'hemisphere Sulcal Lines Rectangular Flat texture', 'Texture' ) ),
#    'white_sulcalines',ReadDiskItem( 'hemisphere Sulcal Lines texture', 'Texture' ),
    'sulcus_labels',ListOf( ReadDiskItem( 'Graph Label Translation', 'Text File') ),
    'model_file',ReadDiskItem( 'HipHop Model', 'Text File'),
#    'model_file_mesh',WriteDiskItem( 'Mesh', 'aims Mesh Formats'),
#    'union_sulcal_lines_mesh',WriteDiskItem( 'Mesh', 'aims Mesh Formats'),
#    'union_sulcal_lines_texture',WriteDiskItem( 'Texture', 'Texture')
)



def initialization( self ):
    def linkSide( proc, dummy ):
        if proc.rectangular_mesh is not None \
                and len( proc.rectangular_mesh ) != 0:
            return proc.rectangular_mesh[0].get( 'side' )
    self.linkParameters( 'side', 'rectangular_mesh', linkSide )
    self.linkParameters( 'flat_white_sulcalines', 'rectangular_mesh')
    self.linkParameters( 'sulcus_labels', 'rectangular_mesh')
    self.setOptional('model_file')


def execution( self, context ):

    re = aims.Reader()
    nb_mesh = len(self.rectangular_mesh)


    model = md.Model()
    if self.model_file is not None:
        model = model.read(self.model_file.fullPath())


    for ind_mesh,r_mesh in enumerate(self.rectangular_mesh):
        context.write('loading mesh nb: ',ind_mesh+1)
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
        full_sulci.sulcalLine2SulcalConstraint(model)
        if ind_mesh==0:
            group_full_sulci = full_sulci
        else:
            group_full_sulci.cat(full_sulci)



    temp_mesh_file = context.temporary(  'GIFTI file' )
    aims.write(group_full_sulci.toMesh(), temp_mesh_file.fullPath())
    temp_tex_file = context.temporary(  'GIFTI file' )
    aims.write(group_full_sulci.toTex(), temp_tex_file.fullPath())


    a = anatomist.Anatomist()
    win = a.createWindow( 'Axial' )
    anamesh = a.loadObject(temp_mesh_file)
    anatex = a.loadObject(temp_tex_file)
    anapalette = a.getPalette('Graph-Label')
    anatex.setPalette( anapalette, minVal = 0, maxVal= 2.04)
    anatex.glSetTexRGBInterpolation(True)
    fusionTexSurf = a.fusionObjects( [anamesh, anatex], method='FusionTexSurfMethod' )
    win.addObjects(fusionTexSurf )
    if self.model_file is not None:
        temp_model_file = context.temporary(  'GIFTI file' )
        model.saveToMesh(temp_model_file.fullPath())
        anamod = a.loadObject(temp_model_file)
        win.addObjects(anamod )
        return [a,win,anamesh,anatex,anamod, fusionTexSurf]
    else:
        return [a,win,anamesh,anatex, fusionTexSurf]

