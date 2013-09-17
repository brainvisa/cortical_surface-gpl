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
    import brainvisa.cortical_surface.parameterization.mapping
  except:
    raise ValidationError( 'brainvisa.cortical_surface.parameterization.mapping module can not be imported.' )
  
from brainvisa.processes import *
import shfjGlobals  
from soma import aims
import numpy as np

try:
  from brainvisa.cortical_surface.parameterization import mapping as map
  from brainvisa.cortical_surface.surface_tools import surface_tools as surfTls
  from brainvisa.cortical_surface.parameterization import model as md
except:
  pass
    
name = 'Harmonic Orthogonal Parameterization (HOP)'

userLevel = 2

signature = Signature(
                      
    'side', Choice('right', 'left'),    
    'rectangular_mesh',ReadDiskItem( 'Rectangular flat mesh', shfjGlobals.aimsMeshFormats),
    'boundary_texture',ReadDiskItem( 'Rectangular boundary texture', 'Texture'),
    'corresp_indices_texture',ReadDiskItem( 'Rectangular flat indices texture', 'Texture'),
    'white_sulcalines',ReadDiskItem( 'hemisphere Sulcal Lines Rectangular Flat texture', 'Texture' ),
    'cstrBalance', Float(),
#    'white_sulcalines',ReadDiskItem( 'hemisphere Sulcal Lines texture', 'Texture' ),
    'sulcus_labels',ReadDiskItem( 'Graph Label Translation', 'Text File'),
    'model_file',ReadDiskItem( 'Graph Label Translation', 'Text File'),
    'unfold_reversed_triangles', Choice('yes','no'),
    'nb_it_local_smoothing_for_unfolding', Integer(),
    'cstr_rectangular_mesh',WriteDiskItem( 'Rectangular flat cstr mesh', shfjGlobals.aimsMeshFormats)
)

def initialization( self ):
    self.linkParameters( 'boundary_texture','rectangular_mesh')
    self.linkParameters( 'corresp_indices_texture','rectangular_mesh')
    self.linkParameters( 'white_sulcalines', 'rectangular_mesh')
    self.cstrBalance = 200
    self.linkParameters( 'sulcus_labels', 'rectangular_mesh')
    self.linkParameters( 'cstr_rectangular_mesh','rectangular_mesh')
    self.unfold_reversed_triangles = 'yes'
    self.nb_it_local_smoothing_for_unfolding = 100
    
def execution( self, context ):
    context.write('Reading model')
    model = md.Model().read(self.model_file.fullPath())
    for line in model.printArgs().splitlines():
        context.write(line)
    re = aims.Reader()
    ws = aims.Writer()
    context.write('Reading textures and mesh')
    mesh = re.read(self.rectangular_mesh.fullPath())
    tex_square_sulci = re.read(self.white_sulcalines.fullPath())
#    tex_sulci = re.read(self.white_sulcalines.fullPath())
    tex_corresp_indices = re.read(self.corresp_indices_texture.fullPath())
    boundary_tex = re.read(self.boundary_texture.fullPath())
    context.write('Reading sulcus-label correspondences file')
    
    sulc_labels_dict = surfTls.readSulcusLabelTranslationFile(self.sulcus_labels.fullPath())
    
    '''
    boundaries (see mapping.path2Boundary for details:"
    boundary[0] == insula_boundary
    boundary[1] == neocortex_poles_path always from insula to cingular pole
    boundary[2] == cingular_boundary
    boundary[3] == new vertices always from insula to cingular pole
    '''
    boundary = []
    for t in  range( boundary_tex.size() ):
        boundary.append(np.where(boundary_tex[t].arraydata()>0)[0])

#     square_mesh_indices = np.where( tex_corresp_indices[0].arraydata() )[0]
# 
#     tex_square_sulci_tmp = tex_sulci[0].arraydata()[square_mesh_indices]
#     tex_square_sulci = np.hstack((tex_square_sulci_tmp, tex_square_sulci_tmp[boundary[1]]))
    square_sulci = tex_square_sulci[0].arraydata()
    labels = np.unique(square_sulci)
    context.write('found the following sulci in the texture :')
    context.write([sulc_labels_dict[lab] for lab in labels])
    context.write('associated to the following labels :')
    context.write(labels)
    context.write('HOP')
    (cstr_mesh) = map.hop(self.cstrBalance, mesh, boundary, square_sulci, sulc_labels_dict, self.side, model)
    (nb_inward, inward) = map.invertedPolygon(cstr_mesh)
    vert = np.array(cstr_mesh.vertex())
    context.write('------------------number of vertices on folded triangles : '+str(nb_inward)+' => '+str(100.0 * nb_inward / vert.shape[0])+' %')

    if self.unfold_reversed_triangles == 'yes':
        context.write('------------------unfolding reversed triangles')
        (cstr_mesh, nb_inward_evol, inward_evol) = map.solveInvertedPolygon(cstr_mesh, boundary, self.nb_it_local_smoothing_for_unfolding)
        context.write('------------------number of vertices on folded triangles : '+str(nb_inward_evol))
#         inward_tex = 'tmp.tex'
#         context.write('------------------writing inward tex in : '+inward_tex)
#         tmp_tex = np.zeros(len(mesh.vertex()))
#         tmp_tex[inward_evol[-1]] = 1
#         tex_unfold = aims.TimeTexture_S16()
#         tex_unfold[0].assign(tmp_tex)
#         ws.write(tex_unfold, inward_tex)

    context.write('Writing meshes and textures')
    
    ws.write( cstr_mesh, self.cstr_rectangular_mesh.fullPath() )
   
    context.write('Done')
      