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
from six.moves import range
def validation():
  try:
    import brainvisa.cortical_surface.parameterization.mapping
  except:
    raise ValidationError( 'brainvisa.cortical_surface.parameterization.mapping module can not be imported.' )
  
from brainvisa.processes import *
from brainvisa.tools import aimsGlobals
from soma import aims
import numpy as np

try:
  from brainvisa.cortical_surface.parameterization import mapping as map
  from brainvisa.cortical_surface.surface_tools import readSulcusLabelTranslationFile as rSLT
  from brainvisa.cortical_surface.parameterization import model as md
  from brainvisa.cortical_surface.parameterization import sulcalLinesSet as slSet
except:
  pass
    
name = 'Harmonic Orthogonal Parameterization (HOP)'

userLevel = 0

signature = Signature(
                      
    'rectangular_mesh',ReadDiskItem( 'Rectangular flat mesh', 'aims mesh formats' ),
    'side', Choice('left', 'right'),
    'boundary_texture',ReadDiskItem( 'Rectangular boundary texture', 'aims Texture formats' ),
    'corresp_indices_texture',ReadDiskItem( 'Rectangular flat indices texture', 'aims Texture formats' ),
    'white_mesh_parts',ReadDiskItem( 'White Mesh Parts', 'aims mesh formats' ),
    'white_sulcalines',ReadDiskItem( 'hemisphere Sulcal Lines Rectangular Flat texture', 'aims Texture formats' ),
    'cstrBalance', Float(),
#    'white_sulcalines',ReadDiskItem( 'hemisphere Sulcal Lines texture', 'Texture' ),
    'sulcus_labels',ReadDiskItem( 'Graph Label Translation', 'Text File'),
    'model_file',ReadDiskItem( 'HipHop Model', 'Text File'),
    'unfold_reversed_triangles', Choice('yes','no'),
    'nb_it_local_smoothing_for_unfolding', Integer(),
    'cstr_rectangular_mesh',WriteDiskItem( 'Rectangular flat cstr mesh',
      'aims mesh formats' )
)

def initialization( self ):
    def linkSide( proc, dummy ):
        if proc.rectangular_mesh is not None:
            return proc.rectangular_mesh.get( 'side' )
    self.linkParameters( 'side', 'rectangular_mesh', linkSide )
    self.linkParameters( 'boundary_texture','rectangular_mesh')
    self.linkParameters( 'corresp_indices_texture','rectangular_mesh')
    self.linkParameters('white_mesh_parts', 'rectangular_mesh')
    self.linkParameters( 'white_sulcalines', 'rectangular_mesh')
    self.cstrBalance = 200
    self.linkParameters( 'sulcus_labels', 'rectangular_mesh')
    self.linkParameters( 'cstr_rectangular_mesh','rectangular_mesh')
    self.linkParameters( 'model_file','rectangular_mesh')
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
    
    sulc_labels_dict = rSLT.readSulcusLabelTranslationFile(self.sulcus_labels.fullPath())
    
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
    labels = labels[labels!=0]
    context.write('found the following sulci in the texture :')
    #context.write( 'labels:', labels )
    #context.write( 'sulc_labels_dict:', sulc_labels_dict )
    #context.write( 'missing:', [ lab not in sulc_labels_dict for lab in labels ] )
    context.write([sulc_labels_dict[lab] for lab in labels])
    context.write('associated to the following labels :')
    context.write(labels)
    context.write('HOP')
    mesh_parts = re.read(self.white_mesh_parts.fullPath())
    '''
    mesh_parts[0] = neoCortex_open
    mesh_parts[1] = insula
    mesh_parts[2] = cingular pole
    '''
    open_mesh = aims.AimsTimeSurface_3_VOID()
    open_mesh.vertex().assign( mesh_parts.vertex(0) )
    open_mesh.normal().assign( mesh_parts.normal(0) )
    open_mesh.polygon().assign( mesh_parts.polygon(0) )

    full_sulci = slSet.SulcalLinesSet()
    full_sulci.extractFromTexture(square_sulci, mesh, sulc_labels_dict)

    (cstr_mesh) = map.hop(self.cstrBalance, mesh, boundary,open_mesh, full_sulci, self.side, model)
    (nb_inward, inward) = map.invertedPolygon(cstr_mesh)
    #vert = np.array(cstr_mesh.vertex())
    context.write('------------------number of vertices on folded triangles : '+str(nb_inward)+' => '+str(100.0 * nb_inward / len(cstr_mesh.polygon()))+' %')

    if self.unfold_reversed_triangles == 'yes' and nb_inward>0:
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
      
