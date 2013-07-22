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
  from brainvisa.cortical_surface.parameterization.mapping import hop
#  from brainvisa.cortical_surface.surface_tools import surface_tools as surfTls
except:
  pass
    
name = 'Harmonic Orthogonal Parameterization (HOP)'

userLevel = 2

signature = Signature(
                      
    'side', Choice('right', 'left'),    
    'square_mesh',ReadDiskItem( 'Rectangular flat mesh', shfjGlobals.aimsMeshFormats),
    'boundary_texture',ReadDiskItem( 'Rectangular boundary texture', 'Texture'),
    'corresp_indices_texture',ReadDiskItem( 'Rectangular flat indices texture', 'Texture'),
    'white_sulcalines',ReadDiskItem( 'hemisphere Sulcal Lines texture', 'Texture' ),
    'sulcus_labels',ReadDiskItem( 'Graph Label Translation', 'Text File'),
    'cstr_square_mesh',WriteDiskItem( 'Rectangular flat mesh', shfjGlobals.aimsMeshFormats)
)

def initialization( self ):
    self.linkParameters( 'boundary_texture','square_mesh')
    self.linkParameters( 'corresp_indices_texture','square_mesh')
    self.linkParameters( 'white_sulcalines', 'square_mesh')
    self.linkParameters( 'sulcus_labels', 'square_mesh')
    self.linkParameters( 'cstr_square_mesh','square_mesh')

    
def execution( self, context ):

#     lon, lat = hipHop(mesh, insula_pole[0].arraydata(), cing_pole[0].arraydata(), texture_sulci[0].arraydata(), self.side)

    re = aims.Reader()
    ws = aims.Writer()
    context.write('Reading textures and mesh')
    mesh = re.read(self.square_mesh.fullPath())
    tex_sulci = re.read(self.white_sulcalines.fullPath())
    tex_corresp_indices = re.read(self.corresp_indices_texture.fullPath())
    boundary_tex = re.read(self.boundary_texture.fullPath())
    context.write(self.sulcus_labels.fullPath())
    sulc_labels = []
    with open(self.sulcus_labels.fullPath(),'r') as inf:
        for line in inf:
            sulc_labels.append(line.split())    
    sulc_labels_dict = dict((int(value), key) for (key, value) in sulc_labels)
    
    '''
    boundaries (see mapping.path2Boundary for details:"
    boundary[0] == insula_boundary
    boundary[1] == neocortex_poles_path always from insula to cingular pole
    boundary[2] == cingular_boundary
    boundary[3] == new vertices always from insula to cingular pole
    '''
    boundary = []
    for t in  xrange( boundary_tex.size() ):
        boundary.append(np.where(boundary_tex[t].arraydata()>0)[0])

    square_mesh_indices = np.where( tex_corresp_indices[0].arraydata() )[0]

    tex_square_sulci_tmp = tex_sulci[0].arraydata()[square_mesh_indices]
    tex_square_sulci = np.hstack((tex_square_sulci_tmp, tex_square_sulci_tmp[boundary[1]]))

    labels = np.unique(tex_square_sulci)
    context.write('found the following sulci in the texture :')
    context.write([sulc_labels_dict[lab] for lab in labels])
    context.write('associated to the following labels :')
    context.write(labels)
    context.write('HOP')
 
     
    (cstr_mesh) = hop(mesh, boundary, tex_square_sulci, sulc_labels_dict, self.side)
    context.write('Writing meshes and textures')
    
    ws.write( cstr_mesh, self.cstr_square_mesh.fullPath() )
   
    context.write('Done')
      