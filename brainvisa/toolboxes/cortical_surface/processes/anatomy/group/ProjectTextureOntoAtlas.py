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
from brainvisa.processes import *
import sys
from soma import aims, aimsalgo
import numpy as np

name = 'Project Texture Onto Atlas From Sphere'

userLevel = 0


signature = Signature(
    'subject_texture', ReadDiskItem('Texture','Aims Texture formats'),
    'subject_spherical_mesh', ReadDiskItem( 'Hemisphere spherical mesh', 'Aims mesh formats' ),
    'atlas_spherical_mesh', ReadDiskItem( 'Hemisphere spherical mesh', 'Aims mesh formats' ),
    'texture_on_atlas', WriteDiskItem('Texture','Aims Texture formats')

)

def initialization( self ):
    self.linkParameters( 'subject_spherical_mesh','subject_texture' )
#    self.findValue( 'spherical_template', {'filename_variable' : 'ico100_7'} )

def nearest_neighbor(vert_template,vert_pits):
    vertex_number1=vert_template.shape[0]
    nn=[]
    for v in vert_pits:
        nn_tmp = np.argmin(np.sum(np.square(np.tile(v,(vertex_number1,1))-vert_template),1))
        nn.append(nn_tmp)
    return nn

def execution( self, context ):

    re = aims.Reader()
    ws = aims.Writer()

    subject_spherical_mesh = re.read(self.subject_spherical_mesh.fullPath())
    subject_texture = re.read(self.subject_texture.fullPath())
    a = subject_texture[0].arraydata()
    # test whether the mesh and texture have the same number of element to prevent the crash from mi.resampleTexture
    vert = np.array(subject_spherical_mesh.vertex())
    if vert.shape[0] != a.shape[0]:
        context.error('subject_spherical_mesh and subject_texture do not have the same number of elements')

    else:
        atlas_spherical_mesh = re.read(self.atlas_spherical_mesh.fullPath())
        isint = a.dtype.type in ( np.int, np.int8, np.int16, np.int32, np.int64 )
        if isint:
            context.write('texture type is interger, using nearest neighbour interpolation')
            #texture_on_atlas = mi.resampleTexture(subject_texture, mi.NearestNeighbour)
            pits=np.where(a==1)[0]
            vert_template=np.array(atlas_spherical_mesh.vertex())
            vert_pits = vert[pits,:]
            nn=nearest_neighbor(vert_template,vert_pits)
            context.write(len(set(nn)), 'of', len(pits), 'points interpolated')

            a_tex_out = np.zeros(vert_template.shape[0])
            a_tex_out[nn]=1
            texture_on_atlas =aims.TimeTexture_S16()
            texture_on_atlas[0].assign(a_tex_out)

        else:
            context.write('texture type is float, using linear interpolation')
            mi = aims.MeshInterpoler(subject_spherical_mesh, atlas_spherical_mesh)
            mi.project() # calcule les correspondances et coord barycentriques
            texture_on_atlas = mi.resampleTexture(subject_texture)
        ws.write( texture_on_atlas, self.texture_on_atlas.fullPath() )
        context.write('Done')
