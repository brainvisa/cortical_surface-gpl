# -*- coding: utf-8 -*-
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
        from brainvisa.cortical_surface.surface_tools import PDE_tools as pdeTls
    except:
        raise ValidationError( 'brainvisa.cortical_surface.parameterization.surface_tools module can not be imported.' )
  

from brainvisa.processes import *
try:
    from brainvisa.cortical_surface.surface_tools import PDE_tools as pdeTls
except:
    pass

name = 'Laplacian Eigen Vectors'
userLevel = 0

# Argument declaration
signature = Signature(
    'input_mesh',ReadDiskItem( 'Hemisphere White Mesh' , 'Aims mesh formats' ),
    'Eigen_vectors_texture',WriteDiskItem( 'Eigen vectors texture',
    'Aims texture formats' ),
    'number_of_vectors', Integer(),
)


# Default values
def initialization( self ):
    self.linkParameters( 'Eigen_vectors_texture', 'input_mesh' )
    self.number_of_vectors =1


def execution( self, context ):
    re = aims.Reader()
    ws = aims.Writer()
    
    mesh = re.read(self.input_mesh.fullPath())
    context.write('computing the '+str(self.number_of_vectors)+' first Laplacian Eigen vector(s) of the mesh')

    # tex_vectors = aims.TimeTexture_FLOAT()
    # for i in range(20):
    #     vectors = pdeTls.meshLaplacianEigenVectors(mesh, self.number_of_vectors)
    #     tex_vectors[i].assign(vectors[:,0])
    vectors = pdeTls.meshLaplacianEigenVectors(mesh, self.number_of_vectors)
    tex_vectors = aims.TimeTexture_FLOAT()
    for ind in range(self.number_of_vectors):
        tex_vectors[ind].assign(vectors[:,ind])

    ws.write(tex_vectors, self.Eigen_vectors_texture.fullPath())
    context.write('... Done')
