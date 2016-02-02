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
    import numpy as np
except:
    pass

name = 'Mesh Fiedler length'
userLevel = 0

# Argument declaration
signature = Signature(
    'input_meshes', ListOf( ReadDiskItem( 'Hemisphere White Mesh' , 'Aims mesh formats' ) ),
    'output_csv_file', WriteDiskItem( 'CSV file', 'CSV file' )
    #'fiedler_texture',WriteDiskItem( 'Texture', 'Aims texture formats' )
    #'distance_type'
)


# Default values
#def initialization( self ):
    #self.setOptional( 'fiedler_texture' )


def execution( self, context ):
    re = aims.Reader()
    nb_mesh = len( self.input_meshes )
    f = open( self.output_csv_file.fullPath(), 'w' )
    f.write( 'subject\tside\tfidler_length\n' )
    context.write('Computing the Fiedler geodesic length')
    group_lengths = []
    for ind_mesh, input_mesh in enumerate(self.input_meshes):
        context.progress( ind_mesh, nb_mesh, process=self )
        subject = input_mesh.get('subject')
        side = input_mesh.get('side')
        mesh = re.read(input_mesh.fullPath())
        (d, fiedler) = pdeTls.meshFiedlerLength(mesh, 'geodesic')
        f.write( subject + '\t' + side + '\t' + str(d) + '\n' )
        group_lengths.append(d)
    context.progress( nb_mesh, nb_mesh, process=self )
    avg = np.mean(group_lengths)
    f.write( 'AVERAGE\tboth\t' + str(avg) + '\n' )
    f.close()
    context.write('average length:')
    context.write(avg)
    # if self.fiedler_texture is not None:
    #     context.write('saving the texture')
    #     ws = aims.Writer()
    #     ws.write(fiedler, self.fiedler_texture.fullPath())
    context.write('Done')
