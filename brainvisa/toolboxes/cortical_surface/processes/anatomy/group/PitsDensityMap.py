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
from soma import aims
import numpy as np

    
name = 'Compute Pits Density From a Set of Subjects'

userLevel = 0

signature = Signature(
                      
    'projected_smoothed_pits',ListOf( ReadDiskItem( 'Texture', 'aims Texture Formats') ),
    'pits_density',WriteDiskItem( 'Texture', 'aims Texture Formats'),
)

#def initialization( self ):

    
def execution( self, context ):

    tex1 = aims.read(self.projected_smoothed_pits[0].fullPath())
    atex1 = tex1[0].arraydata()
    good_shape = atex1.shape
    a_tex_out = np.zeros(good_shape)
    nb_subj = 0
    for proj_pits in self.projected_smoothed_pits:
        tex_pits = aims.read(proj_pits.fullPath())
        atex_pits = np.array(tex_pits[0])
        if atex_pits.shape == good_shape:
            a_tex_out += np.nan_to_num(atex_pits)
            nb_subj += 1
        else:
            pb_subj = proj_pits.get( 'subject' )
            context.error('problem, subject '+pb_subj+' not processed')
    a_tex_out = a_tex_out / nb_subj
    tex_out = aims.TimeTexture_FLOAT()
    tex_out[0].assign(a_tex_out)
    aims.write(tex_out, self.pits_density.fullPath())
    context.write('done')