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

name = 'Sulcal Pits Extraction'
userLevel = 2

# Argument declaration
signature = Signature(
    'input_mesh',ReadDiskItem( 'Hemisphere White Mesh' , 'Aims mesh formats' ),
    'mask_texture',ReadDiskItem( 'Hippocampus pole texture','Aims Texture formats' ),
    'DPF_alpha', Float(),
    'thresh_ridge', Float(),
    'thresh_dist', Float(),
    'thresh_area', Float(),
    'DPF_texture',WriteDiskItem( 'DPF texture',  'Aims texture formats' ),
    'pits_texture',WriteDiskItem( 'pits texture',  'Aims texture formats' ),
    'noisypits_texture',WriteDiskItem( 'noisy pits texture',  'Aims texture formats' ),
    'ridges_texture',WriteDiskItem( 'ridges texture',  'Aims texture formats' ),
    'basins_texture',WriteDiskItem( 'basins texture',  'Aims texture formats' ),
#    'areas_texture',WriteDiskItem( ' texture',  'Aims texture formats' ),
)


# Default values
def initialization( self ):
    self.linkParameters( 'DPF_texture', 'input_mesh' )
    self.linkParameters( 'pits_texture', 'input_mesh' )
    self.linkParameters( 'noisypits_texture', 'input_mesh' )
    self.linkParameters( 'ridges_texture', 'input_mesh' )
    self.linkParameters( 'basins_texture', 'input_mesh' )
    self.DPF_alpha = 0.03
    self.thresh_dist=20
    self.thresh_ridge=1.5
    self.thresh_area=50
    self.setOptional('mask_texture')



def execution( self, context ):
    #compute the DPF
    #apply the watershed
    print 'will compute the DPF followed by the watershed'
