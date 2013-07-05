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
    import brainvisa.cortical_surface.parameterization.mapping
  except:
    raise ValidationError( 'brainvisa.cortical_surface.parameterization.mapping module can not be imported.' )
  
from brainvisa.processes import *
import shfjGlobals   
import sigraph
from soma import aims
from brainvisa.cortical_surface.parameterization import model as md

try :
  from brainvisa.cortical_surface.parameterization.mapping import parcelsFromCoordinates
except:
  pass

name = 'ParcelsTextureFromCoordinates'

userLevel = 2

# def validation():
#     anatomist.validation()
    
signature = Signature(
                      
    'Side', Choice('right', 'left'),
    'latitude',ReadDiskItem( 'Latitude coordinate texture','Texture'),
    'longitude',ReadDiskItem( 'Longitude coordinate texture','Texture'),
#    'model', ReadDiskItem(  ),
    'texture_parcels', WriteDiskItem('hemisphere gyri parcellation texture','Texture')

)

def initialization( self ):
     self.linkParameters( 'longitude','latitude')
#     self.linkParameters( 'remeshed_mesh','white_mesh')
    
def execution( self, context ):
       
    re = aims.Reader()
    ws = aims.Writer()
    
    latitude_texture = re.read(self.latitude.fullPath())
    longitude_texture = re.read(self.longitude.fullPath())
    
    context.write('Using the default model')    
    default_model = md.Model()
    default_model.printArgs()
    default_model.saveToFile('/home/toz/model_default.txt')
    tex_parcels = parcelsFromCoordinates(latitude_texture[0].arraydata(), longitude_texture[0].arraydata(), default_model)
        
    if self.Side =='left':
        tex_parcels = tex_parcels + 100
    context.write('Writing texture')
    aims_tex_parcels = aims.TimeTexture_S16()
    aims_tex_parcels[0].assign(tex_parcels)
    ws.write(aims_tex_parcels, self.texture_parcels.fullPath())

    context.write('Done')
            
      
