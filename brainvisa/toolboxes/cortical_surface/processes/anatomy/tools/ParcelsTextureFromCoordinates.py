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

userLevel = 1

# def validation():
#     anatomist.validation()
    
signature = Signature(
                      
    'latitude',ReadDiskItem( 'Latitude coordinate texture','Texture'),
    'side', Choice('left', 'right'),
    'longitude',ReadDiskItem( 'Longitude coordinate texture','Texture'),
#    'model', ReadDiskItem(  ),
    'texture_parcels', WriteDiskItem('hemisphere gyri parcellation texture','Texture', requiredAttributes={ 'regularized': 'false' }),
    'model_file',ReadDiskItem( 'HipHop Model', 'Text File'),

)

def initialization( self ):
    def linkSide( proc, dummy ):
        if proc.latitude is not None:
            side = proc.latitude.get( 'side' )
            return side
    self.linkParameters( 'side','latitude', linkSide )
    self.linkParameters( 'longitude','latitude')
    self.linkParameters( 'texture_parcels','latitude')
    self.linkParameters( 'model_file', 'latitude' )
    
def execution( self, context ):
       
    re = aims.Reader()
    ws = aims.Writer()
    
    latitude_texture = re.read(self.latitude.fullPath())
    longitude_texture = re.read(self.longitude.fullPath())
    
    model = md.Model().read(self.model_file.fullPath())
    context.write('------------------- model used -------------------')
    for line in model.printArgs().splitlines():
        context.write(line)

    tex_parcels = parcelsFromCoordinates(latitude_texture[0].arraydata(), longitude_texture[0].arraydata(), model)
        
    if self.side =='left':
        tex_parcels = tex_parcels + 100
    context.write('Writing texture')
    aims_tex_parcels = aims.TimeTexture_S16()
    aims_tex_parcels[0].assign(tex_parcels)
    ws.write(aims_tex_parcels, self.texture_parcels.fullPath())

    context.write('Done')
            
      
