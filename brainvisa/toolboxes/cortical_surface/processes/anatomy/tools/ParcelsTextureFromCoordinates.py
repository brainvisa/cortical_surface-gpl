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


from __future__ import absolute_import
def validation():
    try:
        import brainvisa.cortical_surface.parameterization.mapping
    except:
        raise ValidationError( 'brainvisa.cortical_surface.parameterization.mapping module can not be imported.' )


from brainvisa.processes import *
import sigraph
from soma import aims
from brainvisa.cortical_surface.parameterization import model as md
import numpy as np

try:
    from brainvisa.cortical_surface.parameterization.mapping import parcelsFromCoordinates
except:
    pass

name = 'Parcels Texture From Coordinates'
userLevel = 0

# def validation():
#     anatomist.validation()
    
signature = Signature(
    'latitude', ReadDiskItem('Latitude coordinate texture',
                             'aims Texture formats'),
    'longitude', ReadDiskItem('Longitude coordinate texture',
                              'aims Texture formats'),
    'side', Choice('left', 'right'),
    'model_file',ReadDiskItem('HipHop Model', 'Text File'),
    'between_poles_parcel_anterior_coord', Float(),
    'between_poles_parcel_posterior_coord', Float(),
    'texture_model_parcels',
    WriteDiskItem('hemisphere model parcellation texture',
                  'aims Texture formats',
                  requiredAttributes={'regularized': 'false'}),
    'texture_marsAtlas_parcels',
    WriteDiskItem('hemisphere marsAtlas parcellation texture',
                  'aims Texture formats',
                  requiredAttributes={'regularized': 'false'}),
    'texture_lobes_parcels',
    WriteDiskItem('hemisphere lobes parcellation texture',
                  'aims Texture formats',
                  requiredAttributes={'regularized': 'false'}),
    'texture_gyri_parcels',
    WriteDiskItem('hemisphere gyrus parcellation texture',
                  'aims Texture formats',
                  requiredAttributes={'regularized': 'false'}),
)


def initialization( self ):
    def linkSide( proc, dummy ):
        if proc.latitude is not None:
            side = proc.latitude.get( 'side' )
            return side
    self.linkParameters( 'longitude','latitude')
    self.linkParameters('side', 'latitude', linkSide)
    self.between_poles_parcel_anterior_coord = 130.0
    self.between_poles_parcel_posterior_coord = 250.0
    self.linkParameters( 'model_file', 'latitude' )
    self.linkParameters( 'texture_model_parcels', 'latitude')
    self.linkParameters( 'texture_marsAtlas_parcels', 'latitude')
    self.linkParameters( 'texture_lobes_parcels', 'latitude')
    self.linkParameters( 'texture_gyri_parcels', 'latitude')


def execution( self, context ):
       
    re = aims.Reader()
    ws = aims.Writer()
    
    latitude_texture = re.read(self.latitude.fullPath())
    longitude_texture = re.read(self.longitude.fullPath())
    
    model = md.Model().read(self.model_file.fullPath())
#    context.write('------------------- model used -------------------')
#    for line in model.printArgs().splitlines():
#        context.write(line)
    context.warning('NOTE: in every parcellation textures, '
                    'CINGULAR POLE = 0 and PARCEL BETWEEN POLES = 255')

    (tex_parcels, nb_parcels) = \
        parcelsFromCoordinates(latitude_texture[0].arraydata(),
                               longitude_texture[0].arraydata(),
                               model,
                               'model',
                               self.between_poles_parcel_anterior_coord,
                               self.between_poles_parcel_posterior_coord)
#   context.write('----------------------------------------------------------')
    context.write('number of parcels created for the model-based parcellation '
                  '(including the cingular pole) :')
    context.write(nb_parcels)
    if self.side == 'right':
        tex_parcels = tex_parcels + 100
        tex_parcels[tex_parcels == 100] = 0  # cingular pole
        tex_parcels[tex_parcels == 355] = 255  # path between poles

    context.write('Writing texture')
    aims_tex_parcels = aims.TimeTexture_S32()
    aims_tex_parcels[0].assign(tex_parcels.astype(np.int32))
    ws.write(aims_tex_parcels, self.texture_model_parcels.fullPath())

    (tex_parcels, nb_parcels) = \
        parcelsFromCoordinates(latitude_texture[0].arraydata(),
                               longitude_texture[0].arraydata(),
                               model,
                               'marsAtlas',
                               self.between_poles_parcel_anterior_coord,
                               self.between_poles_parcel_posterior_coord)
#   context.write('----------------------------------------------------------')
    context.write('number of parcels created for the marsAtlas parcellation '
                  '(including the cingular pole) :')
    context.write(nb_parcels)
    if self.side == 'right':
        tex_parcels = tex_parcels + 100
        tex_parcels[tex_parcels == 100] = 0  # cingular pole
        tex_parcels[tex_parcels == 355] = 255  # path between poles

    context.write('Writing texture')
    aims_tex_parcels = aims.TimeTexture_S32()
    aims_tex_parcels[0].assign(tex_parcels.astype(np.int32))
    ws.write(aims_tex_parcels, self.texture_marsAtlas_parcels.fullPath())

    (tex_parcels, nb_parcels) = \
        parcelsFromCoordinates(latitude_texture[0].arraydata(),
                               longitude_texture[0].arraydata(),
                               model,
                               'lobes',
                               self.between_poles_parcel_anterior_coord,
                               self.between_poles_parcel_posterior_coord)
#   context.write('----------------------------------------------------------')
    context.write('number of parcels created for the lobes parcellation '
                  '(including the cingular pole) :')
    context.write(nb_parcels)
    if self.side == 'right':
        tex_parcels = tex_parcels + 100
        tex_parcels[tex_parcels == 100] = 0  # cingular pole
        tex_parcels[tex_parcels == 355] = 255  # path between poles

    context.write('Writing texture')
    aims_tex_parcels = aims.TimeTexture_S32()
    aims_tex_parcels[0].assign(tex_parcels.astype(np.int32))
    ws.write(aims_tex_parcels, self.texture_lobes_parcels.fullPath())

    (tex_parcels, nb_parcels) = \
        parcelsFromCoordinates(latitude_texture[0].arraydata(),
                               longitude_texture[0].arraydata(),
                               model,
                               'gyri',
                               self.between_poles_parcel_anterior_coord,
                               self.between_poles_parcel_posterior_coord)
#   context.write('----------------------------------------------------------')
    context.write('number of parcels created for the gyri parcellation '
                  '(including the cingular pole) :')
    context.write(nb_parcels)
    if self.side == 'right':
        tex_parcels = tex_parcels + 100
        tex_parcels[tex_parcels == 100] = 0  # cingular pole
        tex_parcels[tex_parcels == 355] = 255  # path between poles

    context.write('Writing texture')
    aims_tex_parcels = aims.TimeTexture_S32()
    aims_tex_parcels[0].assign(tex_parcels.astype(np.int32))
    ws.write(aims_tex_parcels, self.texture_gyri_parcels.fullPath())

    context.write('Done')
