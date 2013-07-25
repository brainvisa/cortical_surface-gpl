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
from soma import aims
import numpy as np

try:
  from brainvisa.cortical_surface.parameterization.mapping import computeCoordinates
except:
  pass

#from brainvisa import anatomist

name = 'CoordinatesFromHiphopMapping'

userLevel = 2

# def validation():
#     anatomist.validation()
    
signature = Signature(
    'input_texture',ReadDiskItem( 'hemisphere Sulcal Lines texture', 'Texture' ),                     
#    'input_texture',ReadDiskItem('Texture', 'Texture'),
    'corresp_indices_texture',ReadDiskItem( 'Rectangular flat indices texture', 'Texture'),
    'boundary_texture',ReadDiskItem( 'Rectangular boundary texture', 'Texture'),
    'output_texture',WriteDiskItem( 'hemisphere Sulcal Lines Rectangular Flat texture', 'Texture' )
    #'output_texture',WriteDiskItem( 'Rectangular flat texture', 'Texture')
)

def initialization( self ):
    self.linkParameters( 'boundary_texture','input_texture')
    self.linkParameters( 'corresp_indices_texture','input_texture')
    self.linkParameters( 'boundary_texture','input_texture')
    self.linkParameters( 'output_texture','input_texture')

    
def execution( self, context ):
    
    lon, lat = computeCoordinates(mesh, neocortex_indices, neoCortex_square_cstr, neoCortex_open_boundary, poles_lat_insula, poles_lat_cingular)

    print 'param de l insula'
    insula_lon = texture2ROI(lon, insula_indices)
    insula_bound_rad = np.pi * (insula_lon[insula_boundary] - 180) / 180 
    circle = np.array([np.cos(insula_bound_rad), np.sin(insula_bound_rad)])
    insula_disk = diskConformalMapping(insula_mesh, insula_boundary, circle)
    insula_lon, insula_lat = coordinatesFromDisk(insula_disk, poles_lat_insula)
    lon[insula_indices] = insula_lon
    lat[insula_indices] = insula_lat

    cingular_lon = texture2ROI(lon, cingular_indices)
    cingular_bound_rad = np.pi * (cingular_lon[cingular_boundary] - 180) / 180 
    circle = np.array([np.cos(cingular_bound_rad), np.sin(cingular_bound_rad)])
    cingular_disk = diskConformalMapping(cingular_mesh, cingular_boundary, circle)
    cingular_lon, cingular_lat = coordinatesFromDisk(cingular_disk, poles_lat_cingular)
    lon[cingular_indices] = cingular_lon
    lat[cingular_indices] = 180 - cingular_lat
    if write_all_steps_to_disk:
        s_tex_cstr_square = aims.TimeTexture_S16()
        s_tex_cstr_square[0].assign(tex_cstr_square)
        for pt in neocortex_poles_path:
            s_tex_cstr_square[0].append(0)
        ws.write(s_tex_cstr_square, '/home/toz/ammon_Lwhite_neocortex_cstr.tex')
        full_sulci.save('/home/toz/ammon_Lwhite_square_full_sulci')
        model.save('/home/toz/model.mesh')
        ws.write(neoCortex_square_cstr, '/home/toz/ammon_Lwhite_square_cstr_'+str(cstrBalance)+'.mesh')
        ws.write(insula_disk, '/home/toz/ammon_Lwhite_insula_disk.mesh')
        ws.write(cingular_disk, '/home/toz/ammon_Lwhite_cingular_disk.mesh')
