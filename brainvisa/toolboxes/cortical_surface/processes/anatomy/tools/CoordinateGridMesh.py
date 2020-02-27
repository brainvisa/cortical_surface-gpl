# Copyright CEA and IFR 49 (2000-2005)
#
#  This software and supporting documentation were developed by
#      CEA/DSV/SHFJ and IFR 49
#      4 place du General Leclerc
#      91401 Orsay cedex
#      France
#
# This software is governed by the CeCILL license version 2 under
# French law and abiding by the rules of distribution of free software.
# You can  use, modify and/or redistribute the software under the
# terms of the CeCILL license version 2 as circulated by CEA, CNRS
# and INRIA at the following URL "http://www.cecill.info".
#
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
from six.moves import range
def validation():
    try:
        from brainvisa.cortical_surface.surface_tools import basic_tools as bsTls
    except:
        raise ValidationError( 'brainvisa.cortical_surface.parameterization.surface_tools module can not be imported.' )


from brainvisa.processes import *
try:
    from brainvisa.cortical_surface.surface_tools import basic_tools as bsTls
    from brainvisa.cortical_surface.parameterization import model as md
    import numpy as np
except:
    pass


name = 'Coordinate Grid Mesh'

userLevel = 0

signature = Signature(
    'mesh', ReadDiskItem( 'Mesh', 'Aims mesh formats'),
    'latitude', ReadDiskItem( 'Latitude coordinate texture','aims Texture formats'),
    'longitude', ReadDiskItem( 'Longitude coordinate texture','aims Texture formats'),
    'grid', WriteDiskItem( 'Coordinate grid', 'Aims mesh formats' ),
    'mode', Choice('model', 'Regular'),
    'model_file',ReadDiskItem( 'HipHop Model', 'Text File'),
    'interval', Integer(),
    'tube_size', Float()
    )

def initialization( self ):
    self.linkParameters( 'latitude', 'mesh' )
    self.linkParameters( 'model_file','mesh' )
    self.linkParameters( 'longitude', 'mesh' )
    self.linkParameters( 'grid', 'mesh' )
    self.tube_size = 0.25
    self.interval = 20
    self.setOptional('model_file')
  
  

def execution( self, context ):


    mesh = aims.read(self.mesh.fullPath())
    tex_lon = aims.read(self.longitude.fullPath())
    tex_lat = aims.read(self.latitude.fullPath())

    if  (self.mode == "model"):
        model = md.Model().read(self.model_file.fullPath())
        (longitude_axis_coords, latitude_axis_coords) = model.axisCoordToDegree()
        latitude_axis_coords = np.append(latitude_axis_coords, [180-model.cingularPoleBoundaryCoord, model.insularPoleBoundaryCoord])
    else:
        longitude_axis_coords = list(range(1, 360, self.interval))
        latitude_axis_coords =  list(range(1, 180, self.interval))
    context.write('isocoordinates to be extracted:')
    context.write(longitude_axis_coords)
    context.write(latitude_axis_coords)
    all_tubes = aims.AimsSurfaceTriangle()
    line = bsTls.meshIsoLine(mesh, tex_lon, longitude_axis_coords)
    all_tubes += bsTls.linesToTubes(line, self.tube_size)
    line = bsTls.meshIsoLine(mesh, tex_lat, latitude_axis_coords)
    all_tubes += bsTls.linesToTubes(line, self.tube_size)
    aims.write(all_tubes, self.grid.fullPath())
    context.write('done')
  # command = [ 'AimsCoordinateGridMesh',
  #               '-m', self.mesh,
  #               '-x', self.latitude,
  #               '-y', self.longitude,
  #               '-o', self.grid,
  #               '-d', self.tube_size]
  #               #'-c', "c" ]
  # if (self.mode == "Constraints"):
  #   command += ['-c', 'c']
  # if (self.mode == "Regular"):
  #   command += ['-c', 'r']
  # if (self.mode == "Sulcus"):
  #   command += ['-c', 's']
  #
  # context.write('Generating grid')
  # context.system(*command)
  # context.write('Finished')

