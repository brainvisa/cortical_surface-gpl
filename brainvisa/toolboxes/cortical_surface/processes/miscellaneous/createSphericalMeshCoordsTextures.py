from __future__ import absolute_import
from brainvisa.processes import *
from soma import aims
import numpy

name = 'Create Spherical Mesh Coordinates Textures'

signature = Signature(
  'sphere_mesh', ReadDiskItem( 'Mesh', 'Aims mesh formats' ),
  'longitude', WriteDiskItem( 'Longitude Coordinate Texture', 'aims texture formats' ),
  'latitude', WriteDiskItem( 'Latitude Coordinate Texture', 'aims texture formats' ),
)


def execution( self, context ):
  mesh = aims.read( self.sphere_mesh.fullPath() )
  vert = mesh.vertex()
  pos = numpy.array( vert )
  ray = numpy.sqrt( pos[:,0] * pos[:,0] + pos[:,1] * pos[:,1] + pos[:,2] * pos[:,2] )
  plray = numpy.sqrt( pos[:,0] * pos[:,0] + pos[:,1] * pos[:,1] )
  lat = numpy.arccos( pos[:,2] / ray ) / numpy.pi * 180.
  lattex = aims.TimeTexture_FLOAT()
  lattex[0].assign( lat )
  aims.write( lattex, self.latitude.fullPath() )
  lon = numpy.arcsin( pos[:,1] / plray ) / numpy.pi * 180.
  lon[ pos[:,0] > 0 ] = 180. - lon[ pos[:,0] > 0 ]
  lon += 90.
  lontex = aims.TimeTexture_FLOAT()
  lontex[0].assign( lon )
  aims.write( lontex, self.longitude.fullPath() )

