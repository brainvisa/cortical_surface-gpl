
from neuroProcesses import *
import shfjGlobals   

name = 'Project Texture Onto Atlas'

userLevel = 2


signature = Signature(
    'white_mesh',ReadDiskItem( 'Hemisphere White Mesh' , shfjGlobals.aimsMeshFormats),
    'texture_on_mesh', ReadDiskItem('Texture','Texture'),
    'longitude',ReadDiskItem( 'Longitude coordinate texture', 'Texture'),
    'latitude',ReadDiskItem( 'Latitude coordinate texture', 'Texture'),
    'atlas',ReadDiskItem('Mesh' , shfjGlobals.aimsMeshFormats),
    'atlas_longitude', ReadDiskItem('Texture','Texture'),
    'atlas_latitude', ReadDiskItem('Texture','Texture'),
    'texture_on_atlas', WriteDiskItem('Texture','Texture')
)

def initialization( self ):
    self.linkParameters( 'white_mesh','longitude')
    self.linkParameters( 'white_mesh','latitude')

def execution( self, context ):
     context.system('AimsTextureToAtlas', '-i', self.white_mesh.fullPath(), '-t', self.texture_on_mesh.fullPath(), '-o', self.texture_on_atlas.fullPath(), '-a', self.atlas.fullPath(), '-ix', self.longitude.fullPath(), '-iy', self.latitude.fullPath(), '-ax', self.atlas_longitude.fullPath(), '-ay',  self.atlas_latitude.fullPath(), '-px', 360, '-py', 0)
     context.write('Done')
