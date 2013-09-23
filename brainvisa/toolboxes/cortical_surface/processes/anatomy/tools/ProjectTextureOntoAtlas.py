
from brainvisa.processes import *
from brainvisa.tools import aimsGlobals

name = 'Project Texture Onto Atlas'

userLevel = 2


signature = Signature(
    'white_mesh',ReadDiskItem( 'Hemisphere White Mesh' , aimsGlobals.aimsMeshFormats),
    'texture_on_mesh', ReadDiskItem('Texture','Texture'),
    'longitude',ReadDiskItem( 'Longitude coordinate texture', 'Aims Texture formats'),
    'latitude',ReadDiskItem( 'Latitude coordinate texture', 'Aims Texture formats'),
    'atlas',ReadDiskItem('Mesh' , aimsGlobals.aimsMeshFormats),
    'atlas_longitude', ReadDiskItem('Texture','Aims Texture formats'),
    'atlas_latitude', ReadDiskItem('Texture','Aims Texture formats'),
    'texture_on_atlas', WriteDiskItem('Texture','Aims Texture formats')
)

def initialization( self ):
    self.linkParameters( 'white_mesh','longitude')
    self.linkParameters( 'white_mesh','latitude')

def execution( self, context ):
     context.system('AimsTextureToAtlas', '-i', self.white_mesh.fullPath(), '-t', self.texture_on_mesh.fullPath(), '-o', self.texture_on_atlas.fullPath(), '-a', self.atlas.fullPath(), '-ix', self.longitude.fullPath(), '-iy', self.latitude.fullPath(), '-ax', self.atlas_longitude.fullPath(), '-ay',  self.atlas_latitude.fullPath(), '-px', 360)
     context.write('Done')
