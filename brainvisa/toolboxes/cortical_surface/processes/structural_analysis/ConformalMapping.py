from neuroProcesses import *
import shfjGlobals     

name = 'Conformal Mapping'
userLevel = 3

signature = Signature(
  'input', ReadDiskItem( 'Hemisphere White Mesh', 'MESH Mesh'),
  #'mesh',  ReadDiskItem( 'Mesh', 'MESH Mesh'),
  'output', WriteDiskItem( 'Conformal White Mesh', 'MESH Mesh'),
  'latitude', WriteDiskItem( 'Conformal latitude texture', 'Texture'),
  'longitude', WriteDiskItem( 'Conformal longitude texture', 'Texture')
  )

def initialization( self ):
  self.linkParameters( 'output', 'input' )
  self.linkParameters( 'latitude', 'input' )
  self.linkParameters( 'longitude', 'input' )
  #self.linkParameters( 'output', 'input' )
  #pass


def execution( self, context ):

  command = [ 'AimsConformalMapping', '-i', self.input, '-o', self.output, '-l', self.latitude, '-L', self.longitude]

  print command

  context.system( *command )

