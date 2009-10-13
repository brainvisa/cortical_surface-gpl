from neuroProcesses import *
import shfjGlobals     

name = 'Dilate Texture'
userLevel = 3

signature = Signature(
  'input', ReadDiskItem( 'Texture', 'Texture'),
  'mesh',  ReadDiskItem( 'Mesh', 'MESH Mesh'),
  'output', WriteDiskItem( 'Texture', 'Texture' ),
  'size', Float()
  )

def initialization( self ):
  self.linkParameters( 'mesh', 'input' )

  self.linkParameters( 'output', 'input' )

  self.size=5.0


def execution( self, context ):

  command = [ 'AimsTextureDilation', '-t', self.input, '-i', self.mesh, '-o', self.output, '-s', self.size]

  print command

  context.system( *command )

