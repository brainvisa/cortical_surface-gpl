from neuroProcesses import *
import shfjGlobals     

name = 'Create Scale Space'
userLevel = 3

signature = Signature(
  'input', ReadDiskItem( 'Texture', 'Texture'),
  'mesh',  ReadDiskItem( 'Mesh', 'MESH Mesh'),
  'output', WriteDiskItem( 'Scale Space White Curvature Texture', 'Texture' ),
  'scalemax', Float(),
  'dt', Float()
  )

def initialization( self ):
  self.linkParameters( 'mesh', 'input' )

  self.linkParameters( 'output', 'input' )

  self.scalemax=64
  self.dt=0.01


def execution( self, context ):

  command = [ 'AimsTextureScaleSpace', '-i1', self.input, '-i2', self.mesh, '-o', self.output, '-tm', self.scalemax, '-dt', self.dt]

  print command

  context.system( *command )

