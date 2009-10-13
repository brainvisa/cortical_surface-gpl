from neuroProcesses import *
import shfjGlobals     

name = 'Texture Connected Components'
userLevel = 3

signature = Signature(
  'input', ReadDiskItem( 'Texture', 'Texture'),
  'mesh',  ReadDiskItem( 'Mesh', 'MESH Mesh'),
  'output', WriteDiskItem( 'Texture', 'Texture' )
  )

def initialization( self ):
  self.linkParameters( 'mesh', 'input' )
  self.linkParameters( 'output', 'input' )



def execution( self, context ):

  command = [ 'AimsMeshConnectedComponent', '-t', self.input, '-i', self.mesh, '-o', self.output, '-m', '0', '-T', '0.05']

  print command

  context.system( *command )

