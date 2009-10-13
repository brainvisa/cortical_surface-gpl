from neuroProcesses import *
import shfjGlobals     

name = 'Convert to FLOAT Texture'
userLevel = 3

signature = Signature(
  'input', ReadDiskItem( 'Texture', 'Texture'),
  'output', WriteDiskItem( 'Texture', 'Texture'),
  )

def initialization( self ):
  self.linkParameters( 'output', 'input' )



def execution( self, context ):
  
  command = [ 'AimsFileConvert', '-i', self.input, '-o', self.output, '-t', 'FLOAT' ]

  print command

  context.system( *command )

