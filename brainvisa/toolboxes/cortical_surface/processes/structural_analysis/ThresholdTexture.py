from neuroProcesses import *
import shfjGlobals     

name = 'Threshold Texture'
userLevel = 3

signature = Signature(
  'input', ReadDiskItem( 'Texture', 'Texture'),
  'output', WriteDiskItem( 'Texture', 'Texture' ),
  'mode', Choice ( ( 'less than', 'lt' ), 
                   ( 'less or equal', 'le' ), 
                   ( 'greater than', 'gt' ),
                   ( 'greater or equal', 'ge' ),
                   ( 'equal', 'eq' ), 
                   ( 'different', 'di' ),
                   ( 'between, include both bounds', 'be' ), 
                   ( 'between, exclude lower bound', 'beel' ),
                   ( 'between, exclude higher bound', 'beeh' ),
                   ( 'between, exclude both bound', 'bee' ),
                   ( 'outside, exclude both bounds', 'ou' ),
                   ( 'outside, include lower bound', 'ouil' ),
                   ( 'outside, include higher bound', 'ouih' ),
                   ( 'outside, include both bound', 'oui' ), ),
  'threshold1', Float(),
  'threshold2', Float(),
  'binary', Boolean(),
  )

def initialization( self ):
  self.linkParameters( 'output', 'input' )
  self.setOptional('threshold2', 'binary')
  self.binary = 1
  self.threshold1=-0.5
  self.mode='le'


def execution( self, context ):

  command = [ 'AimsTextureThreshold', '-i', self.input, '-o', self.output, '-m', self.mode, '-t', self.threshold1 ]

  if self.threshold2 is not None :
    command += [ '-u', self.threshold2]

  if self.binary:
    command += [ '-b']

  print command

  context.system( *command )

