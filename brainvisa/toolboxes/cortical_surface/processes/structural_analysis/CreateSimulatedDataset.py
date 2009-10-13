from neuroProcesses import *
import shfjGlobals     

name = 'Create Simulated Data Sets'
userLevel = 2

signature = Signature(
  'mesh', ReadDiskItem('Mesh', 'MESH Mesh'),
  'texture', WriteDiskItem('Texture', 'Texture'),
  'noise', Float(),
  'variab_loc', Float(),
  'variab_int', Float(),
  'intensities', ListOf(Integer()),
  'nodes', ListOf(Integer()),
  'smoothing', Integer()
  )

def initialization( self ):
  self.noise = 0.0
  self.variab_loc = 19.0
  self.variab_int = 0.0
  self.intensities = [6,6,6,6,6]
  self.nodes = [11238,20573,11710,4253,20866]
  self.smoothing = 128
     

def execution( self, context ):
  assert(len(self.nodes) == len(self.intensities))
    #context.write( 'ok')
  intensities_str = ""
  for i in xrange(len(self.intensities)-1):
    intensities_str += str(self.intensities[i]) + "' '"
  intensities_str += str(self.intensities[len(self.intensities)-1])
  nodes_str = ""
  for x in self.nodes:
    nodes_str += str(x) + " " 
  call_list = ['surfTexActivationSimulation',
    '-m', self.mesh,
    '-o', self.texture,
    '-n', self.noise,
    '-l', self.variab_loc,
    '-i', self.variab_int,
    '-I'] + self.intensities + ['-F'] + self.nodes + ['-s', self.smoothing]

  context.write('Generating a simulated texture')
  apply( context.system, call_list )
  context.write('Finished')