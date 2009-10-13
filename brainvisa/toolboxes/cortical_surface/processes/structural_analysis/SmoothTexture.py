from neuroProcesses import *
import shfjGlobals     
from soma import aims

name = 'Smooth Texture'
userLevel = 2

signature = Signature(
  'texture', ReadDiskItem('White Curvature Texture', 'Texture'),
  'mesh', ReadDiskItem('Mesh','MESH Mesh'),
  'dt', Float(),
  's', Float(),
  'outfile', WriteDiskItem('Test Texture', 'Texture'))


def initialization( self ):
  self.linkParameters( 'mesh', 'texture' )
  self.linkParameters( 'outfile', 'texture' )
  self.dt=0.01
  self.s=4.0
  
     

def execution( self, context ):
  call_list = ['AimsTextureSmoothing',
    '-m', self.mesh,
    '-i', self.texture,
    '-o', self.outfile,
    '-t', self.dt,
    '-s', self.s]


  context.write("Smoothing the texture (see BrainVISA log for more details...)")
  apply(context.system, call_list)
  
  
  context.write("Rescaling...")
  context.write(str(self.outfile))
  texreader = aims.Reader()
  texture = texreader.read(str(self.outfile))
  
  imin = 60000000.0
  imax = -60000000.0
  for i in xrange(len(texture[0])):
    if (texture[0][i] > imax):
      imax = texture[0][i]
    if (texture[0][i] < imin):
      imin = texture[0][i]
  for i in xrange(len(texture[0])):
    texture[0][i] = (texture[0][i]-imin)/(imax-imin) *2 -1
  texwriter = aims.Writer()
  texwriter.write(texture, str(self.outfile))
  
  context.write("Finished")