from neuroProcesses import *
import shfjGlobals
from soma import aims


name = 'Create Average Texture'
userLevel = 2

signature = Signature(
  'textures', ReadDiskItem('Texture', 'Texture'),
  'average', WriteDiskItem('Texture', 'Texture')
    )

def initialization( self ):
  pass
  
     

def execution( self, context ):
  #for i in xrange(
  textures = aims.read(str(self.textures))

  average = aims.TimeTexture_S16(1,len(textures[0]))
  for j in xrange(len(textures[0])):
    p=0
    for i in xrange(textures.size()):
      if (textures[i][j]!=0):
        p=p+1
    average[0][j] = p
  context.write("Writing texture average")
  wtex = aims.Writer()
  wtex.write(average, str(self.average) )