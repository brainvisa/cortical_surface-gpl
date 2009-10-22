from neuroProcesses import *
import shfjGlobals
from soma import aims


name = 'Create Results Labels Textures'
userLevel = 2

signature = Signature(
  'labeled_primalsketches', ListOf(ReadDiskItem('Primal Sketch', 'Graph and data')),
  'meshes', ListOf(ReadDiskItem('Hemisphere White Mesh', 'MESH Mesh')),
  'textures', ListOf(WriteDiskItem('Labeled Functional Blobs Texture', 'Texture'))
    )

def initialization( self ):
  pass
     

def execution( self, context ):
  assert(len(self.labeled_primalsketches) == len(self.meshes))
  assert(len(self.meshes) == len(self.textures))
  
  for j in xrange(len(self.labeled_primalsketches)):
    context.write(self.labeled_primalsketches[j])
    reader = aims.Reader()
    graph = reader.read(str(self.labeled_primalsketches[j]))
    meshreader = aims.Reader()
    mesh = meshreader.read(str(self.meshes[j]))

    texture = aims.TimeTexture_FLOAT(1,mesh.vertex().size())
    for i in xrange(len(texture[0])):
      texture[0][i] = 0
    context.write(graph.order())
    #context.write(graph['sujet'])
    nbblobsnonnuls=0

    for v in graph.vertices():
      nodes_list = v['nodes_list']
      label = int(v['label'])
      for node in nodes_list:
        texture[0][node] = 1
    for v in graph.vertices():
      nodes_list = v['nodes_list']
      label = int(v['label'])
      if (label>0):
        print str(label) + " " + str(v['tValue'])  #+ " " + str(v['lifeTime']) + " " + str(v['t'])
        nbblobsnonnuls = nbblobsnonnuls + 1
        for node in nodes_list:
          texture[0][node] =int(v['label'])+1
    context.write("nb blobs non nuls : " + str(nbblobsnonnuls))
    

    context.write("Texture " + str(j))
    context.write("Writing texture n_" + str(j) +"...")
    wtex = aims.Writer()
    wtex.write(texture, str(self.textures[j]))