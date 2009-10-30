from neuroProcesses import *
import shfjGlobals
from soma import aims


name = '5 - Create Results Labels Textures'
userLevel = 2

signature = Signature(
  'labeled_primalsketches', ListOf(ReadDiskItem('Primal Sketch', 'Graph and data')),
  'meshes', ListOf(ReadDiskItem('Hemisphere White Mesh', 'MESH Mesh')),
  'textures', ListOf(WriteDiskItem('Labeled Functional Blobs Texture', 'Texture')),
  'mode', Choice('all','more than <threshold> subjects','over 95.% significance'),
  'threshold', Integer()  
    )

def initialization( self ):
  self.setOptional('threshold')
  
     

def execution( self, context ):
  assert(len(self.labeled_primalsketches) == len(self.meshes))
  assert(len(self.meshes) == len(self.textures))
  
  for j in xrange(len(self.labeled_primalsketches)):
    context.write(self.labeled_primalsketches[j])
    reader = aims.Reader()
    graph = reader.read(str(self.labeled_primalsketches[j]))
    meshreader = aims.Reader()
    mesh = meshreader.read(str(self.meshes[j]))

    if (self.mode == 'all'):
      texture = aims.TimeTexture_FLOAT(len(self.labeled_primalsketches)+1,mesh.vertex().size())
    if (self.mode == 'more than <threshold> subjects'):
      texture = aims.TimeTexture_FLOAT(1,mesh.vertex().size())
    for k in xrange(texture.size()):
      for i in xrange(len(texture[k])):
        texture[k][i] = 0
    print "ok"
    context.write(graph.order())
    #context.write(graph['sujet'])
    nbblobsnonnuls=0
  
    for v in graph.vertices():
      nodes_list = v['nodes_list']
      label = int(v['label'])
      for node in nodes_list:
        texture[0][node] = 0
    for v in graph.vertices():
      nodes_list = v['nodes_list']
      label = int(v['label'])
      lon = int(v['label_occur_number'])
      ok = 0
      if (lon>len(self.labeled_primalsketches)):
        lon = len(self.labeled_primalsketches)
      if (self.mode == 'more than <threshold> subjects'):
        if (lon > int(self.threshold)):
          lon=0
          ok = 1
      if (self.mode == 'all'):
        ok = 1
      if (ok == 1 and label>0):
        for node in nodes_list:
          texture[lon][node] =int(v['label_occur_number'])
        #print str(label) + " " + str(v['tValue'])  #+ " " + str(v['lifeTime']) + " " + str(v['t'])
        nbblobsnonnuls = nbblobsnonnuls + 1

    context.write("nb blobs non nuls : " + str(nbblobsnonnuls))
    

    context.write("Texture " + str(j))
    context.write("Writing texture n_" + str(j) +"...")
    wtex = aims.Writer()
    wtex.write(texture, str(self.textures[j]))