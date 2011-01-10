from neuroProcesses import *
import shfjGlobals
from soma import aims


name = '5 - Create Results Labels Textures'
userLevel = 2

signature = Signature(
    'labelled_primalsketches', ListOf(ReadDiskItem('Primal Sketch', 'Graph and data')),
    'textures', ListOf(WriteDiskItem('Labelled Functional Blobs Texture', 'Texture')),
    'mode', Choice('all','more than <threshold> subjects','over 95.% significance'),
    'threshold', Integer()  
    )

def initialization( self ):
  self.setOptional('threshold')
  self.linkParameters("textures", "labelled_primalsketches")

def execution( self, context ):
    graphs = [aims.read(each) for each in self.labelled_primalsketches]
    meshes = [ each['mesh'] for each in graphs ]

    for j in xrange(len(self.labelled_primalsketches)):
        context.write(self.labelled_primalsketches[j])
        if (self.mode == 'all'):
            texture = aims.TimeTexture_S16( len(self.labelled_primalsketches) + 1, meshes[j].vertex().size() )
        elif (self.mode == 'more than <threshold> subjects' or self.mode == 'over 95.% significance'):
            texture = aims.TimeTexture_S16( 1, meshes[j].vertex().size() )
        for k in xrange( texture.size() ):
            for i in xrange( len(texture[k]) ):
                texture[k][i] = 0
        context.write(graphs[j].order())

        nbblobsnonnuls = 0
  
        for v in graphs[j].vertices():
            nodes_list = v['nodes']
            label = int(v['label'])
            for node in nodes_list:
                texture[0][node] = 0
        for v in graphs[j].vertices():
            if ( set(['label_occur_number', 'significance']).issubset(set(v.keys())) ) :
                nodes_list = v['nodes']
                label = int(v['label'])
                lon = int(v['label_occur_number'])
                signif = float(v['significance'])
                ok = 0
                if ( lon > len(self.labelled_primalsketches) ):
                    lon = len(self.labelled_primalsketches)
                if ( self.mode == 'more than <threshold> subjects' ):
                    if ( lon > int(self.threshold) ):
                        lon = 0
                        ok = 1
                elif ( self.mode == 'all' ):
                    ok = 1
                elif ( self.mode == 'over 95.% significance' ):
                    if (signif > 95.0):
                        lon = 0
                        ok = 1
                if ( ok == 1 and label > 0 ):
                    for node in nodes_list:
                        texture[lon][node] = int ( v['label_occur_number'] )
                #print str(label) + " " + str(v['tValue'])  #+ " " + str(v['lifeTime']) + " " + str(v['t'])
                nbblobsnonnuls = nbblobsnonnuls + 1

        context.write("nb blobs non nuls : " + str(nbblobsnonnuls))
        context.write("Texture " + str(j))
        context.write("Writing texture n_" + str(j) +"...")
        wtex = aims.Writer()
        wtex.write(texture, str(self.textures[j]))