from __future__ import print_function

from brainvisa.cortical_surface.structural.blobManip import Node
from brainvisa.cortical_surface.shell import db
from soma import aims
import os

def Graph ( db_path, nodes, hemis_side = 'L', data_type = 'curv' ):
    ''' Creates an Aims Graph given a list of nodes (vertex indices + subjects).
    In this function a node is a Node. Allows to build a quick generic Aims graph
    out from a nodes vector.'''
    #from soma import aims
    import numpy as np 
    sujets = db.getSubjects( db_path )
    pathes = db.getPathes( db_path, hemis_side, data_type )
    meshes = {}
    textures = {}
    for sujet in sujets:
        meshes[sujet] = aims.read(pathes[sujet]['mesh'])
        textures[sujet] = aims.read(pathes[sujet]['tex'])
        
    graph = aims.Graph('BlobsArg')
    graph['boundingbox_min'] = [0.0, 0.0, 0.0]
    graph['boundingbox_max'] = [10.0, 10.0, 10.0]
    graph['voxel_size'] = [1.0, 1.0, 1.0]
    graph['filename_base'] = "*"
    
    sujets = [n.subject for n in nodes]
    set_sujets = list( set(sujets) )
    if ( len(set_sujets) == 1 ):
        sujet = set_sujets[0]
        graph['mesh'] = pathes[sujet]['mesh']
        graph['texture'] = pathes[sujet]['tex']
        graph['subject'] = sujet
        if (os.path.exists(pathes[sujet]['lat'])):
            graph['latitude'] = pathes[sujet]['lat']
        if (os.path.exists(pathes[sujet]['lon'])):
            graph['longitude'] = pathes[sujet]['lon']
    else :
        graph['meshes'] = [pathes[sujet]['mesh'] for sujet in set_sujets]
        graph['textures'] = [pathes[sujet]['tex'] for sujet in set_sujets]
        graph['subjects'] = set_sujets
        if (os.path.exists(pathes[set_sujets[0]]['lat'])):
            graph['latitudes'] = [ pathes[sujet]['lat'] for sujet in set_sujets]
        if (os.path.exists(pathes[set_sujets[0]]['lon'])):
            graph['longitude'] = [ pathes[sujet]['lon'] for sujet in set_sujets]
    
    for i,n in enumerate(nodes) :
        v = graph.addVertex ('glb')
        v['subject'] = n.subject        
        v['nodes'] = n.nodes
        v['scale'] = n.scale
        v['name'] = str(sujets.index(n.subject) + 1)        
        v['t'] = n.t
        v['x'] = list(np.ones(len(n.nodes)))
        v['y'] = list(np.ones(len(n.nodes)))
        #assert( len(n.nodes) == 1 )
        c = meshes[ n.subject ].vertex()[ n.maxnode ]
        mesh = aims.SurfaceGenerator.sphere ( c, 0.3, 10, True )
        aims.GraphManip.storeAims ( graph, v, 'glb', mesh )
        
    #aims.write( graph, graph_path )
    return graph

def BuildAimsGroupGraph ( blobs, cliques, graphs ):
    ''' Creates an Aims Group Graph given a list of nodes (vertex indices + subjects)
    and cliques previously computed by ComputeBliques. In this function a node is a Node.
    The Graphs are given so as to store file paths in the group graph.'''
    #from soma import aims
    graph = aims.Graph('BlobsArg')
    graph['boundingbox_min'] = [0.0, 0.0, 0.0]
    graph['boundingbox_max'] = [10.0, 10.0, 10.0]
    graph['voxel_size'] = [3.0, 3.0, 3.0]
    graph['filename_base'] = "*"
    sujets = blobs.keys()
    
    graph['subjects'] = sujets

    #FIXME: add the filepaths from the graphs
    
    glb_indices = {}
    index = 0
    for sujet in sujets : 
        for each in blobs[sujet] :
            v = graph.addVertex ('ssb')
            for i in ['t', 'subject', 'tmin', 'tmax', 'index']:
                v[i] = each[i]
            bucket = getAimsBucketFromBucket(each['bucket'])
            aims.GraphManip.storeAims( graph, v, 'aims_ssb', bucket )
            glb_indices[each['index']] = v
    overlaps = []
    for each in cliques:
        e = graph.addEdge ( glb_indices[each.node1['index']],
                            glb_indices[each.node2['index']] , 'b2b')
        e['similarity'] = each.overlap
        overlaps.append(each.overlap)
    import numpy as np
    print(min(overlaps), max(overlaps), np.mean(overlaps))
    return graph

def groupGraphInfo ( graph ) :
    '''Returns distances, overlaps and activation ('t' measures) contained in a
    previously built groupgraph'''
    how_many_similarity_cliques = 0
    distances = []
    overlaps = []
    measures = []
    for e in graph.edges() :
        if e.getSyntax() == 'b2b' :
            how_many_similarity_cliques += 1
            if e.has_key('distance') and e['distance'] > 0.0 :
                distances.append(e['distance'])
            elif e.has_key('similarity') and e['similarity'] > 0.0 :
                overlaps.append ( e['similarity'] )
    for v in graph.vertices() :
        if v.getSyntax() == 'ssb' :
            measures.append ( v['t'] )
    return distances, overlaps, measures

def getBucketFromVertex ( v ) :
    '''Returns the bucket associated with a given vertex plus its bounding box
    (minimum point and maximum point)'''
    if v.getSyntax() == 'glb' :
        bucketmap = v['aims_glb']
    elif v.getSyntax() == 'ssb' :
        bucketmap = v['aims_ssb']
    print(v.getSyntax(), v.keys())
       
    bucket = {}
    bucket['voxel_list'] = []
    assert(len(bucketmap[0].keys())>0)
    for each in bucketmap[0].keys() :
        bucket['voxel_list'].append([int(each[x]) for x in xrange(3)])
    assert(len(bucket['voxel_list']) > 0 )
    bucket['voxel_size'] = [bucketmap.sizeX(), bucketmap.sizeY(), bucketmap.sizeZ(), bucketmap.sizeT()]
    #assert(bucket['voxel_size'] == [3.0,3.0,3.0,1.0]
    maxpoints = [ int(max([each[x] for each in bucket['voxel_list']])) for x in xrange(3)]
    minpoints = [ int(min([each[x] for each in bucket['voxel_list']])) for x in xrange(3)]
    return bucket, minpoints, maxpoints
    
def getAimsBucketFromBucket ( b ) :
    bucketMap = aims.BucketMap_VOID()
    bucket = bucketMap[0]
    for each in b['voxel_list']:
        bucket[each] = 1
    assert(len(bucket.keys())>0)
    bucketMap.setSizeXYZT (*b['voxel_size'])
    return bucketMap
    
def getRepresentationFromSSB ( ssb ) :
    glb = []
    for n in ssb.neighbours():
        if n.getSyntax() == 'glb':
            node = Node(n)
            node.bucket = getBucketFromVertex(n)
            glb.append ( node )
    scales = list(set([each.scale for each in glb]))
    for each in glb:
        if each.scale == scales[len(scales)/3]:
            print(each.scale)
            return getAimsBucketFromBucket (each.bucket)
            
    assert(False)
    return None

def getSSBFromGraph ( graph ) :
    ssb = []
    print('beware the result is a list of Vertex')
    for v in graph.vertices():
        if v.getSyntax() == 'ssb':
            ssb.append(v)
    return ssb

def AddSSBBuckets ( graph ) :
    for v in graph.vertices():
        if v.getSyntax() == 'ssb':
            glb = []
            for n in v.neighbours():
                if n.getSyntax() == 'glb':
                    node = Node()
                    node.defineFromVertex(n)
                    bucketmap = n['aims_glb']
                    node.bucket = {}
                    node.bucket['voxel_list'] = []
                    assert( len(bucketmap[0].keys()) > 0 )
                    for each in bucketmap[0].keys() :
                        node.bucket['voxel_list'].append(each)
                    assert(len(node.bucket['voxel_list']) > 0 )
                    node.bucket['voxel_size'] = [bucketmap.sizeX(), bucketmap.sizeY(), bucketmap.sizeZ(), bucketmap.sizeT()]
                    glb.append ( node )
            scales = list(set([each.scale for each in glb]))
            for each in glb:
                if each.scale == scales[len(scales)/3]:
                    bucketMap = getAimsBucketFromBucket(each.bucket)
                    
            assert( len(bucketMap[0].keys())>0)
            aims.GraphManip.storeAims( graph, v, 'aims_ssb', aims.rc_ptr_BucketMap_VOID(bucketMap) )

            assert(v.has_key('aims_ssb'))

def ScaleSpaceBlobsFromOneSubject ( graph, sujet ) :
    blobs = []
    print('Adding bucket to every ssblob from the graph...')
    AddSSBBuckets( graph )
    for v in graph.vertices():
        if v.getSyntax() == 'ssb':
            n = {}
            for each in ['t', 'subject', 'tmin', 'tmax', 'index']:
                n[each] = v[each]
            n['bucket'], n['bbmin'], n['bbmax'] = getBucketFromVertex(v)
            print(n['bbmin'], n['bbmax'])
            blobs.append(n)
    return blobs