#!/usr/bin/env python

from __future__ import print_function

from __future__ import absolute_import
import os,sys
from soma import aims
import numpy as np
import math
import brainvisa.cortical_surface.structural.blobManip as oo
from brainvisa.cortical_surface.multiprocessing.mproc import MultiProcExecute
from brainvisa.cortical_surface.shell.inputParameters import *
import six
from six.moves import input



def getfMRIPathes ( db, contrast, subjects = None ) :
    """Builds and returns the pathes according to the data base path and a
        contrast"""
    pathes = {}
    if subjects is None :
        subjects = oo.getSubjects ( db )
    for sujet in subjects :
        pathes[sujet] = {}
        pathes[sujet]['fMRI'] = os.path.join ( db, sujet, 'fMRI/default_acquisition/loc1/spm_analysis/spm_analysis_Norm_S/spmT_00%i.img' % contrast )
        pathes[sujet]['mask'] = os.path.join ( db, sujet, 'fMRI/default_acquisition/loc1/spm_analysis/spm_analysis_Norm_S/mask.img' )
        pathes[sujet]['primal'] = os.path.join ( db, sujet, 't1mri/t1/default_analysis/segmentation/%s_spmT_00%i_primal.arg' % (sujet, contrast) )
        pathes[sujet]['scale'] = os.path.join ( db, sujet, 't1mri/t1/default_analysis/segmentation/%s_spmT_00%i_ss.ima' % (sujet, contrast) )

    return pathes

def computeOverlap ( blob1, blob2 ) :
    bucket1 = blob1['bucket']['voxel_list']
    bucket2 = blob2['bucket']['voxel_list']
    intersection = []
    div = len(bucket1) + len(bucket2)
    for each1 in bucket1:
        for each2 in bucket2:
            if (each1[0] == each2[0] and each1[1] == each2[1] and each1[2] == each2[2]):
                intersection.append(each1)
    return float( 2.0 * float(len(intersection)) / float(div) )

def convertCliques( rawcliques ) :
    cliques = []
    for each in rawcliques:
        c = Clique()
        c.define(each)
        cliques.append(c)
    return cliques

class Clique:
    def __init__( self ):
        pass
    def define ( self, c ):
        self.node1 = c['node1']
        self.node2 = c['node2']
        self.overlap = c['overlap']

def getBucketFromVertex ( v ):
    if v.getSyntax() == 'glb':        
        bucketmap = v['aims_glb']
    elif v.getSyntax() == 'ssb' :
        bucketmap = v['aims_ssb']
    print(v.getSyntax(), list(v.keys()))
       
    bucket = {}
    bucket['voxel_list'] = []
    assert(len(list(bucketmap[0].keys()))>0)
    for each in bucketmap[0].keys() :
        bucket['voxel_list'].append([int(each[x]) for x in six.moves.xrange(3)])
    assert(len(bucket['voxel_list']) > 0 )
    bucket['voxel_size'] = [bucketmap.sizeX(), bucketmap.sizeY(), bucketmap.sizeZ(), bucketmap.sizeT()]

    maxpoints = [ int(max([each[x] for each in bucket['voxel_list']])) for x in six.moves.xrange(3)]
    minpoints = [ int(min([each[x] for each in bucket['voxel_list']])) for x in six.moves.xrange(3)]
    return bucket, minpoints, maxpoints
    
def getAimsBucketFromBucket ( b ):
    bucketMap = aims.BucketMap_VOID()
    bucket = bucketMap[0]
    for each in b['voxel_list']:
        bucket[each] = 1
    assert(len(list(bucket.keys()))>0)
    bucketMap.setSizeXYZT (*b['voxel_size'])
    return bucketMap
    
#def getRepresentationFromSSB ( ssb ):
    #glb = []
    #for n in ssb.neighbours():
        #if n.getSyntax() == 'glb':
            #node = oo.Node(n)
            #node.bucket = getBucketFromVertex(n)
            #glb.append ( node )
    #scales = list(set([each.scale for each in glb]))
    #for each in glb:
        #if each.scale == scales[len(scales)/2]:
            #print(each.scale)
            #return getAimsBucketFromBucket (each.bucket)
    #return None
    
        
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
                    node = oo.Node()
                    node.defineFromVertex(n)
                    bucketmap = n['aims_glb']
                    node.bucket = {}
                    node.bucket['voxel_list'] = []
                    assert( len(list(bucketmap[0].keys())) > 0 )
                    for each in bucketmap[0].keys() :
                        node.bucket['voxel_list'].append(each)
                    assert(len(node.bucket['voxel_list']) > 0 )
                    node.bucket['voxel_size'] = [bucketmap.sizeX(), bucketmap.sizeY(), bucketmap.sizeZ(), bucketmap.sizeT()]
                    glb.append ( node )
            scales = list(set([each.scale for each in glb]))
            for each in glb:
                if each.scale == scales[len(scales)/2]:
                    bucketMap = getAimsBucketFromBucket(each.bucket)
                    
            assert ( len(list(bucketMap[0].keys())) > 0 )
            aims.GraphManip.storeAims( graph, v, 'aims_ssb', aims.rc_ptr_BucketMap_VOID(bucketMap) )
            assert ( 'aims_ssb' in v )

def ScaleSpaceBlobsFromOneSubject ( graph, sujet ) :
    blobs = []
    AddSSBBuckets( graph )
    for v in graph.vertices():
        if v.getSyntax() == 'ssb':
            n = {}
            for each in ['t', 'subject', 'tmin', 'tmax']:#, 'index']:
                n[each] = v[each]
            n['bucket'], n['bbmin'], n['bbmax'] = getBucketFromVertex(v)
            print(n['bbmin'], n['bbmax'])
            blobs.append(n)
    return blobs

def GreyLevelBlobsFromOneSubject ( graph, sujet ) :

    blobs = []
    AddSSBBuckets( graph )
    for v in graph.vertices():
        if v.getSyntax() == 'glb':
            n = {}
            for each in ['t', 'subject', 'scale'] : #, 'index']:
                n[each] = v[each]
            n['bucket'], n['bbmin'], n['bbmax'] = getBucketFromVertex(v)
            print(n['t'], n['bbmin'], n['bbmax'])
            blobs.append(n)
    return blobs
    
def ScaleSpaceBlobs ( db, contrast ) :
    ''' Returns a dictionary containing for each subject a list of blobs, a blob
        being a dictionary with some info like t, subject, tmin, tmax, index
        and a bucket '''
    subjects = oo.getSubjects ( db )
    pathes = getfMRIPathes ( db, contrast )
    
    graphs_pathes = {}
    graphs = {}
    blobs = {}
    
    for sujet in subjects:
        graphs[sujet] = aims.read(pathes[sujet]['primal'])
    
    for sujet in subjects :
        blobs[sujet] = ScaleSpaceBlobsFromOneSubject ( graphs[sujet], sujet )
    return blobs

def GreyLevelBlobs ( db, contrast ) :
    ''' Returns a dictionary containing for each subject a list of blobs, a blob
        being a dictionary with some info like t, subject, tmin, tmax, index
        and a bucket '''
    subjects = oo.getSubjects ( db )
    pathes = getfMRIPathes ( db, contrast )
    
    graphs_pathes = {}
    graphs = {}
    blobs = {}
    
    for sujet in subjects:
        graphs[sujet] = aims.read(pathes[sujet]['primal'])
    
    for sujet in subjects :
        blobs[sujet] = GreyLevelBlobsFromOneSubject ( graphs[sujet], sujet )
    return blobs

def resetIndices ( blobs ):
    sujets = list(blobs.keys())
    i=0
    for sujet in sujets:
        for each in blobs[sujet]:
            each['index'] = i
            i = i + 1

def FilterBlobsOnT ( blobs, threshold ) :
    filtered_blobs = {}
    subjects = list(blobs.keys())
    for sujet in subjects:
        filtered_blobs[sujet] = []
        for blob in blobs[sujet] :
            if blob['t'] > threshold :
                filtered_blobs[sujet].append(blob)
    return filtered_blobs


def ComputeCliquesBetweenTwoSubjects ( blobs, sujet1, sujet2 ) :
    cliques = []
    for b1, blob1 in enumerate(blobs[sujet1]):
        for b2, blob2 in enumerate(blobs[sujet2]):
            print(b1, len(blob1['bucket']['voxel_list']), b2, len(blob2['bucket']['voxel_list']))
            if ( not ( blob1['bbmin'][0] > blob2['bbmax'][0] \
                or blob1['bbmin'][1] > blob2['bbmax'][1] \
                or blob1['bbmin'][2] > blob2['bbmax'][2] \
                or blob2['bbmin'][0] > blob1['bbmax'][0] \
                or blob2['bbmin'][1] > blob1['bbmax'][1] \
                or blob2['bbmin'][2] > blob1['bbmax'][2] ) ) :

                c = {}
                c['node1'] = blob1
                c['node2'] = blob2
                c['overlap'] = computeOverlap( blob1, blob2 )

                if ( c['overlap'] > 0.0 ):
                    cliques.append(c)

    return cliques

#def ComputeCliques ( blobs, number_of_proc = 2 ):
    #subjects = blobs.keys()
    #cliques = []
    #jobs = []
    
    #for i, sujet1 in enumerate(subjects):
        #for j, sujet2 in enumerate(subjects):
            #if i < j :
                #print(sujet1,sujet2)
                #jobs.append( (blobs, sujet1, sujet2) )

    #print('JOBS:', len(jobs), 'proc:', number_of_proc)
    #results = MultiProcExecute ( ComputeCliquesBetweenTwoSubjects, jobs, number_of_proc )

    #for each in results :
        #(sujet1, sujet2, resultat) = each
        #print(len(resultat))
        #cliques.extend( convertCliques(resultat) )
                ##cliques_sujets = ComputeCliquesBetweenTwoSubjects ( blobs, sujet1, sujet2 )
                ##print(len(cliques_sujets))
                ##cliques.extend ( cliques_sujets )
    #return cliques

def ComputeCliques ( blobs, number_of_proc = 2 ):
    subjects = list(blobs.keys())
    cliques = []
    jobs = []
    
    for i, sujet1 in enumerate(subjects):
        for j, sujet2 in enumerate(subjects):
            if i < j :
                print(sujet1,sujet2)
                #jobs.append( (blobs, sujet1, sujet2) )
                cliques.extend( convertCliques(ComputeCliquesBetweenTwoSubjects( blobs, sujet1, sujet2 ) ) )

    #print('JOBS:', len(jobs), 'proc:', number_of_proc)
    #results = MultiProcExecute ( ComputeCliquesBetweenTwoSubjects, jobs, number_of_proc )

    #for each in results :
        #(sujet1, sujet2, resultat) = each
        #print(len(resultat))
        #cliques.extend( convertCliques(resultat) )

    return cliques

def BuildAimsGroupGraph ( blobs, cliques ):
    graph = aims.Graph('BlobsArg')
    graph['boundingbox_min'] = [0.0, 0.0, 0.0]
    graph['boundingbox_max'] = [10.0, 10.0, 10.0]
    graph['voxel_size'] = [3.0, 3.0, 3.0]
    graph['filename_base'] = "*"
    sujets = list(blobs.keys())
    
    graph['subjects'] = sujets
    
    glb_indices = {}
    index = 0
    for sujet in sujets : 
        for each in blobs[sujet] :
            v = graph.addVertex('ssb')
            for i in ['t', 'subject', 'index']:
                v[i] = each[i]
            for i in ['tmin', 'tmax']:
                v[i] = each['scale']
            #each['index'] = index
            #index = index + 1
            bucket = getAimsBucketFromBucket(each['bucket'])
            aims.GraphManip.storeAims( graph, v, 'aims_ssb', bucket )
            glb_indices[each['index']] = v
    overlaps = []
    for each in cliques:
        e = graph.addEdge(glb_indices[each.node1['index']], glb_indices[each.node2['index']] , 'b2b')
        e['similarity'] = each.overlap
        overlaps.append(each.overlap)
    import numpy as np
    print(min(overlaps), max(overlaps), np.mean(overlaps))
    return graph
    
messages_defaults = [ ('db path', '/home/go224932/data/structural/database_miccai09'),
                      ('contrast', '24'),
                      ('number of proc', '8'),
                      ('threshold on t', '7.0') ]

if __name__ == '__main__' :
    
    params = InputParameters( sys.argv[1:], messages_defaults )
    if (input('OK? y/n') != 'y'):
        sys.exit(0)

    db = str(params[0])
    contrast = int(params[1])
    number_of_proc = int(params[2])
    threshold = float(params[3])

    sujets = oo.getSubjects( db )

    pathes = getfMRIPathes ( db, contrast )
    #blobs = ScaleSpaceBlobs ( db, contrast )
    blobs = GreyLevelBlobs ( db, contrast )
    sum = 0
    for sujet in sujets :
        sum = sum + len(blobs[sujet])
    print(sum, ' blobs')
    filteredblobs = FilterBlobsOnT ( blobs, threshold )
    resetIndices ( filteredblobs )
    sum = 0
    for sujet in sujets :
        print(sujet, len( filteredblobs[sujet]))
        sum = sum + len( filteredblobs[sujet])
    print(sum, ' filteredblobs')
    

    #distmatrix = getDistanceMatrix ( aims.read(pathes[sujets[0]]['fMRI']) )
    cliques = ComputeCliques ( filteredblobs, number_of_proc )
    graph = BuildAimsGroupGraph ( filteredblobs, cliques )
    groupgraphpath = os.path.join( os.path.split( db )[0], 'groupgraphs', 'groupgraph_%s.arg' %(os.path.split( db )[1]) )
    print(groupgraphpath)
    aims.write ( graph, groupgraphpath )

        
    
    