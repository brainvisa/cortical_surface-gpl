#!/usr/bin/env python

import os,sys
from soma import aims
import numpy as np
import math
import brainvisa.cortical_surface.structural.blobManip as oo
from brainvisa.cortical_surface.multiprocessing.mproc import MultiProcExecute
from brainvisa.cortical_surface.shell.inputParameters import *



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

def resetIndices ( blobs ):
    sujets = blobs.keys()
    i = 0
    for sujet in sujets:
        for each in blobs[sujet]:
            each['index'] = i
            i = i + 1

def FilterSSBOnT ( blobs, threshold ) :
    filtered_blobs = {}
    subjects = blobs.keys()
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
            print b1, len(blob1['bucket']['voxel_list']), b2, len(blob2['bucket']['voxel_list'])
            if ( not ( blob1['bbmin'][0] > blob2['bbmax'][0] \
                or blob1['bbmin'][1] > blob2['bbmax'][1] \
                or blob1['bbmin'][2] > blob2['bbmax'][2] \
                or blob2['bbmin'][0] > blob1['bbmax'][0] \
                or blob2['bbmin'][1] > blob1['bbmax'][1] \
                or blob2['bbmin'][2] > blob1['bbmax'][2] ) ) :

                c = {}
                c['node1'] = blob1
                c['node2'] = blob2
                c['overlap'] = oo.computeOverlap( blob1, blob2 )

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
                #print sujet1,sujet2
                #jobs.append( (blobs, sujet1, sujet2) )

    #print 'JOBS:', len(jobs), 'proc:', number_of_proc
    #results = MultiProcExecute ( ComputeCliquesBetweenTwoSubjects, jobs, number_of_proc )

    #for each in results :
        #(sujet1, sujet2, resultat) = each
        #print len(resultat)
        #cliques.extend( convertCliques(resultat) )
                ##cliques_sujets = ComputeCliquesBetweenTwoSubjects ( blobs, sujet1, sujet2 )
                ##print len(cliques_sujets)
                ##cliques.extend ( cliques_sujets )
    #return cliques

def ComputeCliques ( blobs, number_of_proc = 2 ):
    subjects = blobs.keys()
    cliques = []
    jobs = []
    
    for i, sujet1 in enumerate(subjects):
        for j, sujet2 in enumerate(subjects):
            if i < j :
                print sujet1,sujet2
                #jobs.append( (blobs, sujet1, sujet2) )
                cliques.extend( oo.convertCliques(ComputeCliquesBetweenTwoSubjects( blobs, sujet1, sujet2 ) ) )

    #print 'JOBS:', len(jobs), 'proc:', number_of_proc
    #results = MultiProcExecute ( ComputeCliquesBetweenTwoSubjects, jobs, number_of_proc )

    #for each in results :
        #(sujet1, sujet2, resultat) = each
        #print len(resultat)
        #cliques.extend( convertCliques(resultat) )

    return cliques

def BuildAimsGroupGraph ( blobs, cliques ):
    graph = aims.Graph('BlobsArg')
    graph['boundingbox_min'] = [0.0, 0.0, 0.0]
    graph['boundingbox_max'] = [10.0, 10.0, 10.0]
    graph['voxel_size'] = [3.0, 3.0, 3.0]
    graph['filename_base'] = "*"
    sujets = blobs.keys()
    
    graph['subjects'] = sujets
    
    glb_indices = {}
    index = 0
    for sujet in sujets : 
        for each in blobs[sujet] :
            v = graph.addVertex('ssb')
            for i in ['t', 'subject', 'tmin', 'tmax', 'index']:
                v[i] = each[i]
            #each['index'] = index
            #index = index + 1
            bucket = oo.getAimsBucketFromBucket(each['bucket'])
            aims.GraphManip.storeAims( graph, v, 'aims_ssb', bucket )
            glb_indices[each['index']] = v
    overlaps = []
    for each in cliques:
        e = graph.addEdge(glb_indices[each.node1['index']], glb_indices[each.node2['index']] , 'b2b')
        e['similarity'] = each.overlap
        overlaps.append(each.overlap)
    import numpy as np
    print min(overlaps), max(overlaps), np.mean(overlaps)
    return graph
    
messages_defaults = [ ('db path', '/home/go224932/data/structural/database_miccai09'),
                      ('contrast', '24'),
                      ('number of proc', '8'),
                      ('threshold on t', '7.0') ]

if __name__ == '__main__' :
    
    params = InputParameters( sys.argv[1:], messages_defaults )
    if (raw_input('OK? y/n') != 'y'):
        sys.exit(0)

    db = str(params[0])
    contrast = int(params[1])
    number_of_proc = int(params[2])
    threshold = float(params[3])

    sujets = oo.getSubjects( db )

    pathes = getfMRIPathes ( db, contrast )
    blobs = oo.ScaleSpaceBlobs ( db, contrast )
    sum = 0
    for sujet in sujets :
        sum = sum + len(blobs[sujet])
    print sum, ' blobs'
    filteredblobs = FilterSSBOnT ( blobs, threshold )
    resetIndices ( filteredblobs )
    sum = 0
    for sujet in sujets :
        print sujet, len( filteredblobs[sujet])
        sum = sum + len( filteredblobs[sujet])
    print sum, ' filteredblobs'
    

    #distmatrix = getDistanceMatrix ( aims.read(pathes[sujets[0]]['fMRI']) )
    cliques = ComputeCliques ( filteredblobs, number_of_proc )
    graph = BuildAimsGroupGraph ( filteredblobs, cliques )
    groupgraphpath = os.path.join( os.path.split( db )[0], 'groupgraphs', 'groupgraph_%s.arg' %(os.path.split( db )[1]) )
    print groupgraphpath
    aims.write ( graph, groupgraphpath )

        
    
    