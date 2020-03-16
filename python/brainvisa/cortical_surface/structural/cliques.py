from __future__ import print_function


def computeOverlap ( blob1, blob2 ) :
    bucket1 = blob1['bucket']['voxel_list']
    bucket2 = blob2['bucket']['voxel_list']
    intersection = []
    div = len(bucket1) + len(bucket2)
    for each1 in bucket1:
        for each2 in bucket2:
            if (each1[0] == each2[0] and each1[1] == each2[1] and each1[2] == each2[2]):
                intersection.append(each1)
    return float ( 2.0 * float(len(intersection)) / float(div) )

def convertCliques( rawcliques ) :
    cliques = []
    for each in rawcliques:
        c = Clique()
        c.define(each)
        cliques.append(c)
    return cliques

class Clique(object):
    def __init__( self ):
        pass
    def define ( self, c ):
        self.node1 = c['node1']
        self.node2 = c['node2']
        self.overlap = c['overlap']

def ComputeCliquesBetweenTwoSubjects ( blobs, sujet1, sujet2 ) :
    cliques = []
    for b1, blob1 in enumerate ( blobs[sujet1] ) :
        for b2, blob2 in enumerate ( blobs[sujet2] ):
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

def ComputeCliques ( blobs, number_of_proc = 2 ):
    subjects = list(blobs.keys())
    cliques = []
    jobs = []
    
    for i, sujet1 in enumerate ( subjects ):
        for j, sujet2 in enumerate ( subjects ):
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
