
from brainvisa.cortical_surface.structural.graphManip import Graph
from brainvisa.cortical_surface.shell import db
from brainvisa.cortical_surface.structural.blobManip import Node
import os


def Clusters ( db, raw_max_nodes, clustering_distance, hemis_side = 'L', data_type = 'curv' ) :
    """ Calls the surfCreateClustersFromGLB command, writes the output in a temporary graph_path
          and returns the produced blobs and a list of clusters (list int) 
          Works with Nodes From One Subject Only !"""
    from soma import aims
    import shutil, os 
    sujets = [node.subject for node in raw_max_nodes]
    set_sujets = list( set(sujets) )
    assert ( len(set_sujets) == 1 )
    
    tmpdir = os.path.join(db, 'tmp')
    if( os.path.exists( tmpdir ) ):
        shutil.rmtree ( tmpdir )
    os.mkdir( tmpdir )
    
    graph = Graph ( db, raw_max_nodes, hemis_side, data_type )
    tex = aims.read( graph['texture'] )
    clustergraphinpath = os.path.join ( tmpdir, 'clustergraph_in.arg' )
    clustergraphoutpath = os.path.join ( tmpdir, 'clustergraph_out.arg' )
    aims.write ( graph, clustergraphinpath )
    
    textoutputpath = '/tmp/blobsCountTable_' + str(set_sujets[0]) + '.py'
    cluster_command = 'surfCreateClustersFromGLB -i  '+ str(clustergraphinpath) +'  -o '+ str(clustergraphoutpath) +' --dist '+ str(clustering_distance) + ' --textOutput ' + str(textouputpath)
    os.system( cluster_command )
    print cluster_command
    
    graph_out = aims.read ( clustergraphoutpath )
    nodes_out = []
    
    for v in graph_out.vertices():
        if v.getSyntax() == 'ssb':
            node = Node()
            #print v
            node.defineFromArgs( nodes = v['nodes'], t = v['t'], subject = v['subject'], scale = max(v['scales']) - min(v['scales']) )
            node.maxnode = node.MaximumNode( tex )
            neigh = []
            for n in v.neighbours():
                if n.getSyntax() == 'glb':
                    glb_node = Node ()
                    #print n
                    glb_node.defineFromArgs ( nodes = n['nodes'], t = n['t'], subject = n['subject'], scale = n['scale'] )
                    glb_node.maxnode = glb_node.MaximumNode ( tex )
                    neigh.append( glb_node  )
            nodes_out.append( {'ssb' : node, 'glbs' : neigh } )

    return nodes_out

def StatsClusters( clusters, db_path, hemis_side = 'L', data_type = 'curv' ):
    """ returns a count table (how many glbs per ssb?) and a distances table (how
           compact are the clusters?"""
    from soma import aims
    count = [ len(cluster['glbs'] ) for cluster in clusters]
    sujets = [ cluster['ssb'].subject for cluster in clusters ]
    set_sujets = list( set(sujets) )
    assert( len( set_sujets ) == 1 )
    sujet = set_sujets[0]
    pathes = db.getPathes ( db_path, hemis_side, data_type )
    mesh = aims.read ( pathes[sujet]['mesh'] )
    distances = []
    for cluster in clusters :
        dist_moy = 0.0
        if ( len( cluster['glbs'] ) > 1 ) :
            n = 0.0
            for i in range( 0, len( cluster['glbs'] ) - 1 ):
                for j in range( i + 1, len( cluster['glbs'] ) ):
                    a = mesh.vertex()[ cluster['glbs'][i].maxnode ]
                    b = mesh.vertex()[ cluster['glbs'][j].maxnode ]
                    c = a - b
                    dist_moy = dist_moy + c.norm()
                    n = n + 1.0
            dist_moy = dist_moy / float(n)
        distances.append( dist_moy )
        
    return count, distances

def ClustersOneSubjectAtVariousDistances ( db_path, raw_max_subject_nodes, hemis_side = 'L', data_type = 'curv' ) :
    '''The nodes taken in input are the result of MaximaList() that is a list
        of raw local maxima extracted from specific data
        return clusters, count, dist '''
    import numpy as np
    clusters = {}
    stats = {}
    for distance in list ( np.arange(0.5, 25.0, 0.5) ):
        clusters[distance] = Clusters ( db_path, raw_max_subject_nodes, distance, hemis_side, data_type )
        stats[distance] = StatsClusters ( clusters[distance], db_path, hemis_side, data_type )
    count = {}
    dist = {}
    for distance in list( np.arange(0.5, 25.0, 0.5) ):
        count[distance] = stats[distance][0]
        dist[distance] = stats[distance][1]
    return clusters, count, dist
    
def ClustersOneSubjectAtVariousDistances ( db_path, raw_max_subject_nodes, hemis_side = 'L', data_type = 'curv' ) :
    """ Calls the surfCreateClustersFromGLB command, writes the output in a temporary graph_path
          and returns the produced blobs and a list of clusters (list int) 
          Works with Nodes From One Subject Only !"""
    from soma import aims
    import shutil, os, uuid
    sujets = [node.subject for node in raw_max_subject_nodes]
    set_sujets = list( set(sujets) )
    assert ( len(set_sujets) == 1 )
    
    tmpdir = os.path.join( db, 'tmp' )
    if( not os.path.exists( tmpdir ) ):
        os.mkdir( tmpdir )
    
    graph = Graph ( db_path, raw_max_subject_nodes, hemis_side, data_type )
    tex = aims.read( graph['texture'] )
        
    u = str(uuid.uuid1())[:6]
    clustergraphinpath = os.path.join ( tmpdir, 'clustergraph_in_%s_%s.arg'%(str(set_sujets[0]), u ) )
    clustergraphoutpath = os.path.join ( tmpdir, 'clustergraph_out_%s_%s.arg'%(str(set_sujets[0]), u ) )
    aims.write ( graph, clustergraphinpath )
    
    textoutputpath = '/tmp/blobsCountTable_' + str(set_sujets[0]) + '_' + str(u) +'.py'
    
    cluster_command = 'surfCreateClustersFromGLB -i  '+ str(clustergraphinpath) +'  -o '+ str(clustergraphoutpath) + ' --textOutput ' + str(textoutputpath)
    os.system ( cluster_command )
    print cluster_command
    charac_clusters = {}
    count_glb = {}
    execfile ( str(textoutputpath), locals(), globals() )

    return charac_clusters, count_glb


def ClustersOnSimulatedRegions ( db_path, radius = 25.0, hemis_side = 'L', data_type = 'curv' ):

    def closestNode ( mesh, point ):
        import math
        distance = 0.0
        mini = 0
        distmini = 1000.0
        for v in xrange ( len( mesh.vertex()) ):
            vert = mesh.vertex()[v]
            d = [ math.pow(each, 2) for each in [vert[i] - point[i] for i in xrange(3)] ]
            distance = math.sqrt( d[0] + d[1] + d[2] )

            if ( distance < distmini ):
                distmini = distance
                mini = v
        return mini

    clusters = {}
    count = {}
    dist = {}
    import os, shutil, random
    from soma import aims
    pathes = db.getPathes( db_path, hemis_side, data_type )
    sujets = db.getSubjects ( db_path )
    mesh = aims.read ( pathes[sujets[0]]['mesh'] )
    #mesh = aims.read( pathes['Bor1']['mesh'] )
    node = random.randint ( 0, len(mesh.vertex()) )
    #node = 303 #to aim at the central sulcus roughly
    point = aims.Point3df ( mesh.vertex()[node] )

    for sujet in sujets:
        mesh = aims.read( pathes[sujet]['mesh'] )
        node = closestNode( mesh, point )
        #node = random.randint( 0, len(mesh.vertex()) ) # to renew the choice of noce from one subjet to another : the regions then become completely unrelated

        #radius = 25.0
        random_command = 'surfRandomROI -m ' + str(pathes[sujet]['mesh'])  +' -o '+ str(pathes[sujet]['lat']) + ' --node ' + str(node) +' --radius '+ str(radius)
        os.system ( random_command )
        print random_command
        shutil.copy(pathes[sujet]['lat'], pathes[sujet]['lon'])
        
        primal_command = 'surfMesh2Graph -m  '+ str(pathes[sujet]['mesh']) +'  -t '+ str(pathes[sujet]['tex']) +' -g '+ str(pathes[sujet]['primal']) +' -s '+str(sujet) + \
            ' --ss '+ str(pathes[sujet]['scalespace']) +' --recover False --blobs '+ str(pathes[sujet]['blobs']) 
        if ( os.path.exists(pathes[sujet]['lat']) and os.path.exists(pathes[sujet]['lon']) ):
            primal_command = primal_command + ' --lat '+ str(pathes[sujet]['lat']) +' --lon '+ str(pathes[sujet]['lon'])
        os.system( primal_command )
        
    all_nodes = MaximaList ( db, hemis_side, data_type )
    for sujet in sujets:
        nodes = SubjectNodes ( all_nodes, [sujet])
        #clusters[sujet], count[sujet], dist[sujet] = ClustersOneSubjectAtVariousDistances ( db, nodes, hemis_side, data_type )
        print sujet
        clus, coun = ClustersOneSubjectAtVariousDistances2 ( db, nodes, hemis_side, data_type )        
        clusters[sujet] = clus[sujet]
        count[sujet] = coun[sujet]
    return clusters, count