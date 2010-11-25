import os
from soma import aims

class Node:
    ''' Node class (with node, subject and t fields)
    remember : t is supposed to be the value at the max node on the original data'''
    def MaximumNode ( self, data ) :        
        maxi = self.nodes[0];
        for i in self.nodes :
            if ( data[0][i] > data[0][maxi] ):
                maxi = i
        return maxi
        
    def __repr__ ( self ) :
        s = ' nodes: %s,\n maxnode: %s,\n subject: %s,\n t: %s,\n scale: %s' \
            %(str(self.nodes), self.maxnode, self.subject, float(self.t), self.scale)
        return s
    
    def __init__ ( self ):
        self.nodes = None
        self.maxnode  = None
        self.subject = None
        self.t = None
        self.scale = None
        
    def defineFromArgs ( self, nodes, t, subject, scale ):
        self.nodes = list(nodes)
        self.subject = str(subject)
        self.t = float(t)
        self.scale = float(scale)
        
    def defineFromVertex ( self, v ) :
        if v.has_key('nodes'):
            self.nodes = list(v['nodes'])
        self.subject = str(v['subject'])
        self.t = float(v['t'])
        self.scale = float(v['scale'])

def SubjectsSubset(data, nb):  
    from random import randrange
    new_data = [d for d in data]
    assert(nb<=len(data))
    comp = []
    for n in xrange(nb):
        pos = randrange( len(new_data) )
        elem = new_data[pos]
        new_data[pos] = new_data[-1]
        del new_data[-1]
        comp.append( elem )
    return comp

def getSubjects ( db ) :
    """Returns the list of subjects (dirnames) contained in the data base directory
         (ignores any existing 'tmp')"""
    sujets = []
    for direc in os.listdir ( db ):
        if os.path.isdir( os.path.join(db, direc) ) and direc != 'tmp':
            sujets.append ( direc )
    return sujets

def getPathes ( db, hemis_side = 'L' , data_type = 'curv', subjects = None ) :
    """Builds and returns the pathes according to the data base path, hemisphere side
         and data type (curvature or depth)"""
    pathes = {}
    if subjects is None :
        subjects = getSubjects ( db )
    
    for sujet in subjects :
        pathes[sujet] = {}
        
        pathes[sujet]['mesh'] = os.path.join ( db, sujet, 't1mri/t1/default_analysis/segmentation/mesh/%s_%swhite.mesh' % ( sujet, hemis_side ))
        pathes[sujet]['tex'] = os.path.join ( db, sujet, 't1mri/t1/default_analysis/segmentation/%s_%swhite_%s.tex' % ( sujet, hemis_side, data_type ) )
        pathes[sujet]['primal'] = os.path.join( db, sujet, 'surface/%s_%swhite_primal_%s.arg' % ( sujet, hemis_side, data_type ) )
        
        pathes[sujet]['blobs'] = os.path.join( db, sujet, 'surface/%s_%swhite_%s_blobs.tex' % ( sujet, hemis_side, data_type ) )
        pathes[sujet]['scalespace'] = os.path.join ( db, sujet, 'surface/%s_%swhite_%s_ss.tex' % ( sujet, hemis_side, data_type ) )
        pathes[sujet]['lat'] = os.path.join ( db, sujet, 'surface/%s_%s_lat.tex' % ( sujet, hemis_side ) )
        pathes[sujet]['lon'] = os.path.join ( db, sujet, 'surface/%s_%s_lon.tex' % ( sujet, hemis_side ) )
        pathes[sujet]['gyri'] = os.path.join ( db, sujet, 'surface/%s_%s_gyri.tex' % ( sujet, hemis_side ) )
        
    return pathes
    
def MaximaList ( db, hemis_side = 'L' , data_type = 'curv' ) :
    """ Computes a list of maxima based on a type of data (curv/depth)
        The list is a Node list """
    sujets = getSubjects( db )
    pathes = getPathes( db, hemis_side, data_type )
    graphspathes = [pathes[sujet]['primal'] for sujet in sujets]
    nodes = []
    
    textures = {}
    for sujet in sujets:
        textures[sujet] = aims.read(pathes[sujet]['tex'])
        
    for graphpath in graphspathes :
        graph = aims.read(graphpath)
        for v in graph.vertices():
            if v.getSyntax() == 'glb':
                node = Node()
                node.defineFromVertex(v)
                node.maxnode = node.MaximumNode(textures[ v['subject'] ] )
                nodes.append( node )
    return nodes

def Distrib ( db, hemis_side = 'L' , data_type = 'curv' ) :
    """ Computes a distribution based on local maxima extracted from specific data 
        return distrib """
    nodes = MaximaList ( db, hemis_side, data_type )
    distrib = [ n.t for n in nodes ]
    subjects = [ n.subject for n in nodes ]
    return distrib, subjects
            
def spotOutliers ( db, hemis_side = 'L' , data_type = 'curv' ) :
    """ 14/10/10 : just to spot the outlying maxima. Indeed the depth map was 
      wrong with sujet01"""
    import matplotlib.pyplot as plt
    sujets = getSubjects(db)
    pathes = getPathes(db, hemis_side, data_type)
    graphspathes = [pathes[sujet]['primal'] for sujet in sujets]
    distrib = []
    for graphpath in graphspathes :
        graph = aims.read(graphpath)
        for v in graph.vertices():
            if v.getSyntax() == 'glb':
                if float(v['t']) > 50.0:
                    print v['t'], v['subject']
                    
def JointDistrib( db, hemis_side = 'L' ) :
    """ Get a joint distribution between curv and depth data. Are the deep max 
      the same as the curved ones ? 
      return curvlist, depthlist, nodes, sujets"""
    sujets = getSubjects(db)
    curvpathes = getPathes(db, hemis_side, 'curv')
    depthpathes = getPathes(db, hemis_side, 'depth')

    curvgraphspathes = [curvpathes[sujet]['primal'] for sujet in sujets]
    depthgraphspathes = [depthpathes[sujet]['primal'] for sujet in sujets]
    curvlist = []
    depthlist= []
    nodes = []
    sujets = []
    for i, curvgraphpath in enumerate(curvgraphspathes) :
        curvgraph = aims.read(curvgraphpath)
        depthgraph = aims.read(depthgraphspathes[i])
        curvnodes = {}
        depthnodes = {}
        for v in curvgraph.vertices():
            if v.getSyntax() == 'glb':
                for n in v['nodes']:
                    curvnodes[int(n)] = v['t'], v['subject']
        for v in depthgraph.vertices():
            if v.getSyntax() == 'glb':
                for n in v['nodes']:
                    depthnodes[int(n)] = v['t'], v['subject']

        for n, curv in curvnodes.items():
            if depthnodes.has_key(n):
                sujets.append(str(curv[1]))
                nodes.append(n)
                curvlist.append(curv[0])
                depthlist.append(depthnodes[n][0])
    return curvlist, depthlist, nodes, sujets
            
def getAges ( sujets ):
    ages = { 'Vil1' : 26.6,
             'Hus1' : 28.6,
             'Yel1' : 24.4,
             'Eck1' : 35.0,
             'Maz1' : 29.9,
             'Cad1' : 30.4,
             'Kol1' : 31.0,
             'Ozc1' : 27.9,
             'Bar1' : 30.6,
             'Sav1' : 33.6,
             'Del1' : 35.7,
             'DaE1' : 29.6,
             'Mun1' : 33.3,
             'Rem1' : 30.0,
             'Ong1' : 30.3,
             'Fou1' : 31.1,
             'Joh1' : 32.6,
             'Bal1' : 31.1,
             'Bor1' : 33.7,
             'Sha1' : 34.4,
             'Pou1' : 34.4,
             'Roi1' : 30.7,
             'Ben1' : 28.1,
             'Han1' : 32.1,
             'Bur1' : 32.0 }
    res = {}
    for i in sujets :
        if i in ages.keys():
            res[i] = ages[i]
    return res
        
def LocalMaximaCountPerSubject( db, data_type = 'curv' ):
    ''' How many local maxima can we find per subject ? This allows to measure the
    relation between this number and age. Data_type can be curv, depth or joint'''
    import numpy as np
    if ( data_type in ['curv', 'depth'] ):
        distribL, sujetsL =  getDistrib( db, 'L', data_type )
        distribR, sujetsR =  getDistrib( db, 'R', data_type )
    elif ( data_type in ['joint'] ) :
        curvlistL, depthlistL, nodesL, sujetsL =  JointDistrib( db, 'L')
        curvlistR, depthlistR, nodesR, sujetsR =  JointDistrib( db, 'R')
    count = {}
    sujets = getSubjects( db )
    for sujet in sujets:
      count[sujet] = 0
      
    for sujet in sujetsL:
        count[sujet] = count[sujet] + 1
    for sujet in sujetsR:
        count[sujet] = count[sujet] + 1
    
    ages = getAges( sujets )
    x = []
    y = []
    for each in sujets:
        x.append( count[each] )
        y.append( ages[each] )
    np.corrcoef(x,y)
    return x, y

def Graph ( db, nodes, hemis_side = 'L', data_type = 'curv' ):
    ''' Creates an Aims Graph given a list of nodes (vertex indices + subjects)
      and a file path. In this function a node is a Node'''
    from soma import aims
    import numpy as np 
    sujets = getSubjects( db )
    pathes = getPathes( db, hemis_side, data_type )
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
        v = graph.addVertex('glb')
        v['subject'] = n.subject        
        v['nodes'] = n.nodes
        v['scale'] = n.scale
        v['name'] = str(sujets.index(n.subject) + 1)        
        v['t'] = n.t
        v['x'] = list(np.ones(len(n.nodes)))
        v['y'] = list(np.ones(len(n.nodes)))
        #assert( len(n.nodes) == 1 )
        c = meshes[ n.subject ].vertex()[ n.maxnode ]
        mesh = aims.SurfaceGenerator.sphere(c, 0.3, 10, True)
        aims.GraphManip.storeAims( graph, v, 'glb', mesh )
   
    #aims.write( graph, graph_path )
    return graph
    
    
    
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

def StatsClusters( clusters, db, hemis_side = 'L', data_type = 'curv' ):
    """ returns a count table (how many glbs per ssb?) and a distances table (how
           compact are the clusters?"""
    from soma import aims
    count = [ len(cluster['glbs'] ) for cluster in clusters]
    sujets = [ cluster['ssb'].subject for cluster in clusters ]
    set_sujets = list( set(sujets) )
    assert( len( set_sujets ) == 1 )
    sujet = set_sujets[0]
    pathes = getPathes ( db, hemis_side, data_type )
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
                
def SubjectNodes( nodes, subjects ):
    """ returns a list of nodes belonging to a given list of subjects """
    nodes_subject = []
    for each in nodes:
        #print each
        if each.subject in subjects:
            nodes_subject.append(each)
            
    return nodes_subject
    
def ClustersOneSubjectAtVariousDistances ( db, raw_max_subject_nodes, hemis_side = 'L', data_type = 'curv' ) :
    '''The nodes taken in input are the result of MaximaList() that is a list
        of raw local maxima extracted from specific data
        return clusters, count, dist '''
    import numpy as np
    clusters = {}
    stats = {}
    for distance in list ( np.arange(0.5, 25.0, 0.5) ):
        clusters[distance] = Clusters ( db, raw_max_subject_nodes, distance, hemis_side, data_type )
        stats[distance] = StatsClusters ( clusters[distance], db, hemis_side, data_type )
    count = {}
    dist = {}
    for distance in list( np.arange(0.5, 25.0, 0.5) ):
        count[distance] = stats[distance][0]
        dist[distance] = stats[distance][1]
    return clusters, count, dist
    
def ClustersOneSubjectAtVariousDistances ( db, raw_max_subject_nodes, hemis_side = 'L', data_type = 'curv' ) :
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
    
    graph = Graph ( db, raw_max_subject_nodes, hemis_side, data_type )
    tex = aims.read( graph['texture'] )
        
    u = str(uuid.uuid1())[:6]
    clustergraphinpath = os.path.join ( tmpdir, 'clustergraph_in_%s_%s.arg'%(str(set_sujets[0]), u) )
    clustergraphoutpath = os.path.join ( tmpdir, 'clustergraph_out_%s_%s.arg'%(str(set_sujets[0]), u) )
    aims.write ( graph, clustergraphinpath )
    
    textoutputpath = '/tmp/blobsCountTable_' + str(set_sujets[0]) + '_' + str(u) +'.py'
    
    cluster_command = 'surfCreateClustersFromGLB -i  '+ str(clustergraphinpath) +'  -o '+ str(clustergraphoutpath) + ' --textOutput ' + str(textoutputpath)
    os.system ( cluster_command )
    print cluster_command
    charac_clusters = {}
    count_glb = {}
    execfile ( str(textoutputpath), locals(), globals() )
    
    return charac_clusters, count_glb


def SimpleStatsFromCountTable ( count, subjects = None ) :
    import numpy as np
    ct_ave = []
    ct_std = []
    if subjects is None :
        sujets = count.keys()
    else :
        sujets = subjects
    distances = list(np.sort(count[sujets[0]].keys()))
    for dist in distances :
        ct_ave.append ( np.mean( [len(count[sujet][dist]) for sujet in sujets]) )
        ct_std.append ( np.std( [len(count[sujet][dist]) for sujet in sujets]) )
        
    return ct_ave, ct_std, distances
    

def closestNode ( mesh, point ):
  import math
  distance = 0.0
  mini = 0
  distmini = 1000.0
  for v in xrange( len( mesh.vertex()) ):
      vert = mesh.vertex()[v]
      d = [ math.pow(each, 2) for each in [vert[i] - point[i] for i in xrange(3)] ]
      distance = math.sqrt( d[0] + d[1] + d[2] )

      if ( distance < distmini ):
          distmini = distance
          mini = v
  return mini
      
 
  
def ClustersOnSimulatedRegions ( db, radius = 25.0, hemis_side = 'L', data_type = 'curv' ):
    
    clusters = {}
    count = {}
    dist = {}
    import os, shutil, random
    from soma import aims
    pathes = getPathes( db, hemis_side, data_type )
    sujets = getSubjects ( db )
    mesh = aims.read( pathes[sujets[0]]['mesh'] )
    #mesh = aims.read( pathes['Bor1']['mesh'] )
    node = random.randint( 0, len(mesh.vertex()) )
    #node = 303 #to aim at the central sulcus roughly
    point = aims.Point3df(mesh.vertex()[node])

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
     
        

    