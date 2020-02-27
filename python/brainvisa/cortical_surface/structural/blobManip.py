from __future__ import print_function

from __future__ import absolute_import
import os
from soma import aims
from brainvisa.cortical_surface.shell import db


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
        self.index = None
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
        if 'nodes' in v:
            self.nodes = list(v['nodes'])
        self.index = int( v['index'] )
        self.subject = str(v['subject'])
        self.t = float(v['t'])
        self.scale = float(v['scale'])
    
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

def Distrib ( db_path, hemis_side = 'L' , data_type = 'curv' ) :
    """ Computes a distribution based on local maxima extracted from specific data 
        return distrib """
    nodes = MaximaList ( db_path, hemis_side, data_type )
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
                    print(v['t'], v['subject'])

def JointDistrib ( db_path, hemis_side = 'L' ) :
    """ Get a joint distribution between curv and depth data. Are the deep max 
      the same as the curved ones ? 
      return curvlist, depthlist, nodes, sujets"""
    sujets = db.getSubjects ( db_path )
    curvpathes = db.getPathes ( db_path, hemis_side, 'curv' )
    depthpathes = db.getPathes ( db_path, hemis_side, 'depth' )

    curvgraphspathes = [curvpathes[sujet]['primal'] for sujet in sujets]
    depthgraphspathes = [depthpathes[sujet]['primal'] for sujet in sujets]
    
    curvlist = []
    depthlist= []
    nodes = []
    sujets = []
    
    for i, curvgraphpath in enumerate ( curvgraphspathes ) :
        curvgraph = aims.read ( curvgraphpath )
        depthgraph = aims.read ( depthgraphspathes[i] )
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
            if n in depthnodes:
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
        
def LocalMaximaCountPerSubject ( db_path, data_type = 'curv' ) :
    ''' How many local maxima can we find per subject ? This allows to measure the
    relation between this number and age. Data_type can be curv, depth or joint'''
    import numpy as np
    if ( data_type in ['curv', 'depth'] ):
        distribL, sujetsL =  Distrib ( db_path, 'L', data_type )
        distribR, sujetsR =  Distrib ( db_path, 'R', data_type )
    elif ( data_type in ['joint'] ) :
        curvlistL, depthlistL, nodesL, sujetsL =  JointDistrib( db_path, 'L')
        curvlistR, depthlistR, nodesR, sujetsR =  JointDistrib( db_path, 'R')
    count = {}
    sujets = db.getSubjects( db_path )
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
        x.append ( count[each] )
        y.append ( ages[each] )
    np.corrcoef (x, y )
    return x, y

def SubjectNodes ( nodes, subjects ):
    """ returns a list of nodes belonging to a given list of subjects """
    nodes_subject = []
    for each in nodes:
        if each.subject in subjects:
            nodes_subject.append(each)

    return nodes_subject

def SimpleStatsFromCountTable ( count, subjects = None ) :
    import numpy as np
    ct_ave = []
    ct_std = []
    if subjects is None :
        sujets = list(count.keys())
    else :
        sujets = subjects
    distances = list(np.sort(count[sujets[0]].keys()))
    for dist in distances :
        ct_ave.append ( np.mean( [len(count[sujet][dist]) for sujet in sujets]) )
        ct_std.append ( np.std( [len(count[sujet][dist]) for sujet in sujets]) )
        
    return ct_ave, ct_std, distances

def ScaleSpaceBlobs ( db_path, contrast ) :
    ''' Returns a dictionary containing for each subject a list of blobs, a blob
        being a dictionary with some info like t, subject, tmin, tmax, index
        and a bucket '''
    subjects = db.getSubjects ( db_path )
    pathes = db.getfMRIPathes ( db_path, contrast )
    
    graphs_pathes = {}
    graphs = {}
    blobs = {}
    
    for sujet in subjects:
        graphs[sujet] = aims.read(pathes[sujet]['primal'])
    
    for sujet in subjects :
        blobs[sujet] = ScaleSpaceBlobsFromOneSubject ( graphs[sujet], sujet )
    return blobs