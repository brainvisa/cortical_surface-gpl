#  This software and supporting documentation are distributed by
#      Institut Federatif de Recherche 49
#      CEA/NeuroSpin, Batiment 145,
#      91191 Gif-sur-Yvette cedex
#      France
#
# This software is governed by the CeCILL license version 2 under
# French law and abiding by the rules of distribution of free software.
# You can  use, modify and/or redistribute the software under the 
# terms of the CeCILL license version 2 as circulated by CEA, CNRS
# and INRIA at the following URL "http://www.cecill.info". 
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license version 2 and that you accept its terms.
from neuroProcesses import *
import shfjGlobals     

name = '1 - Create Surface-Based Statistical Parametric Maps'
userLevel = 2

signature = Signature(  'intmesh', ListOf(ReadDiskItem( 'Hemisphere White Mesh', 'MESH mesh' )), 
        'time_texture', ListOf(ReadDiskItem('Functional Time Texture', 'Texture')),
        'protocol_text', ReadDiskItem( 'Text File', 'Text File' ),
        'contrast', String(),
        'contrast_name', String(),
        'beta', ListOf(WriteDiskItem('Surface-Based Beta Map', 'Texture')),
        'spmt_maps', ListOf(WriteDiskItem( 'Surface-Based SPMt Map', 'Texture')),
  )

def getContrastName(self, data):
  if (self.contrast_name is not None and len(self.intmesh) == len(self.time_texture) and len(self.intmesh)>0):
      result = []
      for mesh in self.intmesh:
         attributes = mesh.hierarchyAttributes()
         attributes[ 'contrast' ] = str(self.contrast_name)
         print attributes[ 'contrast' ]
         res = self.signature[ 'spmt_maps' ].findValue( attributes )
         print res[0]
         result.append( res[0] )
      return result
  return None
      
def getBetaName(self, data):
  if (self.contrast_name is not None and len(self.intmesh) == len(self.time_texture) and len(self.intmesh)>0):
      result = []
      for mesh in self.intmesh:
         attributes = mesh.hierarchyAttributes()
         attributes[ 'contrast' ] = str(self.contrast_name)
         print attributes[ 'contrast' ]
         res = self.signature[ 'beta' ].findValue( attributes )
         print res[0]
         result.append( res[0] )
      return result
  return None
  

def initialization( self ):

    self.setOptional('beta')
    self.linkParameters('time_texture','intmesh')
    self.addLink('spmt_maps','intmesh',self.getContrastName)
    self.addLink('spmt_maps','contrast_name',self.getContrastName)
    self.addLink('beta','intmesh',self.getBetaName)
    self.addLink('beta','contrast_name',self.getBetaName)
    
    
#def processtag(reader):
   #testa = str(reader.Value())
   #j=0
   #k=1
   #for test in string.split(testa,'\n'):
      #if test == "" or test == "None":
         #continue
      #for i in test:
         #if i not in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '0']:
            #j=j+1
            #continue
         #if i in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '0']:
            #k=k+1
            #continue
         #break
      #testi = int((test[j:k+j]))
      ##print reader.Name()
      #return testi

#def streamFile(filename):
    #import libxml2
    #try:
        #reader = libxml2.newTextReaderFilename(filename)
    #except Exception, e:
        #print e
        #print "unable to open %s" % (filename)
        #return
    #TR=24
    #t_time = []
    #t_type = []
    #hrf = []
    #ret = reader.Read()
    #while ret == 1:
        ##processNode(reader)
        #if reader.Name() == "time" and int(reader.NodeType()) == 1:
           #ret = reader.Read()
           #t_time.append(int(processtag(reader)))
        #elif reader.Name() == "type" and int(reader.NodeType()) == 1:
           #ret = reader.Read()
           #t_type.append(int(processtag(reader)))
        #elif reader.Name() == "tr" and int(reader.NodeType()) == 1:
          #ret = reader.Read()
          #TR = int(processtag(reader))
        #else :
           #ret = reader.Read()
    
    #if ret != 0:
        #print "%s : failed to parse" % (filename)
        
    #return (TR,t_type,t_time)
    
def streamFile(filename):
  from xml.dom import utils,core
  import string
  
  # Read an XML document into a DOM object
  reader = utils.FileReader(filename)
  
  # Retrieve top level DOM document object
  doc = reader.document
  TR=24
  t_time = []
  t_type = []
  #nodes = doc.documentElement.childNodes
  # Walk over the nodes
  for n in doc.documentElement.childNodes:
    if (n.nodeName == "type"):
      # Accumulate contents of text nodes
      t_type.append(int(n.xml))
      print "type:" + str(n.xml)
    elif (n.nodeName == "time"):
      t_time.append(int(n.xml))
      print "time:" + str(n.xml)
    elif (n.nodeName == "tr"):
      TR = int(n.xml)
      print "tr:" + str(n.xml)
    
  return (TR,t_type,t_time)

    
def execution( self, context ):
  assert (len(self.intmesh) == len(self.time_texture) and len(self.time_texture) == len(self.spmt_maps) and len(self.time_texture) == len(self.beta))
  from neurospy.bvfunc.surface_glm import CorticalGLM
  from soma import aims
  import numpy as N
  import sys,os,string

  reader2 = aims.Reader()
  texfilename = self.time_texture[0]
  texture = reader2.read(str(texfilename))
  nb_nodes = int(texture[0].nItem())
  nb_scans = int(texture.size())
 
  tab = N.arange(float(nb_nodes*nb_scans))
  k=0
  baseline = N.zeros(nb_scans)
  
  for i in range(0,nb_nodes):
     for j in range(0,nb_scans):
        tab[k] = texture[j][i]
        k=k+1 
        baseline[j] += texture[j][i]
     
  baseline /= nb_nodes
  
  tab = tab.reshape(nb_nodes,nb_scans)
  
  (TR,t_type,t_time) = streamFile(str(self.protocol_text))

  context.write("Finished")