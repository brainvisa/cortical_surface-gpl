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

from brainvisa.processes import *
import re, os, sys
from soma import aims


name = 'Get Sulci Mesh From Graph'
userLevel = 1
	
################################################################################
# Signature
################################################################################
     
signature = Signature(
     'graphname', ReadDiskItem( 'Cortical folds graph', 'Graph'),
     'sulci', ListOf( String() ),
     'sulcus_identification', Choice( 'label', 'name' ),
     'sulci_mesh_pattern', WriteDiskItem( 'Sulcus mesh', 'Aims mesh formats' ),
     )

################################################################################
# Init of the process
################################################################################

def initialization( self ):
     self.sulcus_identification = 'label'


################################################################################
# Functions specific to this process 
################################################################################

TOP_STATE = 0
GRAPH_STATE = 1
NODE_STATE = 2
UNKNOWN_STATE = 255


def read_graph(context, graphname, sulci, ident):
     context.write('Entering read_graph with ', graphname, ' and ', sulci)
     tag_map = {
          'GRAPH' : GRAPH_STATE,
          'NODE' : NODE_STATE,
     }
     states = [TOP_STATE]
     fd = open(graphname)
     sulcus = 'none'
     ind = -1
     indices = []
     for line in fd.readlines():
          m = re.match('\*BEGIN (\w*) (.*)\n', line)
          if m:
               category, dummy = m.groups()
               states.append(tag_map.get(category, UNKNOWN_STATE))
               continue
          m = re.match('\*END', line)
          if m:
               if states[-1] == NODE_STATE:
                    if sulcus in sulci: indices.append(int(ind))
                    sulcus = 'none'
               states.pop()
               continue
          if states[-1] == GRAPH_STATE:
               m = re.match('filename_base[ \t]*([\*\w\.]*)', line)
               if m:
                    data_dir = m.groups()[0]
                    continue
               m = re.match('fold.global.tri[ \t]*([\w]*) (\w*) (\w*)',
                                             line)
               if m:
                    dummy1, tmtktri, dummy2 = m.groups()
                    continue
          elif states[-1] == NODE_STATE:
               m = re.match('Tmtktri_label[ \t]*(\w*)', line)
               if m:
                    ind = m.groups()[0]
                    continue
               if ident == 'name':
                    m = re.match('name[ \t]*([\w.]*)', line)
               if ident == 'label':
                    m = re.match('label[ \t]*([\w.]*)', line)
               if m:
                    sulcus = m.groups()[0]
                    continue
     fd.close()
     if data_dir == '*':
          data_dir = re.sub('\.arg$', '.data',
               os.path.basename(graphname))
     meshname = os.path.join(data_dir, tmtktri + '.mesh')
     context.write('Exiting read_graph with ', meshname, ' and ', indices)
     return meshname, indices



################################################################################
# Execution of the process 
################################################################################

def execution( self, context ):
     sulcus_mesh=aims.AimsSurfaceTriangle()
     for sul in self.sulci:
          meshname, indices = read_graph(context, self.graphname.fullPath(), sul, self.sulcus_identification)
          meshname = os.path.join(os.path.dirname(self.graphname.fullName()), meshname)
          context.write(self.graphname.get( 'subject' ))
          #context.write(self.graphname.subject)
          context.write(meshname)
          graphmesh = aims.Reader().read(meshname)
          for ind in indices:
               new_mesh = aims.AimsSurfaceTriangle()
               v = graphmesh.vertex(ind)
               new_mesh.vertex().assign(v)
               n = graphmesh.normal(ind)
               new_mesh.normal().assign(n)
               p = graphmesh.polygon(ind)
               new_mesh.polygon().assign(p)
               aims.SurfaceManip.meshMerge(sulcus_mesh, new_mesh)
     writer = aims.Writer()
     fileName='%s_%s.mesh' %(self.sulci_mesh_pattern.fullName(), self.graphname.get( 'subject' ) )
     writer.write(sulcus_mesh, fileName)
          
          
          
          
          
          
          
          
          
          
 
