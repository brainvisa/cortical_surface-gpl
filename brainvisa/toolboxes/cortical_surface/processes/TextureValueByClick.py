# Copyright CEA and IFR 49 (2000-2005)
#
#  This software and supporting documentation were developed by
#      CEA/DSV/SHFJ and IFR 49
#      4 place du General Leclerc
#      91401 Orsay cedex
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
from math import *

name = 'Texture Value By Click'
userLevel = 2

signature = Signature(
  'mesh', ReadDiskItem( 'Mesh', 'Aims mesh formats' ),
  'texture', ReadDiskItem('Texture', 'Texture'),
  'point', Point3D()
  )
def initialization( self ):
  self.signature[ 'point' ].add3DLink( self, 'mesh' )
  

def execution( self, context ):
  from soma import aims
  import sys,os,string
  point = []
  point = self.point
  reader = aims.Reader()
  texture = reader.read(str(self.texture))
  context.write(str(texture.size()) + ' ' + str(len(texture[0])))
  mesh = reader.read(str(self.mesh))
  distance = 0.0
  mini = 0
  distmini=1000.0
  for t in xrange(mesh.size()):
    for v in xrange(len(mesh.vertex(t))):
      vert = mesh.vertex(t)[v]
      distance = sqrt(pow(vert[0]-point[0],2)+pow(vert[1]-point[1],2)+pow(vert[2]-point[2],2))
      if (distance < distmini):
        distmini = distance
        mini = v
  context.write("L'indice du noeud le plus proche du point clique est : " + str(mini))
  context.write("La distance calculee entre ce noeud et le point clique est : " + str(distmini))
  context.write("Les valeurs de ce noeud dans la texture donnee sont : ")
  for t in xrange(texture.size()):
    context.write(texture[t][mini])
  
  

      
