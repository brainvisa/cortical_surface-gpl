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

name = 'Random Effect'
userLevel = 2

signature = Signature(
  'individual_maps', ReadDiskItem( 'Texture', 'Texture'),
  #'white_mesh', ReadDiskItem( 'Hemisphere White Mesh', 'MESH mesh' ),
  'rfx_map', WriteDiskItem('Texture', 'Texture')
  )
def initialization( self ):
  pass
     

def execution( self, context ):
  import fff.glm as G
  import numpy as N
  from soma import aims
  import sys,os,string 

  #########################
  reader = aims.Reader()
  texture = reader.read(str(self.individual_maps))
  context.write(texture.size())

  nb_nodes = int(texture[0].nItem())
  nb_sujets = int(texture.size())
 
  tab = N.arange(float(nb_nodes*nb_sujets))
  k=0
  baseline = N.zeros(nb_sujets)
  #moyenne = 0.0
  
  for i in range(0,nb_nodes):
     for j in range(0,nb_sujets):
        tab[k] = texture[j][i]
        k=k+1 
        baseline[j] += texture[j][i]
        #moyenne += texture[j][i]
     
  baseline /= nb_nodes
  #moyenne /= nb_nodes
  #moyenne /= nb_sujets
  tab = tab.reshape(nb_nodes,nb_sujets)
  nb_cond = int(0)

  reg = N.zeros(nb_sujets, float)
  reg = reg.reshape(nb_sujets, 1)
   
  reg[:,0] = baseline
  #print reg[:,0]

  #B, nVB, s2, dof = G.KF_fit(tab, reg, axis=1)

  #model = G.glm(tab, reg, axis=1)

  B, nVB, s2, a, dof = G.RKF_fit(tab, reg, axis=1)
  #print "dofDOOA"
  context.write(dof)
  #model.fit(method='kalman.ar1') ## default is 'ols'
  c = N.zeros(int(nb_cond+1))
  j=0
  c[j] = 1
#c = array([1,0])
  cB, VcB, t = G.t_contrast(c, B, nVB, s2, axis=1)
  #tcon = model.contrast(c) ## recognizes a t contrast
  #t, p, z = tcon.test(zscore=True)
  writer = aims.Writer()
  textur = aims.TimeTexture_FLOAT()
  textur[0] = aims.Texture_FLOAT(int(t.size))
  
  for i in range(0,t.size):
     textur[0][i] = float(t[i])
  writer.write(textur, str(self.rfx_map))

      
