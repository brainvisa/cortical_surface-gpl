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

from __future__ import absolute_import
from brainvisa.processes import *
from soma import aims
import sys
import os
import numpy
import operator
import six

name = 'Average Sulci'
userLevel = 0

##################################################################
def sulcus2vector(mesh) :
     # met les noeuds du maillage dans un vecteur de type 
     # [x1,x2, .., xN, y1, ..., yN, z1, ... , zN]
     nodes=mesh.vertex(0)
     v=[]
     nn=nodes.size()
     for i in six.moves.xrange(3) :
          for j in six.moves.xrange(nn) :
               v.append(nodes[j][i])
     return v
          
##################################################################
def vector2vertices(coord) :

     # fait l'inverse de sulcus2vector et renvoie le vecteur de vertices 
     # d'un maillage
     N=len(coord)/3
     vertOut=[ ]
     for i in six.moves.xrange(N) :
          p=aims.vector_POINT3DF()
          p=[coord[i], coord[i+N], coord[i+2*N]]
          vertOut.append(p)
     return vertOut
     
##################################################################
def pcaSulci(sulci) :
     from math import sqrt
     #ici sulci doit etre un vecteur de vecteurs definis par sulcus2vector
     s=len(sulci) #nb sulci
     N=len(sulci[0]) #nb de points dans un sulcus (x3)
     mean=[]
     
     # on calcule le vecteur moyen et on centre les vecteurs de sillons
     for i in six.moves.xrange(N) :
          a=0
          for j in six.moves.xrange(s) :
               a=a+sulci[j][i]
          m=a/s
          mean.append(m)
          for j in six.moves.xrange(s) :
               sulci[j][i]=sulci[j][i] - m

     D=numpy.asmatrix(sulci)
     
     # la matrice de covariance est (Dt.D)/s qui est de taille enorme.
     # Ici on va calculer D.Dt qui nous permettra de diagonaliser plus facilement
     # car de taille s.s, puis de recuperer les valeurs et vecteurs propres de la
     # matrice de covariance
     T=numpy.dot(D, D.T)
     T=T/s
     val, vect=numpy.linalg.eig(T)
     
     # val contient maintenant les s premieres valeurs propres de la matrice de covariance
     # (dont les autres sont nulles car de rang <=s). Les vecteurs propres sont dans
     # Dt.vector
     
     eigenv=numpy.dot(D.T, vect)
     evect=numpy.asarray(eigenv.T)
     
     for i in six.moves.xrange(len(evect)) :
          tot=0
          for j in six.moves.xrange(len(evect[i])) :
               tot += evect[i][j]*evect[i][j]
          if (tot>0) :
               for j in six.moves.xrange(len(evect[i])) :
                    evect[i][j] /= sqrt(tot)
     
     return mean, val, evect
               
####################################################################

signature = Signature(
    'sulci', ListOf(ReadDiskItem('Mesh', 'MESH mesh' )),
    'mean', WriteDiskItem('Mesh', 'MESH mesh' ),
    'dev', WriteDiskItem('Texture', 'Texture' ),
    'save_sequence', Boolean(),
)

def initialization( self ):
     pass

def execution ( self, context ):

     from math import sqrt
     reader = aims.Reader()
     ns=len(self.sulci)
     sulcus=[]
     v=[]
     av=[]
     context.write("Found ", ns, " sulci")
     for i in six.moves.xrange(ns) :
          context.write("Reading ", self.sulci[i].fullPath())
          sulcus.append(reader.read(self.sulci[i].fullPath()))
          v.append(sulcus2vector(sulcus[i]))
     context.write("Calling PCA")
     meanS, valS, vectS=pcaSulci(v)
     context.write("PCA done, reconstructing vectors and saving")
     #m = aims.AimsTimeSurface_3_VOID()
     #vertM=vector2vertices(meanS)
     #m.vertex(0).assign(vertM)
     #m.polygon(0).assign(sulcus[0].polygon(0))
     #m.updateNormals()
     vertM=vector2vertices(meanS)
     nv=len(vertM)
     vmap=sulcus[0].vertex(0)
     for i in six.moves.xrange(nv) :
          vmap[i][0]=vertM[i][0]
          vmap[i][1]=vertM[i][1]
          vmap[i][2]=vertM[i][2]
     sulcus[0].vertex(0).assign(vmap)
     sulcus[0].updateNormals()
     W = aims.Writer()
     W.write(sulcus[0], self.mean.fullPath())
     context.write("Computing first mode of variation")
     imax=0
     lmax=0
     ltot=0
     valS2=[ ]
     for i in six.moves.xrange(len(valS)) :
          valS2.append((i,valS[i]))
     valS2_sorted=sorted(valS2, key=operator.itemgetter(1), reverse=True)
     context.write('List of sorted eigenvalues with sorted ', valS2_sorted)
     valTot=0
     for i in six.moves.xrange(len(valS)) :
          valTot+=valS[i]
     i=0;
     eigval=0.0
     percent=0.0
     while (percent<80.0) :
          eigval+=valS2_sorted[i][1]
          percent=eigval*100.0/valTot
          i+=1
     imax=i;
     vmax=vectS[imax]
     vertMax=vector2vertices(vmax)
     modeS=aims.AimsTimeSurface_3_VOID()
     context.write("Computing variations")

     for k in six.moves.xrange(imax) :
          context.write("mode ", k+1)
          eigval=valS2_sorted[k][1]/valTot*100.0
          context.write("Eigenvalue :", valS2_sorted[k][1], " -> ", valS2_sorted[k][0],  " -> ", eigval, "%")
          vmax=vectS[valS2_sorted[k][0]]
          vertMax=vector2vertices(vmax)
          for i in six.moves.xrange(200) :
               dt=(i-100)*5.0
               modeS.polygon(i).assign(sulcus[0].polygon(0))
               modeS.vertex(i).assign(sulcus[0].vertex(0))
               vmap=modeS.vertex(i)
               for j in six.moves.xrange(nv) :
                    vmap[j][0]=vmap[j][0] + dt*vertMax[j][0]
                    vmap[j][1]=vmap[j][1] + dt*vertMax[j][1]
                    vmap[j][2]=vmap[j][2] + dt*vertMax[j][2]
          modeS.updateNormals()
          nameOut="/home/olivier/variation_%i.mesh"%(k+1)
          context.write("Saving", nameOut)
          WM = aims.Writer()
          WM.write(modeS, nameOut)
     context.write("OK")

     



