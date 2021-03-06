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
from numpy import *
from six.moves import range


name = 'Average Functional information in Gyri'

userLevel = 0

signature = Signature(
     'GyriTexture', ReadDiskItem( 'Texture', 'Texture' ),
     'FunctionalTexture',ReadDiskItem( 'Texture', 'Texture' ),
#     'ResultTextFile', WriteDiskItem( 'Text File', 'Text File')
)

def initialization( self ):
     pass

def execution( self, context ):
     nameOut='%s_gyri.txt' %(self.FunctionalTexture.fullName())
     context.write('The results will be written in ',nameOut)
     readerF = aims.Reader()
     texF=readerF.read(str(self.FunctionalTexture) )

     ni=texF.nItem()
     nt=texF.header()['nb_t_pos']

     mat=array(texF[0])
     for i in range(1,nt):
          mat=vstack((mat, array(texF[i])))

     matT=mat.transpose()

     readerG = aims.Reader()
     gyriTex=readerG.read(str(self.GyriTexture) )
     gyri=array(gyriTex[0])
     nlab=gyri.max()

     fileO=open(nameOut ,'w')
#     fileO=open(str(self.ResultTextFile) ,'w')
     for i in range(0, int(nlab+1)):
          gyrus=matT[where( gyri==i )]
          if (gyrus.size != 0):
               #print str(i)+' '
               context.write(str(i))
               meanG=mean(gyrus, axis=0)
               fileO.write(str(i)+'\n')
               sig = ""
               for j in meanG:
                    sig = sig + str(j) + ' '
               fileO.write(sig+'\n')

     fileO.close()
