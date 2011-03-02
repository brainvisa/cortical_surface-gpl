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
from soma import aims
from numpy import *


name = 'Get Functional Signals in Gyrus'

userLevel = 1

signature = Signature(
     'GyriTexture', ReadDiskItem( 'Texture', 'Texture' ), 
     'FunctionalTexture',ReadDiskItem( 'Texture', 'Texture' ),
     'GyrusIndex', Integer(),
     'Output', WriteDiskItem('Text file', 'Text file')
)

def initialization( self ):
     pass
     
def execution( self, context ): 
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

     context.write('Opening output file ', str(self.Output))
     fileO=open(str(self.Output) ,'w')
#     fileO=open(str(self.ResultTextFile) ,'w')
     gyrus=matT[where( gyri==self.GyrusIndex )]
     listI=where( gyri==self.GyrusIndex )[0]
     
     for i in range(listI.size):
          ind=listI[i]
          fileO.write(str(ind)+'\n')
          sig = ""
          for j in range(nt):
               sig = sig + str(gyrus[i][j]) + ' '
          fileO.write(sig+'\n')

     fileO.close()
