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

name = 'Surfacic SPAM'
userLevel = 1

signature = Signature(
    'data_list', ListOf( ReadDiskItem( 'Texture' ,'Texture'  ) ),
    'threshold', Float(),
    'binarise',Choice("No","Yes"),
    'spam', WriteDiskItem( 'Texture' ,'Texture'),
    )

def initialization( self ):
   self.setOptional('threshold')


def execution( self, context ):

   nbTex = len(self.data_list)
   tex_temp = context.temporary( 'Texture')
   tex_sum = context.temporary( 'Texture')
   #Initialisation
   context.runProcess('TexLinearComb',
                      texture1 = (self.data_list)[0],
                      texture2 = (self.data_list)[0],
                      num1 = 0,
                      num2 = 0,
                      den1 = 0,
                      den2 = 0,
                      cst = 0,
                      output = tex_sum)
   
   #Summation of all the textures
   for tex in self.data_list:
      if self.binarise == 'Yes':
         context.runProcess('ThresholdTexture',
                            Texture = tex,
                            Thresholded_texture = tex_temp,
                            Mode = 'EQUAL_TO',
                            Binary = 'Yes',
                            threshold1 = self.threshold)
         context.runProcess('TexLinearComb',
                            texture1 = tex_sum,
                            texture2 = tex_temp,
                            num1 = 1,
                            num2 = 1,
                            den1 = 0,
                            den2 = 0,
                            cst = 0,
                            output = tex_sum)
      else:
         context.runProcess('TexLinearComb',
                            texture1 = tex_sum,
                            texture2 = tex,
                            num1 = 1,
                            num2 = 1,
                            den1 = 0,
                            den2 = 0,
                            cst = 0,
                            output = tex_sum)
   #Division by the number of tecture (mean texture)
   context.runProcess('TexLinearComb',
                        texture1 = tex_sum,
                        texture2 = tex_sum,
                        num1 = 1,
                        num2 = 1,
                        den1 = nbTex,
                        den2 = 0,
                        cst = 0,
                        output = self.spam )
   

                      
