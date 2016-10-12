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

name = 'Create an Average Volume from a 4D-Volume'

userLevel = 0

signature = Signature(
	  'input', ReadDiskItem('4D Volume', 'BrainVISA volume formats'),
	  'output', WriteDiskItem('3D Volume', 'BrainVISA volume formats')
)
def initialization( self ):
  pass

def execution( self, context ):
     import numpy as np
     from soma import aims
     ima = aims.read(self.input.fullPath())
     arr = np.array(ima, order='F')
     arr3 = np.mean(arr, 3)
     sx=ima.getSizeX()
     sy=ima.getSizeY()
     sz=ima.getSizeZ()
	
     ima_out = aims.Volume_FLOAT(sx,sy,sz,1)
     h=ima.header()
     ho=ima_out.header()

     ho['voxel_size']=[h['voxel_size'][0],h['voxel_size'][1],h['voxel_size'][2]]
     ho['volume_dimension']=[sx, sy, sz, 1]
     ho['transformations']=h['transformations']
     ho['referentials']=h['referentials']
	
     for z in range(sz):
	    for y in range(sy):
	         for x in range(sx):
	              ima_out.setValue(arr3[x,y,z], x, y, z)
     
     
     aims.write(ima_out, self.output.fullPath())
     context.write("Finished")
