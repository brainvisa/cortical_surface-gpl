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
import shfjGlobals, neuroDiskItems, neuroHierarchy, neuroConfig

name = 'Creation of Kernels for Volumes Projection'

userLevel = 2

signature = Signature(
      'intmesh', ReadDiskItem( 'Hemisphere White Mesh', 'MESH mesh' ), 
      'output', WriteDiskItem('Projection convolution kernels', 'GIS Image'),
      'size', Integer(),
      'resolutionX', Float(),
      'resolutionY', Float(),
      'resolutionZ', Float(),
      'geod_decay', Float(),
      'norm_decay', Float()
)
def initialization( self ):
      self.linkParameters('output', 'intmesh')
      self.size = 7
      self.resolutionX=3.0
      self.resolutionY=3.0
      self.resolutionZ=3.0
      self.geod_decay=5.0
      self.norm_decay=2.0
     
     

def execution( self, context ):
      from soma import aims
      context.write('Calculating the kernels...')
      projection = [ 'AimsFunctionProjection',
         '-m', self.intmesh,
         '-o', self.output,
         '-i', self.size,
         '-vx', self.resolutionX,
         '-vy', self.resolutionY,
         '-vz', self.resolutionZ,
         '-g', self.geod_decay,
         '-d', self.norm_decay,
         '-op', '0',
         '-t', 0
      ]
      apply( context.system, projection)
      reader = aims.Reader()
      object = reader.read( str(self.output) )
      print 'file:', object
      h = object.header()
      h['geod_decay']=str(self.geod_decay) 
      h['norm_decay']=str(self.norm_decay) 
      writer = aims.Writer()
      writer.write(object, str(self.output))
      context.write('Finished')
      

