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

name = 'Projection of Functional Volumes'

userLevel = 2

signature = Signature(
      'white_mesh', ReadDiskItem( 'Hemisphere White Mesh', 'MESH mesh' ),
      'kernels', ReadDiskItem('Projection convolution kernels', 'GIS Image'),
      'functional_volumes', ListOf(ReadDiskItem('4D Volume', 'BrainVISA volume formats')),
      'projection_texture', ListOf(WriteDiskItem( 'Functional texture', 'Texture'))
)

#def volumemesh( values, process ):
    #result = None
    #if values.white_mesh is not None and values.functional_volumes is not None:
      #result = []
      #for functional in values.functional_volumes:
         #attributes = values.white_mesh.hierarchyAttributes()
         #attributes[ 'volume' ] = functional.name
         #result.append( process.signature[ 'projection_texture' ].contentType.findValue( attributes ) )
    #return result


def initialization( self ):
      self.linkParameters('kernels', 'white_mesh')
      self.linkParameters('projection_texture', 'white_mesh')

def execution( self, context ):
      i=0
      for volume in self.functional_volumes:         
         context.write('Creating a projection texture for ' + str(volume) + "...")
         projection = [ 'AimsFunctionProjection', 
            '-m', self.white_mesh,
            '-d', self.kernels,
            '-d1', volume,
            '-o', self.projection_texture[i],
            '-op', '1'
         ]
         i=i+1
         apply( context.system, projection)
      context.write('Adding mesh info in textures .minf files..')
      #for output in self.projection_texture:
         #reader = aims.Reader()
         #object = reader.read( str(output) )
         #print 'file:', object
         #h = object.header()
         #h['mesh']=str(self.white_mesh)
         #h['functional_volume']=str(self.functional_volumes[i])
         #writer = aims.Writer()
         #writer.write(object, str(output))
      context.write('Finished')

