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
import shfjGlobals     

name = '2 - Create Surface-Based Primal Sketch'
userLevel = 2

signature = Signature(  'intmesh', ListOf(ReadDiskItem( 'Hemisphere White Mesh', 'MESH mesh' )), 
        'surfacebased_activmap', ListOf(ReadDiskItem('Surface-Based SPMt Map', 'Texture')),
        
        
        'primal_sketch', ListOf(WriteDiskItem( 'Primal Sketch', 'Graph')),
          
        #'texture', ReadDiskItem( 'Texture','Texture'),
        #'white', ReadDiskItem( 'Hemisphere White Mesh', 'MESH mesh'),
        #'primalsketch', WriteDiskItem( 'Primal Sketch', 'Graph' ),
        'tMin', Float(),
        'tMax', Float(),
        'whiteAux', ListOf(ReadDiskItem( 'Hemisphere White Mesh', 'MESH mesh')),
        'filterout', Float(),
        'intersectioncriterium', Integer(),
        'latitude', ListOf(ReadDiskItem( 'Latitude coordinate texture','Texture')),
        'longitude', ListOf(ReadDiskItem('Longitude coordinate texture','Texture'))
  )


def initialization( self ):
     self.setOptional('whiteAux')
     self.linkParameters( 'surfacebased_activmap', 'intmesh' )
     self.linkParameters( 'primal_sketch', 'intmesh' )
     self.linkParameters( 'whiteAux', 'intmesh' )
     self.linkParameters( 'latitude', 'intmesh' )
     self.linkParameters( 'longitude', 'intmesh' )
     self.tMin = 1.0
     self.tMax = 64.0
     self.filterout = 2.0
     self.intersectioncriterium = 10    
    
def execution( self, context ):
     scales=context.temporary( 'Texture' )
     blobs=context.temporary( 'Texture' )
     assert(len(self.intmesh)==len(self.surfacebased_activmap) and len(self.intmesh)==len(self.primal_sketch))
     for i in xrange(len(self.intmesh)):
        call_list = [ 'AimsTexturePrimalSketch',
                      '-t', self.surfacebased_activmap[i],
                      '-o', self.primal_sketch[i],
                      '-m', self.intmesh[i],
                      '-os', scales,
                      '-ob', blobs,
                      '-t1', self.tMin,
                      '-t2', self.tMax,
          ]
        if (self.whiteAux[i] is not None):
          call_list += ['-mX', self.whiteAux[i]]
        s = self.intmesh[i].get( 'subject')
        assert(s!=None)
        call_list += ['-sj', s]
            
        if (self.filterout is not ""):
          call_list += ['-f', self.filterout]
        if (self.latitude is not None):
          call_list += ['-l', self.latitude[i]]
        if (self.longitude is not None):
          call_list += ['-L', self.longitude[i]]
        if (self.intersectioncriterium is not ""):
          call_list += ['-iP', self.intersectioncriterium]
    
    
        context.write('Starting primal sketch computation for subject ' + str(i))
        apply( context.system, call_list )
     context.write('Finished')