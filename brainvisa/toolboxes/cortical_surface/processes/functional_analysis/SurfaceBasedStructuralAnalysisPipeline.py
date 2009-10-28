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

name = 'Surface-Based Structural Analysis Pipeline'
userLevel = 2

signature = Signature(  'intmesh', ListOf(ReadDiskItem( 'Hemisphere White Mesh', 'MESH mesh' )), 
        'surfacebased_data', ListOf(ReadDiskItem('Texture', 'Texture')),
        
        
        'primal_sketches', ListOf(WriteDiskItem( 'Graph', 'Graph')),
  )


def initialization( self ):
    eNode = SerialExecutionNode( self.name, parameterized=self )
    #self.setOptional('beta')
    
    eNode.addChild( 'SPMtMaps', ProcessExecutionNode( 'CreateSurfaceBasedSPMtMaps', optional = 1 ) )
    eNode.addChild( 'PrimalSketches', ProcessExecutionNode( 'CreateSurfaceBasedPrimalSketches', optional = 1 ) )
    eNode.addChild( 'GroupAnalysis', ProcessExecutionNode( 'PerformGroupAnalysis', optional = 1 ) )
    eNode.addChild( 'LabelsTexture', ProcessExecutionNode( 'CreateResultsLabelsTexture', optional = 1 ) )
    #eNode.addLink('Kernels.intmesh','intmesh')
    #eNode.addLink('Kernels.output','Projection.kernels')
    #eNode.addLink('FusionTextures.input', 'Projection.projection_texture')
    #eNode.addLink('Projection.white_mesh','intmesh')
    #eNode.addLink('Projection.functional_volumes','functional_volumes')
    #eNode.addLink('FusionTextures.output','timetexture')
    #eNode.addLink('StatisticalAnalysis.projection_texture','FusionTextures.output')
    #eNode.addLink('StatisticalAnalysis.spmt_texture','spmt_texture')
    #eNode.addLink('StatisticalAnalysis.contraste', 'contraste')
    #eNode.addLink('StatisticalAnalysis.beta', 'beta')
    #eNode.addLink('StatisticalAnalysis.protocol_text', 'protocol_text')
    
    self.setExecutionNode( eNode )
    
    
    
def execution (self, context):
  context.write("test")
