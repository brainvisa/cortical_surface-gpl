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
        'surfacebased_data', ListOf(ReadDiskItem('Functional Time Texture', 'Texture')),
        'surfacebased_SPMt_maps', ListOf(WriteDiskItem('Surface-Based SPMt Map', 'Texture')),
        'primal_sketches', ListOf(WriteDiskItem( 'Primal Sketch', 'Graph')),
        'contrast', String(),
        'contrast_name', String(),
        'protocol_text', ReadDiskItem( 'Text File', 'Text File' )
        
          
  )



def initialization( self ):
    eNode = SerialExecutionNode( self.name, parameterized=self )

    
    eNode.addChild( 'SPMtMaps', ProcessExecutionNode( 'CreateSurfaceBasedSPMtMaps', optional = 1 ) )
    eNode.addChild( 'PrimalSketches', ProcessExecutionNode( 'CreateSurfaceBasedPrimalSketches', optional = 1 ) )
    eNode.addChild( 'GroupAnalysis', ProcessExecutionNode( 'PerformGroupAnalysis', optional = 1 ) )
    eNode.addChild( 'Significance', ProcessExecutionNode( 'ResultsSignificance', optional = 1 ) )
    eNode.addChild( 'LabelsTexture', ProcessExecutionNode( 'CreateResultsLabelsTexture', optional = 1 ) )
    
    
    #eNode.addLink('SPMtMaps.meshes', 'intmesh')
    eNode.addLink('surfacebased_SPMt_maps', 'SPMtMaps.spmtmaps')
    eNode.addLink('PrimalSketches.intmesh', 'intmesh')
    eNode.addLink('intmesh', 'PrimalSketches.intmesh')
    eNode.addLink('SPMtMaps.contrast', 'contrast')
    eNode.addLink('contrast','SPMtMaps.contrast')
    eNode.addLink('SPMtMaps.contrast', 'contrast_name')
    eNode.addLink('contrast_name','SPMtMaps.contrast')
    eNode.addLink('surfacebased_data', 'SPMtMaps.boldtextures')
    eNode.addLink('primal_sketches','PrimalSketches.primal_sketch')
    eNode.addLink('GroupAnalysis.primalsketches','PrimalSketches.primal_sketch')
    eNode.addLink('GroupAnalysis.labeled_primalsketches','PrimalSketches.primal_sketch')
    eNode.addLink('SPMtMaps.protocolfile', 'protocol_text')
    eNode.addLink('SPMtMaps.protocolfile', 'protocol_text')
    eNode.PrimalSketches.removeLink( 'surfacebased_activmap', 'intmesh' )
    eNode.addLink('PrimalSketches.surfacebased_activmap','SPMtMaps.spmtmaps')
    eNode.addLink('Significance.ddx1', 'GroupAnalysis.ddx1')
    eNode.addLink('Significance.ddx2', 'GroupAnalysis.ddx2')
    eNode.addLink('Significance.simweight', 'GroupAnalysis.simweight')
    eNode.addLink('Significance.datadrivenweight', 'GroupAnalysis.datadrivenweight')
    eNode.addLink('Significance.intrapsweight', 'GroupAnalysis.intrapsweight')
    eNode.addLink('Significance.lsweight', 'GroupAnalysis.lsweight')
    eNode.addLink('Significance.ddh', 'GroupAnalysis.ddh')
    eNode.addLink('Significance.labeled_primalsketches', 'GroupAnalysis.labeled_primalsketches')
    
    

    self.setExecutionNode( eNode )
    

