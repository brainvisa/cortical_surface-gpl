from neuroProcesses import *
import shfjGlobals     

name = 'Surface-Based Analysis Pipeline'
userLevel = 2

signature = Signature(  'intmesh', ReadDiskItem( 'Hemisphere White Mesh', 'MESH mesh' ), 
        'functional_volumes', ListOf(ReadDiskItem('4D Volume', 'BrainVISA volume formats')),
        'timetexture', WriteDiskItem('Texture', 'Texture'),
        'protocol_text', ReadDiskItem( 'Text File', 'Text File' ),
        'beta', WriteDiskItem('Texture', 'Texture'),
        'spmt_texture', WriteDiskItem( 'Texture', 'Texture'),
        'contraste', String()
  )


def initialization( self ):
    eNode = SerialExecutionNode( self.name, parameterized=self )
    self.setOptional('beta')
    
    eNode.addChild( 'Kernels', ProcessExecutionNode( 'CreateKernelsForProjection', optional = 1 ) )
    eNode.addChild( 'Projection', ProcessExecutionNode( 'FunctionalProjection', optional = 1 ) )
    eNode.addChild( 'FusionTextures', ProcessExecutionNode( 'TexturesAsTimeTexture', optional = 1 ) )
    eNode.addChild( 'StatisticalAnalysis', ProcessExecutionNode( 'SurfaceBasedSPMtMap', optional = 1 ) )
    eNode.addLink('Kernels.intmesh','intmesh')
    eNode.addLink('Kernels.output','Projection.kernels')
    eNode.addLink('FusionTextures.input', 'Projection.projection_texture')
    eNode.addLink('Projection.white_mesh','intmesh')
    eNode.addLink('Projection.functional_volumes','functional_volumes')
    eNode.addLink('FusionTextures.output','timetexture')
    eNode.addLink('StatisticalAnalysis.projection_texture','FusionTextures.output')
    eNode.addLink('StatisticalAnalysis.spmt_texture','spmt_texture')
    eNode.addLink('StatisticalAnalysis.contraste', 'contraste')
    eNode.addLink('StatisticalAnalysis.beta', 'beta')
    eNode.addLink('StatisticalAnalysis.protocol_text', 'protocol_text')
    
    self.setExecutionNode( eNode )