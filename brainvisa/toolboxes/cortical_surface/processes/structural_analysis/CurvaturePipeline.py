from neuroProcesses import *
import shfjGlobals     

name = 'Curvature Pipeline'
userLevel = 2

signature = Signature(  'mesh', ReadDiskItem( 'Mesh', 'MESH mesh' ), 
  'texture', ReadDiskItem( 'Texture', 'Texture' ), 
    'graph', WriteDiskItem( 'Curvature Blobs Graph', 'Graph and data' ), 
  )


def initialization( self ):
    eNode = SerialExecutionNode( self.name, parameterized=self )
    self.setOptional('mesh','texture','graph')
    

    eNode.addChild( 'Smoothing', ProcessExecutionNode( 'SmoothTexture', optional = 1 ) )
    eNode.addChild( 'Thresholding', ProcessExecutionNode( 'ThresholdTexture', optional = 1 ) )
    eNode.addChild( 'Converting', ProcessExecutionNode( 'ConvertToFLOATTexture', optional = 1 ) )
    eNode.addChild( 'Dilating', ProcessExecutionNode( 'DilateTexture', optional = 1 ) )
    eNode.addChild( 'Eroding', ProcessExecutionNode( 'ErodeTexture', optional = 1 ) )
    eNode.addChild( 'ConnectedComponents', ProcessExecutionNode( 'ConnectedComp', optional = 1 ) )
    eNode.addChild( 'GraphCreating', ProcessExecutionNode( 'ExtractBlobsFromLabelTexture', optional = 1 ) )
    eNode.addLink('texture','mesh')
    eNode.addLink('graph','mesh')

    eNode.addLink('Smoothing.texture','texture')
    eNode.addLink('Smoothing.mesh','mesh')
    #eNode.addLink('Smoothing.outfile','texture')
    eNode.addLink('GraphCreating.graph','graph')
    eNode.addLink('Thresholding.input','Smoothing.outfile')
    eNode.addLink('Thresholding.output','Thresholding.input')
    eNode.addLink('Converting.input','Thresholding.output')
    eNode.addLink('Dilating.input','Converting.output')

    eNode.addLink('Eroding.input','Dilating.output')
    eNode.addLink('Dilating.mesh','mesh')
    eNode.addLink('Dilating.output','Dilating.input')
    eNode.addLink('ConnectedComponents.input','Dilating.output')
    eNode.addLink('ConnectedComponents.mesh','mesh')
    eNode.addLink('GraphCreating.texture','Dilating.output')
    eNode.addLink('GraphCreating.white','mesh')
    
    
    

    
    self.setExecutionNode( eNode )