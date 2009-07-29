

from neuroProcesses import *
import shfjGlobals     

name = 'Cortical Surface Parameterization Pipeline'
userLevel = 2

signature = Signature( 
  'Lgraph', ReadDiskItem( 'Cortical folds graph', 'Graph',requiredAttributes={ 'side': 'left' } ),
  'Rgraph', ReadDiskItem( 'Cortical folds graph', 'Graph',requiredAttributes={ 'side': 'right' } ),
  'sulcus_identification',Choice('name','label')
  )


def initialization( self ):
    self.linkParameters( 'Lgraph','Rgraph')
    self.linkParameters( 'Rgraph','Lgraph')
    self.sulcus_identification='label'
    eNode = SerialExecutionNode( self.name, parameterized=self )

    eNode.addChild( 'LeftHemisphere_Process',
                    ProcessExecutionNode( 'LeftHemisphereProcess', optional = 1 ) )
    eNode.addChild( 'RightHemisphere_Process',
                    ProcessExecutionNode( 'RightHemisphereProcess', optional = 1 ) )
    

    eNode.addLink( 'LeftHemisphere_Process.Lgraph', 'Lgraph' )
    eNode.addLink( 'RightHemisphere_Process.Rgraph', 'Rgraph' )

    eNode.addLink( 'LeftHemisphere_Process.ConstraintProjectionLeft.sulcus_identification', 'sulcus_identification')
    eNode.addLink( 'RightHemisphere_Process.ConstraintProjectionRight.sulcus_identification', 'sulcus_identification')
    
    self.setExecutionNode( eNode )
