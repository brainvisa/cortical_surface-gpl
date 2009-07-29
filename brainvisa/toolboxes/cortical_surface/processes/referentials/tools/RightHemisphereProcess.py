from neuroProcesses import *

name = 'Right hemisphere'
userLevel = 2
 
signature = Signature(
  'Rgraph', ReadDiskItem( 'Cortical folds graph', 'Graph', requiredAttributes={ 'side': 'right' } )
)

def initialization( self ):
    eNode = SerialExecutionNode( self.name, parameterized=self )

    eNode.addChild( 'ChangeTemplateReferentialRight',
                    ProcessExecutionNode( 'ChangeTemplateReferentialRight', optional = 1 ) )
    eNode.addChild( 'CingularPoleRight',
                    ProcessExecutionNode( 'CingularPoleProjectionRight', optional = 1 ) )
    eNode.addChild( 'ConstraintProjectionRight',
                    ProcessExecutionNode( 'CorticalConstraintsRight', optional = 1 ) )
    eNode.addChild( 'ConstraintCleanerRight',
                    ProcessExecutionNode( 'ConstraintCleanerRight', optional = 1 ) )
    eNode.addChild( 'CorticalSurfaceParameterizationRight',
                    ProcessExecutionNode( 'ParameterizeHemisphereRight', optional = 1 ) )
    eNode.addChild( 'CorticalSurfaceParcellationRight',
                    ProcessExecutionNode( 'ParcellationRight', optional = 1 ) )

    eNode.addLink( 'ChangeTemplateReferentialRight.mri_corrected', 'Rgraph' )

    eNode.addLink( 'ConstraintProjectionRight.Rgraph', 'Rgraph' )
    
    eNode.addLink( 'ChangeTemplateReferentialRight.transformation_input', 'Rgraph' )
    
    eNode.addLink(  'CingularPoleRight.right_pole_template','ChangeTemplateReferentialRight.output_template' )
    
    eNode.addLink(  'CingularPoleRight.right_white_mesh','ConstraintProjectionRight.Rgraph' )
    
    eNode.addLink(  'CorticalSurfaceParameterizationRight.right_white_mesh','ConstraintProjectionRight.Rgraph' )
    
    eNode.addLink(  'CorticalSurfaceParameterizationRight.right_cingular_pole', 'CingularPoleRight.right_pole')
    
    eNode.addLink(  'ConstraintCleanerRight.right_white_mesh', 'CingularPoleRight.right_white_mesh')
    
    eNode.addLink(  'ConstraintCleanerRight.right_white_sulci_mer', 'ConstraintProjectionRight.right_white_sulci_mer')
    
    eNode.addLink(  'ConstraintCleanerRight.right_white_sulci_par', 'ConstraintProjectionRight.right_white_sulci_par')
    
    eNode.addLink(  'ConstraintCleanerRight.right_cingular_pole', 'CingularPoleRight.right_pole')
    
    eNode.addLink(  'ConstraintCleanerRight.right_sulci_label_to_sulci_name', 'ConstraintProjectionRight.right_sulci_label_to_sulci_name' )
    
    eNode.addLink(  'CorticalSurfaceParameterizationRight.right_poles_texture', 'ConstraintCleanerRight.right_poles_texture')
    
    eNode.addLink( 'CorticalSurfaceParameterizationRight.right_white_sulci_par','ConstraintCleanerRight.right_white_sulci_par_cleaned')
    eNode.addLink( 'CorticalSurfaceParameterizationRight.right_white_sulci_mer','ConstraintCleanerRight.right_white_sulci_mer_cleaned')
    
    eNode.addLink( 'CorticalSurfaceParcellationRight.right_longitude','CorticalSurfaceParameterizationRight.right_longitude')
    
    self.setExecutionNode( eNode )
