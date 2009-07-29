from neuroProcesses import *

name = 'Left hemisphere'
userLevel = 2
 
signature = Signature(
  'Lgraph', ReadDiskItem( 'Cortical folds graph', 'Graph', requiredAttributes={ 'side': 'left' }  )
)

def initialization( self ):
    eNode = SerialExecutionNode( self.name, parameterized=self )

    eNode.addChild( 'ChangeTemplateReferentialLeft',
                    ProcessExecutionNode( 'ChangeTemplateReferentialLeft', optional = 1 ) )
    eNode.addChild( 'CingularPoleLeft',
                    ProcessExecutionNode( 'CingularPoleProjectionLeft', optional = 1 ) )
    eNode.addChild( 'ConstraintProjectionLeft',
                    ProcessExecutionNode( 'CorticalConstraintsLeft', optional = 1 ) )
    eNode.addChild( 'ConstraintCleanerLeft',
                    ProcessExecutionNode( 'ConstraintCleanerLeft', optional = 1 ) )
    eNode.addChild( 'CorticalSurfaceParameterizationLeft',
                    ProcessExecutionNode( 'ParameterizeHemisphereLeft', optional = 1 ) )
    eNode.addChild( 'CorticalSurfaceParcellationLeft',
                    ProcessExecutionNode( 'ParcellationLeft', optional = 1 ) )

    eNode.addLink( 'ChangeTemplateReferentialLeft.mri_corrected', 'Lgraph' )
    
    eNode.addLink( 'ConstraintProjectionLeft.Lgraph', 'Lgraph' )
    
    eNode.addLink( 'ChangeTemplateReferentialLeft.transformation_input', 'Lgraph' )
    
    eNode.addLink(  'CingularPoleLeft.left_pole_template','ChangeTemplateReferentialLeft.output_template' )
    
    eNode.addLink(  'CingularPoleLeft.left_white_mesh','ConstraintProjectionLeft.Lgraph' )
    
    eNode.addLink(  'CorticalSurfaceParameterizationLeft.left_white_mesh','ConstraintProjectionLeft.Lgraph' )
    
    eNode.addLink(  'CorticalSurfaceParameterizationLeft.left_cingular_pole', 'CingularPoleLeft.left_pole')
    
    eNode.addLink(  'ConstraintCleanerLeft.left_white_mesh', 'CingularPoleLeft.left_white_mesh')
    
    eNode.addLink(  'ConstraintCleanerLeft.left_white_sulci_mer', 'ConstraintProjectionLeft.left_white_sulci_mer')
    
    eNode.addLink(  'ConstraintCleanerLeft.left_white_sulci_par', 'ConstraintProjectionLeft.left_white_sulci_par')
    
    eNode.addLink(  'ConstraintCleanerLeft.left_cingular_pole', 'CingularPoleLeft.left_pole')
    
    eNode.addLink(  'ConstraintCleanerLeft.left_sulci_label_to_sulci_name', 'ConstraintProjectionLeft.left_sulci_label_to_sulci_name' )
    
    eNode.addLink(  'CorticalSurfaceParameterizationLeft.left_poles_texture', 'ConstraintCleanerLeft.left_poles_texture')
    
    eNode.addLink( 'CorticalSurfaceParameterizationLeft.left_white_sulci_par','ConstraintCleanerLeft.left_white_sulci_par_cleaned')
    eNode.addLink( 'CorticalSurfaceParameterizationLeft.left_white_sulci_mer','ConstraintCleanerLeft.left_white_sulci_mer_cleaned')
    
    eNode.addLink( 'CorticalSurfaceParcellationLeft.left_longitude','CorticalSurfaceParameterizationLeft.left_longitude')
    
    self.setExecutionNode( eNode )
