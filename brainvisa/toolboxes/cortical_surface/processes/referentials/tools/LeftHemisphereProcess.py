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
    eNode.addChild( 'RegularizeParcellationLeft',
                    ProcessExecutionNode( 'GyriRegularizationLeft', optional = 1 ) )

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
    
    eNode.addLink( 'RegularizeParcellationLeft.left_white_mesh','CorticalSurfaceParameterizationLeft.left_white_mesh')
    
    self.setExecutionNode( eNode )
