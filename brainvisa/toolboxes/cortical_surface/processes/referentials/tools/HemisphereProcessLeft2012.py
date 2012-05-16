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

name = 'Pipeline 2012 parameterization left hemisphere'
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

    eNode.addChild( 'SulcalinesExtractionLeft',
                    ProcessExecutionNode( 'SulcalinesExtractionLeft', optional = 1 ) )

    eNode.addChild( 'ParameterizeUnconstrainedHarmonicLeft',
                    ProcessExecutionNode( 'ParameterizeUnconstrainedHarmonicLeft', optional = 1 ) )
        
    eNode.addChild( 'HarmonicMappingOrthoLeft',
                    ProcessExecutionNode( 'HarmonicMappingOrthoLeft', optional = 1 ) )
        
    eNode.addLink( 'ChangeTemplateReferentialLeft.mri_corrected', 'Lgraph' )
    eNode.addLink( 'ChangeTemplateReferentialLeft.transformation_input', 'Lgraph' )
    
    eNode.addLink(  'CingularPoleLeft.left_pole_template','ChangeTemplateReferentialLeft.output_template' )
    eNode.addLink(  'CingularPoleLeft.left_white_mesh','SulcalinesExtractionLeft.Lwhite_mesh')
    eNode.addLink(  'CingularPoleLeft.left_white_mesh','ParameterizeUnconstrainedHarmonicLeft.Lwhite_mesh')
    eNode.addLink(  'CingularPoleLeft.left_white_mesh','HarmonicMappingOrthoLeft.Lwhite_mesh')

    eNode.addLink(  'SulcalinesExtractionLeft.Lwhite_mesh','CingularPoleLeft.left_white_mesh')
    eNode.addLink(  'SulcalinesExtractionLeft.Lgraph','Lgraph' )        
    
    eNode.addLink(  'ParameterizeUnconstrainedHarmonicLeft.Lwhite_mesh','SulcalinesExtractionLeft.Lwhite_mesh')
    eNode.addLink(  'HarmonicMappingOrthoLeft.Lwhite_mesh','SulcalinesExtractionLeft.Lwhite_mesh')
    

    self.setExecutionNode( eNode )
