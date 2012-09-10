# -*- coding: utf-8 -*-
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

name = 'Pipeline 2012 parameterization right hemisphere'
userLevel = 2
 
signature = Signature(
  'Rgraph', ReadDiskItem( 'Labelled Cortical folds graph', 'Graph', requiredAttributes={ 'side': 'right' }  )
)

def initialization( self ):
    eNode = SerialExecutionNode( self.name, parameterized=self )

    eNode.addChild( 'ChangeTemplateReferential',
                    ProcessExecutionNode( 'ChangeTemplateReferential', optional = 1 ) )
    eNode.addChild( 'CingularPoleRight',
                    ProcessExecutionNode( 'CingularPoleProjectionRight', optional = 1 ) )

    eNode.addChild( 'SulcalinesExtractionRight',
                    ProcessExecutionNode( 'SulcalinesExtractionRight', optional = 1 ) )

    eNode.addChild( 'ParameterizeUnconstrainedHarmonicRight',
                    ProcessExecutionNode( 'ParameterizeUnconstrainedHarmonicRight', optional = 1 ) )
        
    eNode.addChild( 'HarmonicMappingOrthoRight',
                    ProcessExecutionNode( 'HarmonicMappingOrthoRight', optional = 1 ) )
        
    eNode.ChangeTemplateReferential.findValue( 'pole_template',
        { 'side' : 'right' } )

    eNode.addLink( 'ChangeTemplateReferential.mri_corrected', 'Rgraph' )
    eNode.addLink( 'ChangeTemplateReferential.transformation_input', 'Rgraph' )
    
    eNode.addLink(  'CingularPoleRight.right_pole_template','ChangeTemplateReferential.output_template' )
    eNode.addLink(  'CingularPoleRight.right_white_mesh','SulcalinesExtractionRight.Rwhite_mesh')
    eNode.addLink(  'CingularPoleRight.right_white_mesh','ParameterizeUnconstrainedHarmonicRight.Rwhite_mesh')
    eNode.addLink(  'CingularPoleRight.right_white_mesh','HarmonicMappingOrthoRight.Rwhite_mesh')

    eNode.addLink(  'SulcalinesExtractionRight.Rwhite_mesh','CingularPoleRight.right_white_mesh')
    eNode.addLink(  'SulcalinesExtractionRight.Rgraph','Rgraph' )        
    
    eNode.addLink(  'ParameterizeUnconstrainedHarmonicRight.Rwhite_mesh','SulcalinesExtractionRight.Rwhite_mesh')
    eNode.addLink(  'HarmonicMappingOrthoRight.Rwhite_mesh','SulcalinesExtractionRight.Rwhite_mesh')
    

    self.setExecutionNode( eNode )
