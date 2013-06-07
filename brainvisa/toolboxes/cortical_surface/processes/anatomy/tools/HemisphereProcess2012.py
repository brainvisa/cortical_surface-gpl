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

name = 'Hip-Hop Hemispheric  Parameterization'
userLevel = 2
 
signature = Signature(
  'side', Choice('left', 'right'),
  'graph', ReadDiskItem( 'Labelled Cortical folds graph', 'Graph'  ),
  'sulcus_identification', Choice('name','label')
)

def initialization( self ):
    eNode = SerialExecutionNode( self.name, parameterized=self )

    eNode.addChild( 'ChangeTemplateReferential',
                    ProcessExecutionNode( 'ChangeTemplateReferential', optional = 1 ) )
    eNode.addChild( 'CingularPole',
                    ProcessExecutionNode( 'CingularPoleProjection', optional = 1 ) )

    eNode.addChild( 'InsularPole',
                    ProcessExecutionNode( 'InsularPoleProjection', optional = 1 ) )

    eNode.addChild( 'SulcalinesExtraction',
                    ProcessExecutionNode( 'SulcalinesExtraction', optional = 1 ) )

    eNode.addChild( 'ParameterizeUnconstrainedHarmonic',
                    ProcessExecutionNode( 'ParameterizeUnconstrainedHarmonic', optional = 1 ) )
        
    eNode.addChild( 'HarmonicMappingOrtho',
                    ProcessExecutionNode( 'HarmonicMappingOrtho', optional = 1 ) )

    eNode.addLink( 'ChangeTemplateReferential.mri_corrected', 'graph' )
    eNode.addLink( 'ChangeTemplateReferential.transformation_input', 'graph' )
    eNode.addLink( 'ChangeTemplateReferential.side', 'side' )

    eNode.addLink( 'CingularPole.white_mesh', 'graph')
    eNode.addLink(  'SulcalinesExtraction.white_mesh','graph')
    eNode.addLink(  'ParameterizeUnconstrainedHarmonic.white_mesh','graph')
    eNode.addLink(  'HarmonicMappingOrtho.white_mesh','graph')

    eNode.addLink(  'CingularPole.pole_template','ChangeTemplateReferential.output_template' )
#    eNode.addLink(  'CingularPole.white_mesh','SulcalinesExtraction.white_mesh')
#    eNode.addLink(  'CingularPole.white_mesh','ParameterizeUnconstrainedHarmonic.white_mesh')
#    eNode.addLink(  'CingularPole.white_mesh','HarmonicMappingOrtho.white_mesh')
    eNode.addLink(  'CingularPole.side', 'side' )

    eNode.addLink( 'InsularPole.side', 'side')
    eNode.addLink( 'InsularPole.graph', 'graph')
    eNode.addLink( 'InsularPole.sulcus_identification', 'sulcus_identification')


    eNode.addLink(  'SulcalinesExtraction.graph','graph' )        
    eNode.addLink(  'SulcalinesExtraction.side','side' )
    eNode.addLink(  'SulcalinesExtraction.sulcus_identification','sulcus_identification' )

    eNode.addLink(  'ParameterizeUnconstrainedHarmonic.side','side')
    eNode.ParameterizeUnconstrainedHarmonic.removeLink('cingular_pole_texture', 'white_mesh')
    eNode.ParameterizeUnconstrainedHarmonic.removeLink('insular_pole_texture', 'white_mesh')
    eNode.ParameterizeUnconstrainedHarmonic.removeLink('white_sulcalines', 'white_mesh')
    eNode.ParameterizeUnconstrainedHarmonic.removeLink('sulcus_labels', 'white_mesh')
    eNode.addLink(  'ParameterizeUnconstrainedHarmonic.cingular_pole_texture', 'CingularPole.pole' )
    eNode.addLink(  'ParameterizeUnconstrainedHarmonic.insular_pole_texture', 'InsularPole.pole' )
    eNode.addLink(  'ParameterizeUnconstrainedHarmonic.cingular_pole_texture', 'CingularPole.pole' )
    eNode.addLink(  'ParameterizeUnconstrainedHarmonic.white_sulcalines','SulcalinesExtraction.white_sulcalines')
    eNode.addLink(  'ParameterizeUnconstrainedHarmonic.sulcus_labels','SulcalinesExtraction.graph_label_basins')

    eNode.addLink(  'HarmonicMappingOrtho.side','side')


    self.setExecutionNode( eNode )
