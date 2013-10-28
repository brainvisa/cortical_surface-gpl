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
userLevel = 0
 
signature = Signature(
  'side', Choice('left', 'right'),
  'graph', ReadDiskItem( 'Labelled Cortical folds graph', 'Graph and data'  ),
  'sulcus_identification', Choice('name','label')
)

def initialization( self ):
    eNode = SerialExecutionNode( self.name, parameterized=self )

#    eNode.addChild( 'ChangeTemplateReferential',
#                    ProcessExecutionNode( 'ChangeTemplateReferential', optional = 1 ) )
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

    eNode.addChild( 'CoordinatesFromHipHopMapping',
                    ProcessExecutionNode( 'CoordinatesFromHipHopMapping', optional=1) )
                    
    eNode.addChild( 'ParcelsTextureFromCoordinates',
                    ProcessExecutionNode( 'ParcelsTextureFromCoordinates', optional=1) )

    eNode.addDoubleLink( 'CingularPole.side', 'side' )

    eNode.InsularPole.removeLink( 'side', 'graph' )
    eNode.addDoubleLink( 'InsularPole.white_mesh', 'CingularPole.white_mesh' )
    eNode.addDoubleLink( 'InsularPole.side', 'side')
    eNode.addDoubleLink( 'InsularPole.graph', 'graph')
    eNode.addDoubleLink( 'InsularPole.sulcus_identification', 'sulcus_identification')

    eNode.SulcalinesExtraction.removeLink( 'white_mesh', 'graph' )
    eNode.SulcalinesExtraction.removeLink( 'mri', 'white_mesh' )
    eNode.addDoubleLink(  'SulcalinesExtraction.graph', 'graph' )
    eNode.addDoubleLink(  'SulcalinesExtraction.side', 'side' )
    eNode.addDoubleLink(  'SulcalinesExtraction.sulcus_identification', 'sulcus_identification' )
    eNode.addDoubleLink( 'SulcalinesExtraction.white_mesh', 'CingularPole.white_mesh' )
    eNode.addDoubleLink( 'SulcalinesExtraction.mri', 'InsularPole.mri_corrected' )

    eNode.ParameterizeUnconstrainedHarmonic.removeLink( 'cingular_pole_texture', 'white_mesh' )
    eNode.ParameterizeUnconstrainedHarmonic.removeLink('insular_pole_texture', 'white_mesh')
    eNode.ParameterizeUnconstrainedHarmonic.removeLink('white_sulcalines', 'white_mesh')
    eNode.ParameterizeUnconstrainedHarmonic.removeLink('sulcus_labels', 'white_mesh')
    eNode.addDoubleLink( 'ParameterizeUnconstrainedHarmonic.white_mesh', 'CingularPole.white_mesh' )
    eNode.addDoubleLink( 'ParameterizeUnconstrainedHarmonic.side', 'side' )
    eNode.addDoubleLink( 'ParameterizeUnconstrainedHarmonic.cingular_pole_texture', 'CingularPole.pole' )
    eNode.addDoubleLink(  'ParameterizeUnconstrainedHarmonic.insular_pole_texture', 'InsularPole.pole' )
    eNode.addDoubleLink( 'ParameterizeUnconstrainedHarmonic.white_sulcalines', 'SulcalinesExtraction.white_sulcalines')
    eNode.addDoubleLink( 'ParameterizeUnconstrainedHarmonic.sulcus_labels', 'SulcalinesExtraction.graph_label_basins')

    eNode.HarmonicMappingOrtho.removeLink( 'side', 'rectangular_mesh' )
    eNode.HarmonicMappingOrtho.removeLink( 'boundary_texture', 'rectangular_mesh' )
    eNode.HarmonicMappingOrtho.removeLink( 'corresp_indices_texture', 'rectangular_mesh' )
    eNode.HarmonicMappingOrtho.removeLink( 'white_sulcalines', 'rectangular_mesh' )
    eNode.HarmonicMappingOrtho.removeLink( 'sulcus_labels', 'rectangular_mesh' )
    eNode.addDoubleLink( 'HarmonicMappingOrtho.rectangular_mesh', 'ParameterizeUnconstrainedHarmonic.rectangular_mesh' )
    eNode.addDoubleLink(  'HarmonicMappingOrtho.side','side')
    eNode.addDoubleLink( 'HarmonicMappingOrtho.boundary_texture', 'ParameterizeUnconstrainedHarmonic.boundary_texture' )
    eNode.addDoubleLink( 'HarmonicMappingOrtho.corresp_indices_texture', 'ParameterizeUnconstrainedHarmonic.corresp_indices_texture' )
    eNode.addDoubleLink( 'HarmonicMappingOrtho.white_sulcalines', 'ParameterizeUnconstrainedHarmonic.rectangular_white_sulcalines' )
    eNode.addDoubleLink( 'HarmonicMappingOrtho.sulcus_labels', 'SulcalinesExtraction.graph_label_basins' )

    eNode.CoordinatesFromHipHopMapping.removeLink( 'boundary_texture', 'cstr_rectangular_mesh' )
    eNode.CoordinatesFromHipHopMapping.removeLink( 'corresp_indices_texture', 'cstr_rectangular_mesh' )
    eNode.CoordinatesFromHipHopMapping.removeLink( 'white_mesh_parts', 'cstr_rectangular_mesh' )
    #eNode.CoordinatesFromHipHopMapping.removeLink( 'model_file', 'cstr_rectangular_mesh' )
    eNode.addDoubleLink( 'CoordinatesFromHipHopMapping.cstr_rectangular_mesh', 'HarmonicMappingOrtho.cstr_rectangular_mesh' )
    eNode.addDoubleLink( 'CoordinatesFromHipHopMapping.boundary_texture', 'ParameterizeUnconstrainedHarmonic.boundary_texture' )
    eNode.addDoubleLink( 'CoordinatesFromHipHopMapping.corresp_indices_texture', 'ParameterizeUnconstrainedHarmonic.corresp_indices_texture' )
    eNode.addDoubleLink( 'CoordinatesFromHipHopMapping.white_mesh_parts', 'ParameterizeUnconstrainedHarmonic.white_mesh_parts' )
    eNode.addDoubleLink( 'CoordinatesFromHipHopMapping.model_file', 'HarmonicMappingOrtho.model_file')

    eNode.ParcelsTextureFromCoordinates.removeLink( 'side', 'latitude' )
    eNode.ParcelsTextureFromCoordinates.removeLink( 'longitude', 'latitude' )
    eNode.ParcelsTextureFromCoordinates.removeLink( 'model_file', 'latitude' )
    eNode.addDoubleLink( 'ParcelsTextureFromCoordinates.latitude', 'CoordinatesFromHipHopMapping.latitude')
    eNode.addDoubleLink( 'ParcelsTextureFromCoordinates.side', 'side' )
    eNode.addDoubleLink( 'ParcelsTextureFromCoordinates.longitude', 'CoordinatesFromHipHopMapping.longitude')
    eNode.addDoubleLink( 'ParcelsTextureFromCoordinates.model_file', 'HarmonicMappingOrtho.model_file' )

    self.setExecutionNode( eNode )
