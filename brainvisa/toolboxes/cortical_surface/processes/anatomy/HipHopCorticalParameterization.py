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
from __future__ import absolute_import
from brainvisa.processes import *

name = 'Hip-Hop Cortical Parameterization'
userLevel = 0

signature = Signature( 
  'Lgraph', ReadDiskItem( 'Labelled Cortical folds graph', 'Graph and data',requiredAttributes={ 'side': 'left' } ),
  'Rgraph', ReadDiskItem( 'Labelled Cortical folds graph', 'Graph and data',requiredAttributes={ 'side': 'right' } ),
  'sulcus_identification', Choice('name','label')
  )


def initialization( self ):
    self.linkParameters( 'Rgraph','Lgraph')
    self.sulcus_identification='label'

    eNode = ParallelExecutionNode( self.name, parameterized=self )

    eNode.addChild( 'Hemisphere_Process_Left',
                    ProcessExecutionNode( 'HemisphereProcess2012', optional = 1 ) )
    eNode.addChild( 'Hemisphere_Process_Right',
                    ProcessExecutionNode( 'HemisphereProcess2012', optional = 1 ) )

    eNode.Hemisphere_Process_Left.side='left'
    eNode.Hemisphere_Process_Right.side='right'

    eNode.addDoubleLink( 'Hemisphere_Process_Left.graph', 'Lgraph' )
    eNode.addDoubleLink( 'Hemisphere_Process_Right.graph', 'Rgraph' )

    eNode.addDoubleLink( 'Hemisphere_Process_Left.sulcus_identification', 'sulcus_identification')
    eNode.addDoubleLink( 'Hemisphere_Process_Right.sulcus_identification', 'sulcus_identification')

    self.setExecutionNode( eNode )
