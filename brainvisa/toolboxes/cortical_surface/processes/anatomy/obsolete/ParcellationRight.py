
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
import shfjGlobals   

name = 'Right Cortical Surface Parcellation'

userLevel = 2


signature = Signature(
    'Side', Choice("Right"),
    'right_longitude',ReadDiskItem( 'Right hemisphere longitude texture', 'Texture',requiredAttributes={ 'side': 'right' }  ),
    'right_latitude',ReadDiskItem( 'Right hemisphere latitude texture', 'Texture',requiredAttributes={ 'side': 'right' }  ),
    'file_correspondance_constraint',ReadDiskItem( 'Constraint coordinates values', 'Text File' ),
    'right_gyri',WriteDiskItem( 'Right hemisphere gyri parcellation texture','Texture',requiredAttributes={ 'side': 'right' } ),
)

def initialization( self ):
    self.linkParameters( 'right_latitude','right_longitude')
    self.linkParameters( 'right_gyri','right_longitude')
    self.findValue( 'file_correspondance_constraint', {} )
    self.setOptional('file_correspondance_constraint')

def execution( self, context ):
    context.write('Right hemisphere Parcellation')
    context.system('AimsGyriStuff', '-x', self.right_longitude.fullPath() , '-y', self.right_latitude.fullPath() ,  '-a', self.file_correspondance_constraint.fullPath() ,  '-o' , self.right_gyri.fullPath() )
    context.write('Done')