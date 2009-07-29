
from neuroProcesses import *
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
