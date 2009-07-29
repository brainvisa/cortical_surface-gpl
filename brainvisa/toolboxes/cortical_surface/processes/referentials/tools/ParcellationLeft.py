
from neuroProcesses import *
import shfjGlobals   

name = 'Left Cortical Surface Parcellation'

userLevel = 2


signature = Signature(
    'Side', Choice("Left"),
    'left_longitude',ReadDiskItem( 'Left hemisphere longitude texture', 'Texture',requiredAttributes={ 'side': 'left' }  ),
    'left_latitude',ReadDiskItem( 'Left hemisphere latitude texture', 'Texture',requiredAttributes={ 'side': 'left' }  ),
    'file_correspondance_constraint',ReadDiskItem( 'Constraint coordinates values', 'Text File' ),
    'left_gyri',WriteDiskItem( 'Left hemisphere gyri parcellation texture','Texture',requiredAttributes={ 'side': 'left' } ),
)

def initialization( self ):
    self.linkParameters( 'left_latitude','left_longitude')
    self.linkParameters( 'left_gyri','left_longitude')
    self.findValue( 'file_correspondance_constraint', {} )
    self.setOptional('file_correspondance_constraint')

def execution( self, context ):
    context.write('Left hemisphere Parcellation')
    context.system('AimsGyriStuff', '-x', self.left_longitude.fullPath() , '-y', self.left_latitude.fullPath() ,  '-a', self.file_correspondance_constraint.fullPath() ,  '-o' , self.left_gyri.fullPath() )
    context.write('Done')
