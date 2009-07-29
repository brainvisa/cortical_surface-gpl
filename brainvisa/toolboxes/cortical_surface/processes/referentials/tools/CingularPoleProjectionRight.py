
from neuroProcesses import *
import shfjGlobals
from brainvisa import anatomist

name = 'Right Cingular Pole Projection'

userLevel = 2

def validation():
    anatomist.validation()

signature = Signature(
    'Side', Choice("Right"),
    'right_white_mesh',ReadDiskItem( 'Right Hemisphere White Mesh' , shfjGlobals.aimsMeshFormats),
    'right_pole_template',ReadDiskItem( 'Right Cingular Pole Template Subject' , 'Aims readable volume formats' ),
    'right_pole',WriteDiskItem( 'Right hippocampus pole texture','Texture',requiredAttributes={ 'side': 'right' } )
)

def initialization( self ):
    self.linkParameters( 'right_pole', 'right_white_mesh' )
    #self.findValue( 'right_pole_template', {} )
    #self.setOptional('right_pole_template')
     
def execution( self, context ):
    a = anatomist.Anatomist()
    mesh = a.loadObject( self.right_white_mesh.fullPath() )
    vol = a.loadObject( self.right_pole_template.fullPath() )
    fusion = a.fusionObjects( [mesh, vol], method='Fusion3DMethod' )
    fusion.exportTexture(filename=self.right_pole.fullPath())

