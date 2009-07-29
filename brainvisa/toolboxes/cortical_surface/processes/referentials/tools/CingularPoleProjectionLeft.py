
from neuroProcesses import *
import shfjGlobals
from brainvisa import anatomist

name = 'Left Cingular Pole Projection'

userLevel = 2

def validation():
    anatomist.validation()

signature = Signature(
    'Side', Choice("Left"),
    'left_white_mesh',ReadDiskItem( 'Left Hemisphere White Mesh' , shfjGlobals.aimsMeshFormats),
    'left_pole_template',ReadDiskItem( 'Left Cingular Pole Template Subject' , 'Aims readable volume formats' ),
    'left_pole',WriteDiskItem( 'Left hippocampus pole texture','Texture',requiredAttributes={ 'side': 'left' } )
)

def initialization( self ):
    self.linkParameters( 'left_pole', 'left_white_mesh' )
    #self.findValue( 'left_pole_template', {} )
    #self.setOptional('left_pole_template')
     
def execution( self, context ):
    a = anatomist.Anatomist()
    mesh = a.loadObject( self.left_white_mesh.fullPath() )
    vol = a.loadObject( self.left_pole_template.fullPath() )
    fusion = a.fusionObjects( [mesh, vol], method='Fusion3DMethod' )
    fusion.exportTexture(filename=self.left_pole.fullPath())
