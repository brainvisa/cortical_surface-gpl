
from neuroProcesses import *
import shfjGlobals   

name = 'Parameterize Left hemisphere'

userLevel = 2


signature = Signature(
    'Side', Choice("Left"),
    'left_white_mesh',ReadDiskItem( 'Left Hemisphere White Mesh' , shfjGlobals.aimsMeshFormats),
    'left_white_sulci_par',ReadDiskItem( 'Left hemisphere latitude cleaned constraints texture', 'Texture',requiredAttributes={ 'side': 'left' }  ),
    'left_white_sulci_mer',ReadDiskItem( 'Left hemisphere longitude cleaned constraints texture', 'Texture',requiredAttributes={ 'side': 'left' }  ),
    'left_cingular_pole',ReadDiskItem( 'Left hippocampus pole texture'  , 'Texture',requiredAttributes={ 'side': 'left' } ),
    'left_poles_texture',ReadDiskItem( 'Left poles texture'  , 'Texture',requiredAttributes={ 'side': 'left' } ),
    'origin_mer', Integer(),
    'stop_criterium', Float(),
    'time_step', Float(),
    'process', Choice("Both","Longitude","Latitude","None"),
    'constraint_attach', Float(),
    'type_beta', Choice("Map","Value"),
    'left_latitude',WriteDiskItem( 'Left hemisphere latitude texture','Texture',requiredAttributes={ 'side': 'left' } ),
    'left_longitude',WriteDiskItem( 'Left hemisphere longitude texture','Texture',requiredAttributes={ 'side': 'left' })
)

def initialization( self ):
#    self.linkParameters( 'left_cingular_pole','left_white_mesh')
    self.linkParameters( 'left_latitude','left_white_mesh')
    self.linkParameters( 'left_longitude','left_white_mesh')
    self.setOptional( 'origin_mer', 'stop_criterium', 'time_step', 'process', 'constraint_attach', 'left_latitude', 'left_longitude', 'left_white_sulci_par', 'left_white_sulci_mer','type_beta' )
    self.stop_criterium = 1e-6
    self.time_step = 0.2
    self.constraint_attach = 0.2
    self.origin_mer = 0

def execution( self, context ):
    context.write(self.Side)
    if self.process == 'Longitude' :
        process_choice = 2
    else :
        if self.process == 'Latitude' :
            process_choice = 1
        else :
            if self.process == 'None' :
                process_choice = 3
            else :
                process_choice = 0

    if self.type_beta == 'Value' :
        t_beta = 0
    else :
        t_beta = 1

    context.write('Left hemisphere')
    context.system('AimsCorticalReferential', '-i', self.left_white_mesh.fullPath() ,  '-p' , self.left_white_sulci_par.fullPath(), '-m' , self.left_white_sulci_mer.fullPath() , '-r', self.left_poles_texture.fullPath() ,'-l', self.left_cingular_pole.fullPath() , '-c' , self.stop_criterium, '-d' , self.time_step, '-b' , self.origin_mer , '-f' , process_choice , '-t' , self.constraint_attach , '-a', t_beta , '-x' , self.left_latitude.fullPath() , '-y' , self.left_longitude.fullPath() )
    context.write('Done')
