
from neuroProcesses import *
import shfjGlobals   

name = 'Parameterize Right hemisphere'

userLevel = 2


signature = Signature(
    'Side', Choice("Right"),
    'right_white_mesh',ReadDiskItem( 'Right Hemisphere White Mesh' , shfjGlobals.aimsMeshFormats),
    'right_white_sulci_par',ReadDiskItem( 'Right hemisphere latitude cleaned constraints texture', 'Texture',requiredAttributes={ 'side': 'right' }  ),
    'right_white_sulci_mer',ReadDiskItem( 'Right hemisphere longitude cleaned constraints texture', 'Texture',requiredAttributes={ 'side': 'right' }  ),
    'right_cingular_pole',ReadDiskItem( 'Right hippocampus pole texture'  , 'Texture',requiredAttributes={ 'side': 'right' } ),
    'right_poles_texture',ReadDiskItem( 'Right poles texture'  , 'Texture',requiredAttributes={ 'side': 'right' } ),
    'origin_mer', Integer(),
    'stop_criterium', Float(),
    'time_step', Float(),
    'process', Choice("Both","Longitude","Latitude","None"),
    'constraint_attach', Float(),
    'type_beta', Choice("Map","Value"),
    'right_latitude',WriteDiskItem( 'Right hemisphere latitude texture','Texture',requiredAttributes={ 'side': 'right' } ),
    'right_longitude',WriteDiskItem( 'Right hemisphere longitude texture','Texture',requiredAttributes={ 'side': 'right' })
)

def initialization( self ):
#    self.linkParameters( 'right_cingular_pole','right_white_mesh')
    self.linkParameters( 'right_latitude','right_white_mesh')
    self.linkParameters( 'right_longitude','right_white_mesh')
    self.setOptional( 'origin_mer', 'stop_criterium', 'time_step', 'process', 'constraint_attach', 'right_latitude', 'right_longitude', 'right_white_sulci_par', 'right_white_sulci_mer','type_beta' )
    self.stop_criterium = 1e-6
    self.time_step = 0.2
    self.constraint_attach = 0.2
    self.origin_mer = 0

def execution( self, context ):
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

    context.write('Right hemisphere')
    context.system('AimsCorticalReferential', '-i', self.right_white_mesh.fullPath() ,  '-p' , self.right_white_sulci_par.fullPath(), '-m' , self.right_white_sulci_mer.fullPath() ,  '-r', self.right_poles_texture.fullPath(), '-l', self.right_cingular_pole.fullPath() , '-c' , self.stop_criterium, '-d' , self.time_step, '-b' , self.origin_mer , '-f' , process_choice , '-t' , self.constraint_attach  , '-a', t_beta , '-x' , self.right_latitude.fullPath() , '-y' , self.right_longitude.fullPath() )
    context.write('Done')
