
from neuroProcesses import *
import shfjGlobals   

name = 'Constraint Cleaner Right hemisphere'

userLevel = 2


signature = Signature(
    'Side', Choice("Right"),
    'right_white_mesh',ReadDiskItem( 'Right Hemisphere White Mesh' , shfjGlobals.aimsMeshFormats),
    'right_cingular_pole',ReadDiskItem( 'Right hippocampus pole texture'  , 'Texture',requiredAttributes={ 'side': 'right' } ),
    'right_white_sulci_mer',ReadDiskItem( 'Right hemisphere longitude constraints texture', 'Texture',requiredAttributes={ 'side': 'right' }  ),
    'right_white_sulci_par',ReadDiskItem( 'Right hemisphere latitude constraints texture', 'Texture',requiredAttributes={ 'side': 'right' }  ),
    'file_correspondance_constraint',ReadDiskItem( 'Constraint coordinates values', 'Text File' ),
    'right_sulci_label_to_sulci_name',ReadDiskItem( 'Sulci To White Texture Translation', 'Text File' ),
    'constraint_dist_param', Float(),
    'curvature_param', Float(),
    'elasticity_param', Float(),
    'right_poles_texture',WriteDiskItem( 'Right poles texture'  , 'Texture',requiredAttributes={ 'side': 'right' } ),
    'right_white_sulci_mer_cleaned',WriteDiskItem( 'Right hemisphere longitude cleaned constraints texture', 'Texture',requiredAttributes={ 'side': 'right' }  ),
    'right_white_sulci_par_cleaned',WriteDiskItem( 'Right hemisphere latitude cleaned constraints texture', 'Texture',requiredAttributes={ 'side': 'right' }  ),
    'process', Choice("Both","Longitude","Latitude")
)

def initialization( self ):
    #self.linkParameters( 'right_cingular_pole','right_white_mesh')
    #self.linkParameters( 'right_white_sulci_mer','right_cingular_pole')
    #self.linkParameters( 'right_white_sulci_par','right_cingular_pole')
    #self.linkParameters( 'right_cingular_pole','right_white_sulci_mer')
    #self.linkParameters( 'right_cingular_pole','right_white_sulci_par')
    #self.linkParameters( 'right_white_sulci_par','right_white_sulci_mer')
    #self.linkParameters( 'right_white_sulci_mer','right_white_sulci_par')
    self.linkParameters( 'right_white_sulci_mer_cleaned','right_white_mesh')
    self.linkParameters( 'right_white_sulci_par_cleaned','right_white_mesh')
    self.linkParameters( 'right_poles_texture','right_white_mesh')
    #self.linkParameters( 'right_sulci_label_to_sulci_name', 'right_white_mesh' )
    self.setOptional( 'process', 'right_white_sulci_mer_cleaned', 'right_white_sulci_par_cleaned', 'constraint_dist_param', 'curvature_param', 'elasticity_param'  )
    self.findValue( 'file_correspondance_constraint', {} )
    self.setOptional('file_correspondance_constraint')
    
    self.constraint_dist_param = 20
    self.curvature_param = 500
    self.elasticity_param = 800
    
def execution( self, context ):
    if self.process == 'Longitude' :
        process_choice = 0
    else :
        if self.process == 'Latitude' :
            process_choice = 1
        else :
            process_choice = 2

    side = '_right'
    
    context.write('Right hemisphere')
    context.system('AimsConstraintCleaner', '-m', self.right_white_mesh.fullPath() ,  '-t', self.right_cingular_pole.fullPath() ,  '-p' , self.right_poles_texture.fullPath(), '-x' , self.right_white_sulci_mer.fullPath() , '-y', self.right_white_sulci_par.fullPath() , '-a' , self.right_white_sulci_mer_cleaned.fullPath(), '-b', self.right_white_sulci_par_cleaned.fullPath(), '-f' , self.file_correspondance_constraint.fullPath(), '-g', self.right_sulci_label_to_sulci_name.fullPath(), '-c',  process_choice, '-s', side, '-i', self.constraint_dist_param, '-j', self.curvature_param , '-k' , self.elasticity_param )
    context.write('Done')

