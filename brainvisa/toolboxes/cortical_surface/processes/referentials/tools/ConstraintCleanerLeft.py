
from neuroProcesses import *
import shfjGlobals   

name = 'Constraint Cleaner Left hemisphere'

userLevel = 2


signature = Signature(
    'Side', Choice("Left"),
    'left_white_mesh',ReadDiskItem( 'Left Hemisphere White Mesh' , shfjGlobals.aimsMeshFormats),
    'left_cingular_pole',ReadDiskItem( 'Left hippocampus pole texture'  , 'Texture',requiredAttributes={ 'side': 'left' } ),
    'left_white_sulci_mer',ReadDiskItem( 'Left hemisphere longitude constraints texture', 'Texture',requiredAttributes={ 'side': 'left' }  ),
    'left_white_sulci_par',ReadDiskItem( 'Left hemisphere latitude constraints texture', 'Texture',requiredAttributes={ 'side': 'left' }  ),
    'file_correspondance_constraint',ReadDiskItem( 'Constraint coordinates values', 'Text File' ),
    'left_sulci_label_to_sulci_name',ReadDiskItem( 'Sulci To White Texture Translation', 'Text File' ),
    'constraint_dist_param', Float(),
    'curvature_param', Float(),
    'elasticity_param', Float(),
    'left_poles_texture',WriteDiskItem( 'Left poles texture'  , 'Texture',requiredAttributes={ 'side': 'left' } ),
    'left_white_sulci_mer_cleaned',WriteDiskItem( 'Left hemisphere longitude cleaned constraints texture', 'Texture',requiredAttributes={ 'side': 'left' }  ),
    'left_white_sulci_par_cleaned',WriteDiskItem( 'Left hemisphere latitude cleaned constraints texture', 'Texture',requiredAttributes={ 'side': 'left' }  ),
    'process', Choice("Both","Longitude","Latitude" )
)

def initialization( self ):
    #self.linkParameters( 'left_cingular_pole','left_white_mesh')
    #self.linkParameters( 'left_white_sulci_mer','left_white_mesh')
    #self.linkParameters( 'left_white_sulci_par','left_white_mesh')
    #self.linkParameters( 'left_cingular_pole','left_white_sulci_mer')
    #self.linkParameters( 'left_cingular_pole','left_white_sulci_par')
    #self.linkParameters( 'left_white_sulci_par','left_white_sulci_mer')
    #self.linkParameters( 'left_white_sulci_mer','left_white_sulci_par')
    self.linkParameters( 'left_white_sulci_mer_cleaned','left_white_mesh')
    self.linkParameters( 'left_white_sulci_par_cleaned','left_white_mesh')
    self.linkParameters( 'left_poles_texture','left_white_mesh')
    #self.linkParameters( 'left_sulci_label_to_sulci_name', 'left_white_mesh' )
    self.setOptional( 'process', 'left_white_sulci_mer_cleaned', 'left_white_sulci_par_cleaned', 'constraint_dist_param', 'curvature_param', 'elasticity_param'  )
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

    side = '_left'
    
    
    context.write('Left hemisphere')
    context.system('AimsConstraintCleaner', '-m', self.left_white_mesh.fullPath() ,  '-t', self.left_cingular_pole.fullPath() ,  '-p' , self.left_poles_texture.fullPath(), '-x' , self.left_white_sulci_mer.fullPath() , '-y', self.left_white_sulci_par.fullPath() , '-a' , self.left_white_sulci_mer_cleaned.fullPath(), '-b', self.left_white_sulci_par_cleaned.fullPath(), '-f' , self.file_correspondance_constraint.fullPath(), '-g', self.left_sulci_label_to_sulci_name.fullPath(), '-c',  process_choice, '-s', side, '-i', self.constraint_dist_param, '-j', self.curvature_param , '-k' , self.elasticity_param )
    context.write('Done')
