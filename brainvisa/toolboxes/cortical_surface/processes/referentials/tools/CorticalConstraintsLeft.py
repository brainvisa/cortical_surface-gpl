
from neuroProcesses import *
import shfjGlobals     

name = 'Create Left Hemisphere Cortical Constraints Texture'

userLevel = 2

signature = Signature(
    'Side', Choice("Left"),
    'Lgraph', ReadDiskItem( 'Cortical folds graph', 'Graph',requiredAttributes={ 'side': 'left' } ),
    'sulcus_identification',Choice('name','label'),
    'mri_corrected', ReadDiskItem( 'T1 MRI Bias Corrected', 'Aims readable volume formats' ),
    'left_white_mesh',ReadDiskItem( 'Left Hemisphere White Mesh' , shfjGlobals.aimsMeshFormats),
    'left_white_sulci_mer',WriteDiskItem( 'Left hemisphere longitude constraints texture' ,'Texture',requiredAttributes={ 'side': 'left' } ),
    'left_white_sulci_par',WriteDiskItem( 'Left hemisphere latitude constraints texture' ,'Texture',requiredAttributes={ 'side': 'left' } ),
    'translation',ReadDiskItem('Surface Label Translation','Label Translation' ),
    'ParModel',ReadDiskItem('Latitude Constraint Gyri Model','Gyri Model'),
    'MerModel',ReadDiskItem('Longitude Constraint Gyri Model','Gyri Model'),
    'left_sulci_label_to_sulci_name',WriteDiskItem( 'Sulci To White Texture Translation', 'Text File',requiredAttributes={ 'side': 'left' }),
    'coord',Choice("Both","Longitude","Latitude"),
    'Metric',Choice('Euclidean','Internode')
)

def initialization( self ):
    self.linkParameters( 'left_white_mesh', 'Lgraph' )
    self.linkParameters( 'mri_corrected', 'left_white_mesh' )
    self.linkParameters( 'left_white_sulci_mer', 'left_white_mesh' )
    self.linkParameters( 'left_white_sulci_par', 'left_white_mesh' )
    self.linkParameters( 'left_sulci_label_to_sulci_name', 'Lgraph' )
    self.setOptional('Lgraph', 'left_white_mesh', 'left_white_sulci_mer','left_white_sulci_par' )
    self.sulcus_identification = 'label'
    self.findValue( 'translation', {} )
    self.setOptional('translation')
    self.findValue( 'ParModel', {} )
    self.setOptional('ParModel')
    self.findValue( 'MerModel', {} )
    self.setOptional('MerModel')
    
#     self.translation = os.environ['P4'] + '/shared-main/nomenclature/translation/surfaceReferential.trl'
     
def execution( self, context ):
    Affine_estimation_coef = 0.9
    context.write('Left hemisphere')
    if self.coord in ('Longitude','Both'):
        context.write('processing longitude constraints...')
        context.system('siMeshSulciProjection', '-i', self.left_white_mesh.fullPath() ,  '-g' , self.Lgraph.fullPath(), '-l' , self.translation.fullPath() , '-m', self.MerModel.fullPath() , '-s' , self.sulcus_identification, '-v' , self.mri_corrected.fullPath(), '-o' ,self.left_white_sulci_mer.fullPath() ,'-V' , 1,'-M',2,'-n',5,'-a', Affine_estimation_coef, '-e', 20,'-t', self.left_sulci_label_to_sulci_name.fullPath() ,'-p',1 )
    if self.coord in ('Latitude','Both'):
        context.write('processing latitude constraints...')
        context.system('siMeshSulciProjection', '-i', self.left_white_mesh.fullPath() ,  '-g' , self.Lgraph.fullPath(), '-l' , self.translation.fullPath() , '-m', self.ParModel.fullPath() , '-s' , self.sulcus_identification, '-v' , self.mri_corrected.fullPath(), '-o' ,self.left_white_sulci_par.fullPath() ,'-V' , 1,'-M',2,'-n',5,'-a', Affine_estimation_coef, '-e', 20,'-t', self.left_sulci_label_to_sulci_name.fullPath() ,'-p',1 )
    context.write('Done')
