
from neuroProcesses import *
import shfjGlobals     

name = 'Create Right Hemisphere Cortical Constraints Texture'

userLevel = 2

signature = Signature(
    'Side', Choice("Right"),
    'Rgraph', ReadDiskItem( 'Cortical folds graph', 'Graph',requiredAttributes={ 'side': 'right' } ),
    'sulcus_identification',Choice('name','label'),
    'mri_corrected', ReadDiskItem( 'T1 MRI Bias Corrected', 'Aims readable volume formats' ),
    'right_white_mesh',ReadDiskItem( 'Right Hemisphere White Mesh' , shfjGlobals.aimsMeshFormats),
    'right_white_sulci_mer',WriteDiskItem( 'Right hemisphere longitude constraints texture' ,'Texture',requiredAttributes={ 'side': 'right' } ),
    'right_white_sulci_par',WriteDiskItem( 'Right hemisphere latitude constraints texture' ,'Texture',requiredAttributes={ 'side': 'right' } ),
    'translation',ReadDiskItem('Surface Label Translation','Label Translation' ),
    'ParModel',ReadDiskItem('Latitude Constraint Gyri Model','Gyri Model'),
    'MerModel',ReadDiskItem('Longitude Constraint Gyri Model','Gyri Model'),
    'right_sulci_label_to_sulci_name',WriteDiskItem( 'Sulci To White Texture Translation', 'Text File',requiredAttributes={ 'side': 'right' }),
    'coord',Choice("Both","Longitude","Latitude"),
    'Metric',Choice('Euclidean','Internode')
)

def initialization( self ):
    self.linkParameters( 'right_white_mesh', 'Rgraph' )
    self.linkParameters( 'mri_corrected', 'right_white_mesh' )
    self.linkParameters( 'right_white_sulci_mer', 'right_white_mesh' )
    self.linkParameters( 'right_white_sulci_par', 'right_white_mesh' )
    self.linkParameters( 'right_sulci_label_to_sulci_name', 'Rgraph' )
    self.setOptional('Rgraph', 'right_white_mesh', 'right_white_sulci_mer', 'right_white_sulci_par'  )
    self.sulcus_identification = 'label'
    self.findValue( 'translation', {} )
    self.setOptional('translation')
    self.findValue( 'ParModel', {} )
    self.setOptional('ParModel')
    self.findValue( 'MerModel', {} )
    self.setOptional('MerModel')
     
def execution( self, context ):
    Affine_estimation_coef = 0.9
    context.write('Right hemisphere')
    if self.coord in ('Longitude','Both'):
        context.write('processing longitude constraints...')
        context.system('siMeshSulciProjection', '-i', self.right_white_mesh.fullPath() ,  '-g' , self.Rgraph.fullPath(), '-l' , self.translation.fullPath() , '-m', self.MerModel.fullPath() , '-s' , self.sulcus_identification, '-v' , self.mri_corrected.fullPath(), '-o' ,self.right_white_sulci_mer.fullPath() ,'-V', 1,'-M',2,'-n',5,'-a', Affine_estimation_coef, '-e', 20,'-t',self.right_sulci_label_to_sulci_name.fullPath() ,'-p',1 )
    if self.coord in ('Latitude','Both'):
        context.write('processing latitude constraints...')
        context.system('siMeshSulciProjection', '-i', self.right_white_mesh.fullPath() ,  '-g' , self.Rgraph.fullPath(), '-l' , self.translation.fullPath() , '-m', self.ParModel.fullPath() , '-s' , self.sulcus_identification, '-v' , self.mri_corrected.fullPath(), '-o' ,self.right_white_sulci_par.fullPath() ,'-V', 1,'-M',2,'-n',5,'-a', Affine_estimation_coef, '-e', 20,'-t',self.right_sulci_label_to_sulci_name.fullPath() ,'-p',1 )
    context.write('Done')
