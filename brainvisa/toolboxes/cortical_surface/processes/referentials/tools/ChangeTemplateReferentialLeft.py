
from neuroProcesses import *
import shfjGlobals   

name = 'Change Template Referential Left'

userLevel = 2


signature = Signature(
    'Side', Choice("Left"),
    'mri_corrected', ReadDiskItem( 'T1 MRI Bias Corrected', 'Aims readable volume formats' ),
    'transformation_input',ReadDiskItem( 'Transform Raw T1 MRI to Talairach-AC/PC-Anatomist',
    'Transformation matrix' ),
    'left_pole_template',ReadDiskItem( 'Left Cingular Pole Template' , 'Aims readable volume formats' ),
    'talairach_to_subject',WriteDiskItem( 'Talairach To Subject Transformation', 'Transformation matrix' ),
    'subject_to_template',WriteDiskItem( 'Subject To Template Transformation', 'Transformation matrix' ),
    'template_pole_transformation',ReadDiskItem( 'Template Pole To Talairach Tranformation', 'Transformation matrix' ),
    'output_template', WriteDiskItem( 'Left Cingular Pole Template Subject' , 'Aims writable volume formats' ),
)

def initialization( self ):
#    self.linkParameters( 'left_cingular_pole','left_white_mesh')
    self.linkParameters( 'transformation_input','mri_corrected')
    self.linkParameters( 'subject_to_template','mri_corrected')
    self.linkParameters( 'talairach_to_subject','mri_corrected')
    self.linkParameters( 'output_template','mri_corrected')
    self.findValue( 'left_pole_template', {} )
    self.setOptional('left_pole_template')
    self.findValue( 'template_pole_transformation', {} )
    self.setOptional('template_pole_transformation')

def execution( self, context ):
    context.write('Changing Referential...')
    
    param = 0
    context.system('AimsInvertTransformation', '-i', self.transformation_input.fullPath(), '-o', self.talairach_to_subject.fullPath() )
    context.system('AimsComposeTransformation', '-i', self.talairach_to_subject.fullPath(), '-j', self.template_pole_transformation.fullPath(), '-o', self.subject_to_template.fullPath() )
    context.system('VipSplineResamp', '-i', self.left_pole_template.fullPath(), '-o', self.output_template.fullPath(), '-t', self.mri_corrected, '-d', self.subject_to_template.fullPath(), '-ord', param )
    
    context.write('Done')
