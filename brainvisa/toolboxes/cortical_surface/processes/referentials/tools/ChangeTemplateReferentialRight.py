# -*- coding: utf-8 -*-

#  This software and supporting documentation are distributed by
#      Institut Federatif de Recherche 49
#      CEA/NeuroSpin, Batiment 145,
#      91191 Gif-sur-Yvette cedex
#      France
#
# This software is governed by the CeCILL license version 2 under
# French law and abiding by the rules of distribution of free software.
# You can  use, modify and/or redistribute the software under the 
# terms of the CeCILL license version 2 as circulated by CEA, CNRS
# and INRIA at the following URL "http://www.cecill.info". 
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license version 2 and that you accept its terms.
from brainvisa.processes import *
import brainvisa.tools.aimsGlobals as shfjGlobals
from brainvisa import registration

name = 'Change Template Referential Right'

userLevel = 2


signature = Signature(
    'Side', Choice("Right"),
    'mri_corrected', ReadDiskItem( 'T1 MRI Bias Corrected', 'Aims readable volume formats' ),
    'transformation_input',ReadDiskItem( 'Transform Raw T1 MRI to Talairach-AC/PC-Anatomist',
    'Transformation matrix' ),
    'right_pole_template',ReadDiskItem( 'Right Cingular Pole Template' , 'Aims readable volume formats' ),
    'talairach_to_subject',WriteDiskItem( 'Talairach To Subject Transformation', 'Transformation matrix' ),
    'subject_to_template',WriteDiskItem( 'Subject To Template Transformation', 'Transformation matrix' ),
    'template_pole_transformation',ReadDiskItem( 'Template Pole To Talairach Tranformation', 'Transformation matrix' ),
    'output_template', WriteDiskItem( 'Right Cingular Pole Template Subject' , 'Aims writable volume formats' ),
)

def initialization( self ):
#    self.linkParameters( 'left_cingular_pole','left_white_mesh')
    self.linkParameters( 'transformation_input','mri_corrected')
    self.linkParameters( 'subject_to_template','mri_corrected')
    self.linkParameters( 'talairach_to_subject','mri_corrected')
    self.linkParameters( 'output_template','mri_corrected')
    self.findValue( 'right_pole_template', {} )
    self.setOptional('right_pole_template')
    self.findValue( 'template_pole_transformation', {} )
    self.setOptional('template_pole_transformation')

def execution( self, context ):
    context.write('Changing Referential...')
    
    param = 0
    context.system('AimsInvertTransformation', '-i', self.transformation_input.fullPath(), '-o', self.talairach_to_subject.fullPath() )
    context.system('AimsComposeTransformation', '-i', self.talairach_to_subject.fullPath(), '-j', self.template_pole_transformation.fullPath(), '-o', self.subject_to_template.fullPath() )
    context.system('VipSplineResamp', '-i', self.right_pole_template.fullPath(), '-o', self.output_template.fullPath(), '-t', self.mri_corrected, '-d', self.subject_to_template.fullPath(), '-ord', param )
    tm = registration.getTransformationManager()
    tm.copyReferential( self.mri_corrected, self.output_template )

    context.write('Done')
