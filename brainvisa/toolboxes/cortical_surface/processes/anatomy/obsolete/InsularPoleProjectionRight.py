
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
import shfjGlobals
from brainvisa.cortical_surface.surface_tools import texture_tools as textureTls
import numpy as np
#from brainvisa import anatomist


name = 'Right Insular Pole Projection'

userLevel = 2

#def validation():
#    anatomist.validation()

signature = Signature(
    'Side', Choice("Right"),
    'Rgraph', ReadDiskItem( 'Cortical folds graph', 'Graph',requiredAttributes={ 'side': 'right' } ),
    'mri_corrected', ReadDiskItem( 'T1 MRI Bias Corrected', 'aims readable volume formats' ),
    'sulcus_identification',Choice('name','label'),
#    'trl',ReadDiskItem( 'Label Translation' ,'Label Translation'),
    'gyri_model',ReadDiskItem('Gyri Model','Gyri Model' ),
    'transformation', ReadDiskItem('Transformation matrix', 'Transformation matrix'),
    'right_white_mesh',ReadDiskItem( 'Right Hemisphere White Mesh' , shfjGlobals.aimsMeshFormats),
    'right_pole',WriteDiskItem( 'Right insula pole texture','Texture',requiredAttributes={ 'side': 'right' } )
)

def initialization( self ):
    def linkLabelAtt( self, dummy ):
        if self.Rgraph is not None:
            m = self.Rgraph.get( 'manually_labelled' )
            if m and m == 'Yes':
                return 'name'
        return 'label'
    self.linkParameters( 'right_white_mesh','Rgraph' )
    self.linkParameters( 'right_pole', 'right_white_mesh' )
    self.linkParameters( 'mri_corrected', 'right_white_mesh' )
    self.linkParameters( 'sulcus_identification', 'Rgraph', linkLabelAtt )
    try:
      self.gyri_model = databases.getDiskItemFromUuid( '172c4168-a9d3-dc41-464c-1226ad07c19c' )
    except: pass

    #self.findValue( 'right_pole_template', {} )
    #self.setOptional('right_pole_template')
     
def execution( self, context ):
    context.write('TODO project the insula spam')
#    tmp_white = '/tmp/tmp_white.gii'#context.temporary( 'GIFTI file' )
#    context.system('AimsMeshTransform', '-i',self.right_white_mesh.fullPath(),'-o',tmp_white,'-t',self.transformation)
#    context.write('Transormation to Talairach Done')
#    spam_INSULA='/home/toz/BVs/brainvisa-4.3.0/share/brainvisa-share-4.3/models/models_2008/descriptive_models/segments/talairach_spam_right/spam_distribs/bayesian_spam_density_INSULA_right.nii.gz'
#    command = [ 'AimsCreateTexture', '-m',tmp_white, '-v',spam_INSULA, '-t', self.right_pole.fullPath() ]
#    context.write(*command)
#    context.system(*command)
#    context.write('Projection Done')
    tmp_trl = context.temporary(  'GIS image' )
    f = open(tmp_trl.fullPath(),'w')
    f.write('%INSULA\n')
    f.close()
    out_trsl_txt = context.temporary('Text File')
    command = ['siMeshSulciProjection','-i',self.right_white_mesh.fullPath(),'-g',self.Rgraph.fullPath(),'-l',tmp_trl.fullPath(),'-m',self.gyri_model.fullPath(),'-s',self.sulcus_identification,'-v',self.mri_corrected.fullPath(),'-o',self.right_pole.fullPath(),'-V','1','-M','2','-n','5','-a','0.9','-e','10','-t',out_trsl_txt.fullPath(),'-p','1']
    context.system(*command)
   
    from soma import aims
    re = aims.Reader()
    ws = aims.Writer()

    texture_poles = re.read(self.right_pole.fullPath())
    tex_S16 = aims.TimeTexture_S16()
    tex_S16[0].assign(texture_poles[0])
    context.write(max(tex_S16[0].arraydata()))
    ws.write(tex_S16, self.right_pole.fullPath())
    context.system('AimsTextureDilation', '-i',self.right_white_mesh.fullPath(), '-t',self.right_pole.fullPath(),'-o',self.right_pole.fullPath(),'-s','2','--connexity')
    context.system('AimsTextureErosion', '-i',self.right_white_mesh.fullPath(), '-t',self.right_pole.fullPath(),'-o',self.right_pole.fullPath(),'-s','3','--connexity')
    context.system('AimsTextureDilation', '-i',self.right_white_mesh.fullPath(), '-t',self.right_pole.fullPath(),'-o',self.right_pole.fullPath(),'-s','2','--connexity')
    context.write('Dilation Erosion Done')

    mesh = re.read(self.right_white_mesh.fullPath())
    tex = re.read(self.right_pole.fullPath())
    context.write(max(tex[0].arraydata()))
    tmp_tex_value = 2
    cingular_tex_clean, cing_tex_boundary = textureTls.poleTextureClean(mesh, tex[0].arraydata(), tmp_tex_value)
    cingular_tex_clean[np.where(cingular_tex_clean == tmp_tex_value)[0]] = 180
    tex_out = aims.TimeTexture_S16()
    tex_out[0].assign(cingular_tex_clean)
    ws.write(tex_out, self.right_pole.fullPath())
    context.write('cleaning Done')
