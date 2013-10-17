
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
from brainvisa.cortical_surface.surface_tools import surface_tools as surfTls
import numpy as np
from soma import aims
#from brainvisa import anatomist


name = 'Insular Pole Projection'

userLevel = 0

#def validation():
#    anatomist.validation()

signature = Signature(
    'graph', ReadDiskItem( 'Labelled Cortical folds graph', 'Graph' ),
    'side', Choice('left', 'right'),
    'mri_corrected', ReadDiskItem( 'T1 MRI Bias Corrected', 'aims readable volume formats' ),
    'sulcus_identification',Choice('label', 'name'),
#    'trl',ReadDiskItem( 'Label Translation' ,'Label Translation'),
    'gyri_model',ReadDiskItem('Gyri Model','Gyri Model' ),
#    'transformation',ReadDiskItem( 'Transform Raw T1 MRI to Talairach-AC/PC-Anatomist', 'Transformation matrix' ),
    'white_mesh',ReadDiskItem( 'Hemisphere White Mesh' , shfjGlobals.aimsMeshFormats),
    'dilation_1', Integer(),
    'erosion', Integer(),
    'dilation_2', Integer(),
    'pole',WriteDiskItem( 'insula pole texture','Texture' )
)

def initialization( self ):
    def linkLabelAtt( self, dummy ):
        if self.graph is not None:
            m = self.graph.get( 'manually_labelled' )
            if m and m == 'Yes':
                return 'name'
        return 'label'
    def linkSide( proc, dummy ):
        if proc.graph is not None:
            return proc.graph.get( 'side' )
    self.linkParameters( 'side', 'graph', linkSide )
    self.linkParameters( 'white_mesh','graph' )
    self.linkParameters( 'pole', 'white_mesh' )
    self.dilation_1 = 2 
    self.erosion = 3 
    self.dilation_2 = 2 

#    self.linkParameters( 'transformation', 'white_mesh' )
    self.linkParameters( 'mri_corrected', 'white_mesh' )
    self.linkParameters( 'sulcus_identification', 'graph', linkLabelAtt )
    self.findValue('gyri_model', { 'sulci_database' : '2008', 'graph_version': '3.0', 'model' : 'gyrus' })
#    try:
#      self.findValue('gyri_model', {})
#    except: pass

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
    tmp_trl = context.temporary('Text File')
    f = open(tmp_trl.fullPath(),'w')
    f.write('%INSULA\n')
    f.close()
    out_trsl_txt = context.temporary('Text File')
    command = ['siMeshSulciProjection','-i',self.white_mesh.fullPath(),'-g',self.graph.fullPath(),'-l',tmp_trl.fullPath(),'-m',self.gyri_model.fullPath(),'-s',self.sulcus_identification,'-v',self.mri_corrected.fullPath(),'-o',self.pole.fullPath(),'-V','1','-M','2','-n','5','-a','0.9','-e','10','-t',out_trsl_txt.fullPath(),'-p','1']
    context.write('Projecting the insula from sulci graph')
    context.system(*command)
   
    re = aims.Reader()
    ws = aims.Writer()

    texture_poles = re.read(self.pole.fullPath())
    tex_S16 = aims.TimeTexture_S16()
    tex_S16[0].assign(texture_poles[0])
#    context.write(max(tex_S16[0].arraydata()))
    ws.write(tex_S16, self.pole.fullPath())
    if self.dilation_1>0:
        context.system('AimsTextureDilation', '-i',self.white_mesh.fullPath(), '-t',self.pole.fullPath(),'-o',self.pole.fullPath(),'-s',self.dilation_1,'--connexity')
    if self.erosion>0:
        context.system('AimsTextureErosion', '-i',self.white_mesh.fullPath(), '-t',self.pole.fullPath(),'-o',self.pole.fullPath(),'-s',self.erosion,'--connexity')
    if self.dilation_2>0:
        context.system('AimsTextureDilation', '-i',self.white_mesh.fullPath(), '-t',self.pole.fullPath(),'-o',self.pole.fullPath(),'-s',self.dilation_2,'--connexity')
    context.write('Dilation Erosion Done')

    mesh = re.read(self.white_mesh.fullPath())
    tex = re.read(self.pole.fullPath())
#    context.write(max(tex[0].arraydata()))
#    context.write(self.side)
    if self.side == 'right':
        tmp_tex_value = 2
    elif self.side == 'left':
        tmp_tex_value = 1
    else:
        context.write('side must be set to left or right!')     
    context.write('Topological correction...')                           
    cingular_tex_clean, cing_tex_boundary = surfTls.textureTopologicalCorrection(mesh, tex[0].arraydata(), tmp_tex_value)
    cingular_tex_clean[np.where(cingular_tex_clean == tmp_tex_value)[0]] = 180
    tex_out = aims.TimeTexture_S16()
    tex_out[0].assign(cingular_tex_clean)
    ws.write(tex_out, self.pole.fullPath())
    context.write('... Done')