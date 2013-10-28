
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
# from brainvisa import anatomist
from soma import aims
from brainvisa.cortical_surface.surface_tools import surface_tools as surfTls
import numpy as np

name = 'Cingular Pole Projection'

userLevel = 0

# def validation():
#     anatomist.validation()

signature = Signature(
    'white_mesh',ReadDiskItem( 'Hemisphere White Mesh' , shfjGlobals.aimsMeshFormats),
    'side', Choice('left', 'right'),
#    'pole_template',ReadDiskItem( 'Cingular Pole Template Subject' , 'Aims readable volume formats' ),
    'subject_transformation',ReadDiskItem( 'Transform Raw T1 MRI to Talairach-AC/PC-Anatomist',
    'Transformation matrix' ),
    'pole_template',ReadDiskItem( 'Cingular Pole Template', 'Aims readable volume formats' ),
    'template_pole_transformation',ReadDiskItem( 'Template Pole To Talairach Tranformation', 'Transformation matrix' ),
    'dilation_1', Integer(),
    'erosion', Integer(),
    'dilation_2', Integer(),
    'pole',WriteDiskItem( 'Hippocampus pole texture', 'aims Texture formats' )
)

def initialization( self ):
    def linkSide( proc, dummy ):
        if proc.white_mesh is not None:
            return proc.white_mesh.get( 'side' )
    self.linkParameters( 'side', 'white_mesh', linkSide )
    self.linkParameters( 'pole', 'white_mesh' )
    #self.findValue( 'pole_template', { 'side' : self.side , 'filename_variable' : 'Template_icbm', '_ontology' : 'shared'} )
    self.linkParameters( 'pole_template', 'white_mesh' )
    self.findValue( 'template_pole_transformation', {} )
    self.linkParameters( 'subject_transformation','white_mesh')
    self.dilation_1 = 2#7
    self.erosion = 6#11
    self.dilation_2 = 4#2 
 
def execution( self, context ):
    context.write('Changing Referential...')
    cingular_tex_value = 1
    tmp_trm1 = context.temporary(  'Transformation matrix' )
    tmp_trm2 = context.temporary(  'Transformation matrix' )
    tmp_mesh = context.temporary(  'Mesh mesh' )
    context.system('AimsInvertTransformation', '-i', self.template_pole_transformation, '-o', tmp_trm1.fullPath() )
    context.system('AimsComposeTransformation', '-i', tmp_trm1.fullPath(), '-j', self.subject_transformation, '-o', tmp_trm2.fullPath() )
    context.system('AimsMeshTransform', '-i',self.white_mesh.fullPath(),'-o',tmp_mesh.fullPath(),'-t',tmp_trm2.fullPath())
    command = [ 'AimsCreateTexture', '-m',tmp_mesh.fullPath(), '-v',self.pole_template.fullPath(), '-t', self.pole.fullPath() ]

#    command = [ 'AimsCreateTexture', '-m',self.white_mesh.fullPath(), '-v',self.pole_template.fullPath(), '-t', self.pole.fullPath() ]
    context.system(*command)
    context.write('Projection Done')

    re = aims.Reader()
    ws = aims.Writer()

    texture_poles = re.read(self.pole.fullPath())
    atex = np.zeros(texture_poles[0].arraydata().shape)
    atex[texture_poles[0].arraydata() > 0] = cingular_tex_value
    tex_S16 = aims.TimeTexture_S16()
    tex_S16[0].assign(atex)
    ws.write(tex_S16, self.pole.fullPath())
    if self.dilation_1>0:
        context.system('AimsTextureDilation', '-i',self.white_mesh.fullPath(), '-t',self.pole.fullPath(),'-o',self.pole.fullPath(),'-s',self.dilation_1,'--connexity')#10
    if self.erosion>0:
        context.system('AimsTextureErosion', '-i',self.white_mesh.fullPath(), '-t',self.pole.fullPath(),'-o',self.pole.fullPath(),'-s',self.erosion,'--connexity')#10
    if self.dilation_2>0:
        context.system('AimsTextureDilation', '-i',self.white_mesh.fullPath(), '-t',self.pole.fullPath(),'-o',self.pole.fullPath(),'-s',self.dilation_1,'--connexity')#10
    context.write('Dilation Erosion Done')
    context.write('Topological correction...')
    mesh = re.read(self.white_mesh.fullPath())
    tex = re.read(self.pole.fullPath())
    cingular_tex_clean, cing_tex_boundary = surfTls.textureTopologicalCorrection(mesh, tex[0].arraydata(), cingular_tex_value)
    tex_out = aims.TimeTexture_S16()
    tex_out[0].assign(cingular_tex_clean)
    ws.write(tex_out, self.pole.fullPath())
   
#    ws.write(tex, self.pole.fullPath())
    context.write('... Done')