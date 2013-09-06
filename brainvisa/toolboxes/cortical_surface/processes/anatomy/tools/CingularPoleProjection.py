
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

userLevel = 2

# def validation():
#     anatomist.validation()

signature = Signature(
    'side', Choice('left', 'right'),
    'white_mesh',ReadDiskItem( 'Hemisphere White Mesh' , shfjGlobals.aimsMeshFormats),
    'pole_template',ReadDiskItem( 'Cingular Pole Template Subject' , 'Aims readable volume formats' ),
    'dilation_1', Integer(),
    'erosion', Integer(),
    'dilation_2', Integer(),
    'pole',WriteDiskItem( 'Hippocampus pole texture','Texture' )
)

def initialization( self ):
    self.linkParameters( 'pole', 'white_mesh' )
    self.linkParameters( 'pole_template', 'white_mesh' )
    self.dilation_1 = 7
    self.erosion = 11
    self.dilation_2 = 2 
 
def execution( self, context ):
#     a = anatomist.Anatomist()
#     mesh = a.loadObject( self.white_mesh.fullPath() )
#     vol = a.loadObject( self.pole_template.fullPath() )
#     fusion = a.fusionObjects( [mesh, vol], method='Fusion3DMethod' )
#     fusion.exportTexture(filename=self.pole.fullPath())
#    context.system('AimsMeshTransform', '-i',self.left_white_mesh.fullPath(),'-o',tmp_white,'-t',self.transformation)
    command = [ 'AimsCreateTexture', '-m',self.white_mesh.fullPath(), '-v',self.pole_template.fullPath(), '-t', self.pole.fullPath() ]
#    context.write(command)
    context.system(*command)
    context.write('Projection Done')

    re = aims.Reader()
    ws = aims.Writer()

    texture_poles = re.read(self.pole.fullPath())
    tex_S16 = aims.TimeTexture_S16()
    tex_S16[0].assign(texture_poles[0])
    ws.write(tex_S16, self.pole.fullPath())
    context.system('AimsTextureDilation', '-i',self.white_mesh.fullPath(), '-t',self.pole.fullPath(),'-o',self.pole.fullPath(),'-s',self.dilation_1,'--connexity')#10
    context.system('AimsTextureErosion', '-i',self.white_mesh.fullPath(), '-t',self.pole.fullPath(),'-o',self.pole.fullPath(),'-s',self.erosion,'--connexity')#10
    context.system('AimsTextureDilation', '-i',self.white_mesh.fullPath(), '-t',self.pole.fullPath(),'-o',self.pole.fullPath(),'-s',self.dilation_1,'--connexity')#10
    context.write('Dilation Erosion Done')
    context.write('Topological correction...')
    mesh = re.read(self.white_mesh.fullPath())
    tex = re.read(self.pole.fullPath())
    cingular_tex_value = 1
    cingular_tex_clean, cing_tex_boundary = surfTls.textureTopologicalCorrection(mesh, tex[0].arraydata(), cingular_tex_value)
    tex_out = aims.TimeTexture_S16()
    tex_out[0].assign(cingular_tex_clean)
    ws.write(tex_out, self.pole.fullPath())
    context.write('... Done')