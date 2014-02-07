# -*- coding: utf-8 -*-
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
from freesurfer.brainvisaFreesurfer import launchFreesurferCommand
from glob import glob
from brainvisa.cortical_surface.surface_tools import texture_tools as textureTls
import numpy as np

name = 'Cingular Pole Texture From Freesurfer Label Texture'
userlevel = 2

# def validation():
#     testFreesurferCommand()

signature = Signature(
    'white_mesh', ReadDiskItem('White', 'Aims mesh formats', enableConversion=0),  
#    'side', Choice( ( 'left', 'lh' ), ( 'right', 'rh' ), None ),
    'database', ReadDiskItem( 'Directory', 'Directory' ),
    'subject', String(),
    'pole',WriteDiskItem( 'Hippocampus pole texture','Aims Texture formats' ),
    )

def initialization(self):
    def linkDB( self, dummy ):
        if self.white_mesh is not None:
            return self.white_mesh.get('_database')
    def linkSubject( self, dummy ):
        if self.white_mesh is not None:
             return self.white_mesh.get('subject')
#    def linkside( self, dummy ):
#        if self.white_mesh is not None:
#            return self.white_mesh.get( 'side' )
    self.linkParameters('pole', 'white_mesh')
    self.linkParameters('subject', 'white_mesh', linkSubject )
    self.linkParameters('database', 'white_mesh', linkDB )



def execution(self, context):
    if self.white_mesh.get('side')=='left':
        side = 'lh'
    elif self.white_mesh.get('side')=='right':
        side = 'rh'
    else:
        context.write('error in side')
        return
    
    re = aims.Reader()
    ws = aims.Writer()
    context.write('Reading textures and mesh')
    mesh = re.read(self.white_mesh.fullPath())

    context.write('Conversion of freesurfer labels to aims labels.')

    nbv = len(mesh.vertex())
    data = np.ones(nbv)

    fs_cortex_label = self.database.fullPath()+'/'+self.subject+'/label/'+side+'.cortex.label'
    f = open(fs_cortex_label,'r')
    lines = f.readlines()
    for a in lines[2:]:
        a_spl = a.split()
        data[int(a_spl[0])] = a_spl[4]
    f.close()

    context.write(data)
#    context.system('python', '-c', 'from freesurfer.freesurferTexture2Tex import freesurferTexture2TexBrainvisa as f; f(\"%s\", \"%s\", \"%s\");'%(fs_cortex_label, self.white_mesh.fullPath(), tmp_tex.fullPath()))


    context.write('Topological correction...')  
    cingular_tex_clean, cingular_tex_boundary = textureTls.textureTopologicalCorrection(mesh, data, 1)
    
    tex_out = aims.TimeTexture_S16()
    tex_out[0].assign(cingular_tex_clean)
    ws.write(tex_out, self.pole.fullPath())
    context.write('... Done')


