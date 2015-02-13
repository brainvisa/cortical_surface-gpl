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


def validation():
    try:
        from brainvisa.cortical_surface.surface_tools import PDE_tools as pdeTls
    except:
        raise ValidationError( 'brainvisa.cortical_surface.parameterization.surface_tools module can not be imported.' )
  

from brainvisa.processes import *
try:
    from brainvisa.cortical_surface.surface_tools import PDE_tools as pdeTls
except:
    pass

name = 'Depth Potential Function'
userLevel = 0

# Argument declaration
signature = Signature(
    'input_mesh',ReadDiskItem( 'Hemisphere White Mesh' , 'Aims mesh formats' ),
    'DPF_texture',WriteDiskItem( 'DPF texture',
    'Aims texture formats' ),
    'alphas', ListOf( Float() ),
)


# Default values
def initialization( self ):
    self.linkParameters( 'DPF_texture', 'input_mesh' )
    self.alphas = [0.03]


def execution( self, context ):
    re = aims.Reader()
    ws = aims.Writer()
    
    mesh = re.read(self.input_mesh.fullPath())
    tmp_curv_tex = context.temporary(  'Texture' )
    context.write('computing mean curvature of the mesh')
    context.system( 'AimsMeshCurvature', '-i', self.input_mesh, '-o', tmp_curv_tex.fullPath() , '-m', 'fem' )
    curv = re.read(tmp_curv_tex.fullPath())
    k = curv[0].arraydata()

    #alphas=[0.03] #[0.00,0.025,0.035,0.07,0.09]
    tex_dpf = aims.TimeTexture_FLOAT()
    context.write('computing DPF')
    context.write(self.alphas)
    dpf = pdeTls.depthPotentialFunction(mesh, k, self.alphas)
    for ind,alpha in enumerate(self.alphas):
        print 'alpha :',alpha
        tex_dpf[ind].assign(dpf[ind])

    ws.write(tex_dpf, self.DPF_texture.fullPath())
    context.write('... Done')
