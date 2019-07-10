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
from brainvisa import anatomist
from soma import aims
import numpy as np
from brainvisa.cortical_surface.surface_tools import texture_tools as texTls

name = 'Anatomist Show Texture Extrema As Spheres'
roles = ('viewer',)
userLevel = 0
allowed_processes = ()

def validation():
  anatomist.validation()

signature = Signature(
    'texture', ReadDiskItem('texture', 'aims Texture formats'),
    'white_mesh',ReadDiskItem( 'Hemisphere White Mesh', 'aims mesh formats' ),
    'sphere_size', Float(),
)

def initialization( self ):
    self.linkParameters('white_mesh','texture' )
    self.sphere_size = 1.0

def execution( self, context ):
    context.write('use temporary files...')
    spheres_mesh_file = context.temporary(  'GIFTI file' )
    spheres_texture_file = context.temporary(  'GIFTI file' )
    white_mesh = aims.read(self.white_mesh.fullPath())
    texture = aims.read(self.texture.fullPath())
    atex = np.array(texture[0])
    gen = aims.SurfaceGenerator()
    spheres_mesh = aims.AimsSurfaceTriangle()
    vert = np.array(white_mesh.vertex())  # vertex coordinates
    extrema_tex = texTls.TextureExtrema(white_mesh, atex)
    extrema = np.where(extrema_tex)[0]
    spheres_texture = list()
    nb_vert_s = 64
    for ex in extrema:
        spheres_mesh += gen.sphere(vert[ex], self.sphere_size,10)
        spheres_texture.extend(extrema_tex[ex]*np.ones((nb_vert_s,1),np.int16))

    a_tex_out = np.array(spheres_texture)
    tex_out = aims.TimeTexture_S16()
    tex_out[0].assign(a_tex_out)
    a = anatomist.Anatomist()
    win = a.createWindow( 'Axial' )
    objects = [win]
    if len(spheres_mesh.vertex()) != 0:
        anamesh = a.toAObject(spheres_mesh)
        anamesh.releaseAppRef()
        anatex = a.toAObject(tex_out)
        anatex.releaseAppRef()
        anatex.setPalette('GREEN-RED-ufusion')
        fusionTexSurf = a.fusionObjects([anamesh, anatex],
                                        method='FusionTexSurfMethod')
        fusionTexSurf.releaseAppRef()
        win.addObjects( fusionTexSurf )
        objects += [fusionTexSurf]
    anamesh2 = a.toAObject(white_mesh)
    anamesh2.releaseAppRef()
    anatex2 = a.toAObject(texture)
    anatex2.releaseAppRef()
    anatex2.setPalette('Purple-Red + Stripes')
    fusionTexSurf2 = a.fusionObjects([anamesh2, anatex2],
                                     method='FusionTexSurfMethod')
    fusionTexSurf2.releaseAppRef()
    objects.append(fusionTexSurf2)
    win.addObjects( fusionTexSurf2 )
    return objects
    #return [win, anamesh,anatex,fusionTexSurf,anamesh2,fusionTexSurf2]

