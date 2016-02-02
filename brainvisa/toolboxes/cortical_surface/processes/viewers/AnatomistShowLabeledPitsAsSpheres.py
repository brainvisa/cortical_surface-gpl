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

from neuroProcesses import *
import shfjGlobals
from brainvisa import anatomist
from soma import aims
import numpy as np

name = 'Anatomist Show Labeled Pits As Spheres'
roles = ('viewer',)
userLevel = 0

def validation():
  anatomist.validation()

signature = Signature(
    'texture_labeled_pits', ReadDiskItem('labeled pits texture', 'aims Texture formats'),
    'white_mesh',ReadDiskItem( 'Hemisphere White Mesh', 'aims mesh formats' ),
    'sphere_size', Float(),
)

def initialization( self ):
    self.linkParameters('white_mesh','texture_labeled_pits' )
    self.sphere_size = 1.0

def execution( self, context ):
    context.write('use temporary files...')
    spheres_mesh_file = context.temporary(  'GIFTI file' )
    spheres_texture_file = context.temporary(  'GIFTI file' )
    white_mesh = aims.read(self.white_mesh.fullPath())
    pits_texture = aims.read(self.texture_labeled_pits.fullPath())
    apits_texture = np.array(pits_texture[0])
    gen = aims.SurfaceGenerator()
    spheres_mesh = aims.AimsSurfaceTriangle()
    vert = np.array(white_mesh.vertex())  # vertex coordinates
    pits = np.where(apits_texture)[0]
    spheres_texture = list()
    nb_vert_s = 64
    for pit in pits:
        spheres_mesh += gen.sphere(vert[pit], self.sphere_size,10)
        spheres_texture.extend(apits_texture[pit]*np.ones((nb_vert_s,1),np.int16))
    aims.write(spheres_mesh, spheres_mesh_file.fullPath())

    a_tex_out = np.array(spheres_texture)
    tex_out = aims.TimeTexture_S16()
    tex_out[0].assign(a_tex_out)
    aims.write(tex_out, spheres_texture_file.fullPath())
    a = anatomist.Anatomist()
    win = a.createWindow( 'Sagittal' )
    anamesh = a.loadObject( spheres_mesh_file.fullPath() )
    anatex = a.loadObject( spheres_texture_file.fullPath() )
    anapalette = a.getPalette('Graph-Label')
    anatex.setPalette( anapalette, minVal = 0, maxVal= 2.04)
    fusionTexSurf = a.fusionObjects( [anamesh, anatex], method='FusionTexSurfMethod' )
    win.addObjects(fusionTexSurf )
    anamesh2 = a.loadObject( self.white_mesh.fullPath() )
    win.addObjects(anamesh2 )
    return [win, anamesh,anamesh2,anatex,fusionTexSurf]

