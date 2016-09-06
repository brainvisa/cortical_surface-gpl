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
from registration import getTransformationManager

name = 'Anatomist Show Pits As Spheres'
roles = ('viewer',)
userLevel = 0

def validation():
  anatomist.validation()

signature = Signature(
    'texture_pits', ReadDiskItem('pits texture', 'aims Texture formats'),
    'white_mesh',ReadDiskItem( 'Hemisphere White Mesh', 'aims mesh formats' ),
    'sphere_size', Float(),
)

def initialization( self ):
    self.linkParameters('white_mesh','texture_pits' )
    self.sphere_size = 1.0

def execution( self, context ):
    spheres_mesh_file = context.temporary(  'GIFTI file' )
    white_mesh = aims.read(self.white_mesh.fullPath())
    pits_texture = aims.read(self.texture_pits.fullPath())
    gen = aims.SurfaceGenerator()
    spheres_mesh = aims.AimsSurfaceTriangle()
    vert = np.array(white_mesh.vertex())  # vertex coordinates
    pits = np.where(np.array(pits_texture[0]))[0]
    for pit in pits:
        spheres_mesh += gen.sphere(vert[pit], self.sphere_size,10)
    #spheres_mesh.header()[ 'referentials' ] = white_mesh.header()[ 'referentials' ]
    transformManager = getTransformationManager()
    transformManager.copyReferential( self.white_mesh, spheres_mesh_file )
    aims.write(spheres_mesh, spheres_mesh_file.fullPath())
    a = anatomist.Anatomist()
    win = a.createWindow( 'Axial' )
    anamesh = a.loadObject( spheres_mesh_file.fullPath() )
    win.addObjects(anamesh )
    anamesh.setMaterial( a.Material(diffuse = [1.0, 0.0, 0.0, 1]) )
    anamesh2 = a.loadObject( self.white_mesh.fullPath() )
    win.addObjects(anamesh2 )
    #r = a.viewTextureOnMesh( self.white_mesh,
    #                                           self.texture_pits,
     #                                         a.getPalette('Talairach'),
      #                                        interpolation = 'rgb' )
    #context.write(r['window'])
    return [win, anamesh,anamesh2]
