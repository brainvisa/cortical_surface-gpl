#
#  Copyright (C) 2000-2001 INSERM
#  Copyright (C) 2000-2002 CEA
#  Copyright (C) 2000-2001 CNRS
#
#  This software and supporting documentation were developed by
#    INSERM U494
#    Hopital Pitie Salpetriere
#    91 boulevard de l'Hopital
#    75634 Paris cedex 13
#    France
#    --
#    CEA/DSV/SHFJ
#    4 place du General Leclerc
#    91401 Orsay cedex
#    France
#    --
#    CNRS UPR640-LENA
#    Hopital Pitie Salpetriere
#    47 boulevard de l'Hopital
#    75651 Paris cedex 13
#    France
#
#  $Id$
#

from neuroProcesses import *
import shfjGlobals
name = 'Send Texture To Atlas'
userLevel = 2


signature = Signature(
    'texture', ReadDiskItem('Texture', 'Texture' ),
    'mesh', ReadDiskItem('Mesh', 'MESH mesh' ),
    'mesh_x', ReadDiskItem( 'Coordinate Texture', 'Texture' ),
    'mesh_y', ReadDiskItem( 'Coordinate Texture', 'Texture' ),
    'atlasMesh', ReadDiskItem('Mesh', 'MESH mesh' ),
    'atlas_x', ReadDiskItem( 'Coordinate Texture', 'Texture' ),
    'atlas_y', ReadDiskItem( 'Coordinate Texture', 'Texture' ),
    'period_x', Float(),
    'textureOut', WriteDiskItem('Texture', 'Texture' )
)

def initialization( self ):
     self.period_x = 360.0

def execution ( self, context ):
    resample = [ 'AimsTextureToAtlas',
                 '-i', self.mesh.fullPath(),
                 '-t', self.texture.fullPath(),
                 '-o', self.textureOut.fullPath(),
                 '-a', self.atlasMesh.fullPath(),
                 '-ix', self.mesh_x.fullPath(),
                 '-iy', self.mesh_y.fullPath(),
                 '-ax', self.atlas_x.fullPath(),
                 '-ay', self.atlas_y.fullPath(),
                 '-px', self.period_x]
    apply( context.system, resample )

