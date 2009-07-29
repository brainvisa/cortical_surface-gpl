# Copyright CEA and IFR 49 (2000-2005)
#
#  This software and supporting documentation were developed by
#      CEA/DSV/SHFJ and IFR 49
#      4 place du General Leclerc
#      91401 Orsay cedex
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

name = 'Create Cortical Thickness Texture'

userLevel = 2

signature = Signature(
      'white_mesh', ReadDiskItem( 'Hemisphere White Mesh', 'MESH mesh' ),
      'hemi_mesh', ReadDiskItem( 'Hemisphere Mesh', 'MESH mesh'),
      'neighbours_white', ReadDiskItem( 'Any Type', 'Text File'),
      'neighbours_hemi', ReadDiskItem( 'Any Type', 'Text File'),
      'outtex_white', WriteDiskItem( 'Cortical thickness', 'Texture', requiredAttributes={ 'mesh': 'white' }),
      'outtex_hemi', WriteDiskItem( 'Cortical thickness', 'Texture', requiredAttributes={ 'mesh': 'hemi' }),
      'write_median_mesh', Boolean(),
      'median_mesh', WriteDiskItem( 'Median Mesh', 'MESH mesh')
)

def initialization( self ):
      self.setOptional('median_mesh','neighbours_white', 'neighbours_hemi')
      self.linkParameters(  'hemi_mesh', 'white_mesh' )
      self.linkParameters('outtex_white', 'white_mesh' )
      self.linkParameters('median_mesh', 'white_mesh' )
      self.linkParameters( 'outtex_hemi', 'white_mesh' )
      
      self.write_median_mesh = True
     

def execution( self, context ):
      context.write('Creating a cortical thickness texture...')

      voisins_int_path = ""
      voisins_ext_path = ""
      varintmesh = self.white_mesh
      varextmesh = self.hemi_mesh
      medianpath = ""
      if self.write_median_mesh is True :
         medianpath = self.median_mesh
      
      if self.neighbours_white is not None :
         voisins_int_path = self.neighbours_white
         context.write("Using provided neighbours table (for G/W mesh)...")
         
      if self.neighbours_hemi is not None :
         voisins_ext_path = self.neighbours_hemi
         context.write("Using provided neighbours table (for G/CSF mesh)...")
         
      context.write("Process i/e :" , varintmesh , ";" , varextmesh , ";" , self.neighbours_hemi , ";" , self.outtex_white )
      createInt2Extmesh0 = [ 'AimsCorticalThickness',
         '-i', varintmesh,
         '-e', varextmesh,
         '-v', voisins_ext_path,
         '-d', 0,
         '-o', self.outtex_white,
         '-m', medianpath
      ]
      apply( context.system, createInt2Extmesh0 )

      context.write("Process e/i :" , varextmesh , ";" , varintmesh , ";" , self.neighbours_white, ";" , self.outtex_hemi )

      createExt2Intmesh0 = [ 'AimsCorticalThickness',
         '-i', varextmesh,
         '-e', varintmesh,
         '-v', voisins_int_path,
         '-d', 1,
         '-o', self.outtex_hemi
      ]
      apply( context.system, createExt2Intmesh0 )

      
      

#       i=1
#       while (i<self.iterations):
#          j=i+1
#          varintmesh = tmpInt2Ext
#          varextmesh = tmpExt2Int
#          tmpInt2Ext = context.temporary( 'MESH mesh')
#          tmpExt2Int = context.temporary( 'MESH mesh')
#          if i == self.iterations - 1 :
#             tmpInt2Ext=self.outmesh
#          context.write(j , "eme process i/e :" , varintmesh , ";" , varextmesh , ";" , tmpVoisinsExt , ";" , tmpInt2Ext )
#          createInt2Extmesh = [ 'AimsMeshMedianSurface',
#             '-i', varintmesh,
#             '-e', varextmesh,
#             '-v', tmpVoisinsExt,
#             '-d', 0,
#             '-o', tmpInt2Ext
#          ]
#          apply( context.system, createInt2Extmesh )
#          if i < self.iterations - 1 :
#             context.write(j , "eme process e/i :" , varextmesh , ";" , varintmesh , ";" , tmpVoisinsInt , ";" , tmpExt2Int )
#             createExt2Intmesh = [ 'AimsMeshMedianSurface',
#                '-i', varextmesh,
#                '-e', varintmesh,
#                '-v', tmpVoisinsInt,
#                '-d', 1,
#                '-o', tmpExt2Int
#             ]
#             apply( context.system, createExt2Intmesh )
#          i=i+1
      
      context.write('Finished')

