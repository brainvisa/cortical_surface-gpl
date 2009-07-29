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

name = 'Create Mesh Median Surface'

userLevel = 2

signature = Signature(
      'white_mesh', ReadDiskItem( 'Hemisphere White Mesh', 'MESH mesh' ),
      'hemi_mesh', ReadDiskItem( 'Hemisphere Mesh', 'MESH mesh'),
      'neighbours_white', ReadDiskItem( 'Any Type', 'Text File'),
      'neighbours_hemi', ReadDiskItem( 'Any Type', 'Text File'),
      'median_mesh', WriteDiskItem( 'Median Mesh', 'MESH mesh'),
      'iterations', Integer()
)

def initialization( self ):
      self.setOptional('neighbours_white', 'neighbours_hemi', 'iterations')
      self.linkParameters( 'hemi_mesh', 'white_mesh' )
      self.linkParameters(  'median_mesh', 'white_mesh')
      self.iterations = 1
     

def execution( self, context ):
      context.write('Creating a median surface...')
      
      tmpVoisinsInt = ""
      tmpVoisinsExt = ""
      varintmesh = self.white_mesh
      varextmesh = self.hemi_mesh
      voisins_int_path = ""
      voisins_ext_path = ""
      if self.neighbours_white is not None :
         voisins_int_path = self.neighbours_white
         context.write("Using provided neighbours table (for G/W mesh)...")
      else :
         self.neighbours_white = ""
         if self.iterations > 1 :
            tmpVoisinsInt = context.temporary( 'Text File' )
         
      if self.neighbours_hemi is not None :
         voisins_ext_path = self.neighbours_hemi
         context.write("Using provided neighbours table (for G/CSF mesh)...")
      else :
         self.neighbours_hemi = ""
         if self.iterations > 1 :
            tmpVoisinsExt = context.temporary( 'Text File' )

      tmpExt2Int=context.temporary( 'MESH mesh')
      if self.iterations > 1 : 
         tmpInt2Ext=context.temporary( 'MESH mesh')         
      else :
         tmpInt2Ext=self.median_mesh
         
      context.write("Premier process i/e :" , varintmesh , ";" , varextmesh , ";" , self.neighbours_hemi , ";", tmpVoisinsExt , ";" , tmpInt2Ext )
      createInt2Extmesh0 = [ 'AimsMeshMedianSurface',
         '-i', varintmesh,
         '-e', varextmesh,
         '-v', voisins_ext_path,
         '-d', 0,
         '-vout', tmpVoisinsExt,
         '-o', tmpInt2Ext
      ]
      apply( context.system, createInt2Extmesh0 )

      context.write("Premier process e/i :" , varextmesh , ";" , varintmesh , ";" , self.neighbours_white , ";", tmpVoisinsInt , ";" , tmpExt2Int )
      if self.iterations > 1 :
         createExt2Intmesh0 = [ 'AimsMeshMedianSurface',
            '-i', varextmesh,
            '-e', varintmesh,
            '-v', voisins_int_path,
            '-d', 1,
            '-vout', tmpVoisinsInt,
            '-o', tmpExt2Int
         ]
         apply( context.system, createExt2Intmesh0 )
      if tmpVoisinsInt == "" :
         tmpVoisinsInt = self.neighbours_white
      if tmpVoisinsExt == "" :
         tmpVoisinsExt = self.neighbours_hemi
            
      i=1
      while (i<self.iterations):
         j=i+1
         varintmesh = tmpInt2Ext
         varextmesh = tmpExt2Int
         tmpInt2Ext = context.temporary( 'MESH mesh')
         tmpExt2Int = context.temporary( 'MESH mesh')
         if i == self.iterations - 1 :
            tmpInt2Ext=self.median_mesh
         context.write(j , "eme process i/e :" , varintmesh , ";" , varextmesh , ";" , tmpVoisinsExt , ";" , tmpInt2Ext )
         createInt2Extmesh = [ 'AimsMeshMedianSurface',
            '-i', varintmesh,
            '-e', varextmesh,
            '-v', tmpVoisinsExt,
            '-d', 0,
            '-o', tmpInt2Ext
         ]
         apply( context.system, createInt2Extmesh )
         if i < self.iterations - 1 :
            context.write(j , "eme process e/i :" , varextmesh , ";" , varintmesh , ";" , tmpVoisinsInt , ";" , tmpExt2Int )
            createExt2Intmesh = [ 'AimsMeshMedianSurface',
               '-i', varextmesh,
               '-e', varintmesh,
               '-v', tmpVoisinsInt,
               '-d', 1,
               '-o', tmpExt2Int
            ]
            apply( context.system, createExt2Intmesh )
         i=i+1
      
      context.write('Finished')

