# -*- coding: utf-8 -*-
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

name = 'Registration Mesh'

userLevel = 0

signature = Signature(
      'mesh', ReadDiskItem( 'Hemisphere White Mesh', 'BrainVISA mesh formats' ),
      'fMRI_data', ReadDiskItem('4D Volume', 'BrainVISA volume formats'),
      'out_mesh', WriteDiskItem( 'Hemisphere White Mesh', 'BrainVISA mesh formats' )

)

def execution( self, context ):

      from soma import aims
      import numpy as np
      h = aims.Finder();
      meshpath =  self.mesh.fullPath()
      fctimage = self.fMRI_data.fullPath()
      
      # image fonctionnelle format NIFTI, dans son header 
      # y'a la matrice de passage "fct -> anat"
      # on prend de preference l'image moyenne
      
      fct_image = h.check(fctimage);
      assert(fct_image)
      fct_image = h.header();
      
      # image anatomique format NIFTI, on lit uniquement le header
      mesh_header = h.check(meshpath);
      assert(mesh_header)
      mesh_header = h.header();
      
      
      # matrice "mm anat ref. aims -> mm anat ref. spm"      
      trans_mesh = aims.Motion(mesh_header['transformations'][0]);
      
      # matrice "mm fct, ref. aims -> mm anat, ref. spm"
      trans_fct = aims.Motion(fct_image['transformations'][0]);
  
      # ouverture du maillage anatomique
      mesh = aims.read(meshpath);
      
      # inversion de la matrice de transformation "fct -> anat"
      # on obtient donc : "mm anat, ref. spm -> mm fct, ref. aims"
      itrans_fct = trans_fct.inverse()
      
      # taille des voxels fonctionnels, attention ici, il faut le convertir
      # en numpy array pour appliquer un diag plus bas
      fvox = np.array(fct_image['voxel_size']);
      
      # On va avoir besoin de diviser les coordonnees du maillage par la
      # taille des voxels un peu plus loin. On va creer un aims.motion.
      trans_fvox = aims.Motion( np.diag ( np.hstack( ( 1/fvox, [1] ) ) ).flatten() )
      
      # taille de l'image fonctionnelle
      fct_size = np.array( fct_image['volume_dimension'][0:3] ) 
      
      # pt_a_aims est en mm ref aims. faut le passer en ref spm
      aims.SurfaceManip.meshTransform( mesh, trans_mesh )
  
      # passage en mm fonctionnels :
      #aims.SurfaceManip.meshTransform( mesh, itrans_fct )
      
      # passage en voxels fonctionnels :
      aims.SurfaceManip.meshTransform( mesh, trans_fvox )
      
      aims.write ( mesh, self.out_mesh.fullPath() )
      context.write ( 'Finished' )