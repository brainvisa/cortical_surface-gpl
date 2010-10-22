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

from neuroProcesses import *
import shfjGlobals

name = 'Projection using Convolution Kernels'

userLevel = 0

signature = Signature(
      'white_mesh', ReadDiskItem( 'Hemisphere White Mesh', 'BrainVISA mesh formats' ),
      'kernels', ReadDiskItem('Projection convolution kernels', 'BrainVISA volume formats'),
      'fMRI_4D_data', ReadDiskItem('4D Volume', 'BrainVISA volume formats'),
      'fMRI_surface_data', WriteDiskItem( 'Functional Time Texture', 'Texture')
)


def initialization ( self ):
      #self.linkParameters('white_mesh', 'anatomy')

      self.linkParameters('kernels', 'white_mesh')
      

def execution( self, context ):

      context.write ( volume_path )

      projection = [ 
      'AimsFunctionProjection', 
      '-op', '1',
      '-d', self.kernels.fullPath(),
      '-d1', volume_path,
      '-m', self.white_mesh.fullPath(),
      '-o', self.fMRI_surface_data.fullPath()
      ]
      apply(context.system, projection)
  
     
      #[ mesh, kernels, fctimage, outputtexture] =  [self.white_mesh.fullPath(), self.kernels.fullPath(), self.fMRI_4D_data.fullPath(), self.fMRI_surface_data.fullPath()]
      
      
      
      
      #r = aims.Reader();
      #w = aims.Writer();
      #h = aims.Finder();
      ## image fonctionnelle format NIFTI, dans son header 
      ## y'a la matrice de passage "fct -> anat"
      ## on prend de preference l'image moyenne
      
      #fct_image = h.check(fctimage);
      #assert(fct_image)
      #fct_image = h.header();
      
      ## image anatomique format NIFTI, on lit uniquement le header
      #mesh_header = h.check(mesh);
      #assert(mesh_header)
      #mesh_header = h.header();
      
      
      ## matrice "mm anat ref. aims -> mm anat ref. spm"      
      ##trans_mesh = aims.Motion(mesh_header['transformations'][0]);
      
      ## matrice "mm fct, ref. aims -> mm anat, ref. spm"
      ##trans_fct = aims.Motion(fct_image['transformations'][0]);
  
      ## ouverture du maillage anatomique
      ##mesh = r.read(mesh);
      
      ## inversion de la matrice de transformation "fct -> anat"
      ## on obtient donc : "mm anat, ref. spm -> mm fct, ref. aims"
      ##itrans_fct = trans_fct.inverse()
      
      ## taille des voxels fonctionnels, attention ici, il faut le convertir
      ## en numpy array pour appliquer un diag plus bas
      ##fvox = np.array(fct_image['voxel_size']);
      
      ## On va avoir besoin de diviser les coordonnees du maillage par la
      ## taille des voxels un peu plus loin. On va creer un aims.motion.
      #trans_fvox = aims.Motion( np.diag ( np.hstack( ( 1/fvox, [1] ) ) ).flatten() )
      
      ## taille de l'image fonctionnelle
      #fct_size = np.array( fct_image['volume_dimension'][0:3] ) 
      
      ## pt_a_aims est en mm ref aims. faut le passer en ref spm
      #aims.SurfaceManip.meshTransform( mesh, trans_mesh )
  
      ## passage en mm fonctionnels :
      ##aims.SurfaceManip.meshTransform( mesh, itrans_fct )
      
      ## passage en voxels fonctionnels :
      #aims.SurfaceManip.meshTransform( mesh, trans_fvox )
      
      ## chargement des KERNELS
      #kernels = r.read(kernels)
      #k = kernels.arraydata();
      #functional_mesh = np.array(mesh.vertex())
      
      ## parfois l'image fonctionnelle est coupee en haut, le maillage
      ## risque donc de se retrouver en dehors de l'image. on va calculer combien
      ## lignes/collones il manque, et on va les rajouter 
      ## on recupere les valeur max des points pour voir ceux qui sont en dehors
      ## de l'image fonctionnelle
      #maxpoints = [functional_mesh[:,x].max() for x in range(3)]
      #minpoints = [functional_mesh[:,x].min() for x in range(3)]
      #addinf = np.array([abs(min(int(minpoints[x]-4),0)) for x in range(3)])
      #addsup = np.array([int(max(0,maxpoints[x]+3-(fct_size[x]-1))) for x in range(3)])
      
      #nf = []
      #fimg = aims.read(fctimage).arraydata()
      #print "--", fctimage
      #for j, img in enumerate(fimg):
          #print j,"/",len(fimg)
          ##print img.shape
          ## lecture de l'image fonctionnelle a projeter :
          #d = img.transpose()
          ##print d.shape
          ## (le transpose sert a passer de ZYX a XYZ)
          
          ## on encapsule les donnees de l'image fonctionnelle dans un nouveau
          ## tableau plus grand
          ## TODO : bon, pas super optimal ca ... penser a faire mieux.
          #A = (addsup+addinf+fct_size)[0]
          #B = (addsup+addinf+fct_size)[1]
          #C = (addsup+addinf+fct_size)[2]
          #n = np.zeros((addsup+addinf+fct_size),float)
          #n[addinf[0]:A-addsup[0],addinf[1]:B-addsup[1],addinf[2]:C-addsup[2]]=d
          
          ## la projection
          #f = np.zeros(len(mesh.vertex()))
          #for i in range(len(mesh.vertex())):
              ##if (i%2500)==0:
              ##    print i;
              #[x,y,z] = np.array(mesh.vertex()[i]).astype(int)+addinf;
              ## une multiplication de matrices toute bete
              #try:
                  #f[i]=np.sum(np.sum(np.sum(n[x-3:x+4,y-3:y+4,z-3:z+4]*k[i])))
              #except IndexError:
                  #print x, y, z
                  #raise "toto"
            ## TODO: utiliser la taille des matrices des kernels ici !!!
          #nf.append(f)
	  #if j > 5:
	    #break

      
      #nf = np.array(nf)
      #print np.shape(nf)
      #print len(nf)
      #T = aims.TimeTexture_FLOAT(len(nf), len(nf[0]) )
      #for i, n in enumerate(nf):
          #t = aims.Texture_FLOAT(np.array(n).astype(np.float32))
          #for j in xrange(len(t)):
              #T[i][j] = t[j]
      
      #context.write('Finished')
      #aims.write( T, str(self.fMRI_surface_data.fullPath()) )

