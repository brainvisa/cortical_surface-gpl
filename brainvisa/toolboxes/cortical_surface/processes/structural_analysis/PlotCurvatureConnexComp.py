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

name = 'Plot Curvature Connex Components'

userLevel = 3

signature = Signature(
    'white', ReadDiskItem( 'Hemisphere White Mesh', 'MESH mesh'),
    'scale_space', ReadDiskItem( 'Scale Space White Curvature Texture', 'Texture')
    )

def initialization( self ):

def execution( self, context ):
     sujet = self.white.get( 'subject')
     context.write(sujet)

    # on lit la texture d'espace échelle et on crée un fichier temporaire
    scale_space = aims.read(str(self.scale_space))
    curvature_texfile = context.temporary("Texture")

    nb_scales = len(scale_space)
    print str(nb_scales) + " scales"

    # pour chaque niveau d'échelle
    for _ in xrange(nb_scales) :
      # on recupere la texture
      
      # on calcule une carte seuillee


    
     
     #call_list = [ 'surfLabelsTex2Graph',
                #'-m', self.white,
                #'-t', self.texture,                   
                #'-g', self.graph,
                #'--lat', self.latitude,
                #'--lon', self.longitude,
                #'-s', s,
                #'--repM', self.repMesh]
     context.write ( "Executing process.." )
     apply(context.system, call_list)
     context.write('Finished')

