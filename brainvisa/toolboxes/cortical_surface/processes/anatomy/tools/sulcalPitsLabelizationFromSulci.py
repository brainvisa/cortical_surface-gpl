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
from soma import aims
import numpy as np

name = 'Sulcal Pits Labelization From Sulci'

userLevel = 2


signature = Signature(
    'sulci_voronoi',ReadDiskItem( 'sulci voronoi texture', 'aims Texture formats'),
    'basins_texture',ReadDiskItem( 'basins texture',  'Aims texture formats' ),
    'pits_texture',ReadDiskItem( 'pits texture',  'Aims texture formats' ),
    'white_mesh',ReadDiskItem( 'Hemisphere White Mesh' , 'Aims mesh formats' ),
    'labeled_basins_texture',WriteDiskItem( 'labeled basins texture',  'Aims texture formats' ),
    'labeled_pits_texture',WriteDiskItem( 'labeled pits texture',  'Aims texture formats' ),


)

def initialization( self ):
    self.linkParameters( 'white_mesh', 'sulci_voronoi' )
    self.linkParameters( 'basins_texture', 'sulci_voronoi' )
    self.linkParameters( 'pits_texture', 'sulci_voronoi' )
    self.linkParameters( 'labeled_basins_texture', 'sulci_voronoi' )
    self.linkParameters( 'labeled_pits_texture', 'sulci_voronoi' )


def execution( self, context ):
    def valCount(vect,vals):
        count = np.zeros(vals.shape)
        for i,v in enumerate(vals):
            count[i] = np.sum(vect == v)
        return count

    mesh = aims.read(self.white_mesh.fullPath())
    voronoi_tex = aims.read( self.sulci_voronoi.fullPath() )
    voronoi = np.array( voronoi_tex[0] )
    basins_tex = aims.read( self.basins_texture.fullPath() )
    basins = np.array( basins_tex[0] )
    pits_tex = aims.read( self.pits_texture.fullPath() )


    basins_labeled = np.zeros(basins.shape)
    basins_labels = np.unique(basins)
    context.write(basins_labels)
    context.write('labelling the basins based on majority voting')
    for b_lab in basins_labels:
        b_inds = np.where(basins == b_lab)[0]
        b_sulci_lab = voronoi[b_inds]
        u_b_sulci_lab =np.unique(b_sulci_lab)
        context.write(u_b_sulci_lab)
        bins = valCount(b_sulci_lab,u_b_sulci_lab)
        #bins = np.bincount(b_sulci_lab)
        #context.write(bins)
        perc = 100.0 * bins / len(b_inds)
        context.write(perc)
        basins_labeled[b_inds] = u_b_sulci_lab[np.argmax(perc)]
        #for sulci_lab in u_b_sulci_lab:


    tex_out = aims.TimeTexture_S16()
    tex_out[0].assign(basins_labeled)
    aims.write(tex_out, self.labeled_basins_texture.fullPath())

    context.write('labelling the pits')
    pits = np.array( pits_tex[0] )
    pits_ind = np.where(pits)[0]
    pits[pits_ind] = np.uint16(basins_labeled[pits_ind])
    tex_out = aims.TimeTexture_S16()
    tex_out[0].assign(pits)
    aims.write(tex_out, self.labeled_pits_texture.fullPath())

    # context.write('Done')










