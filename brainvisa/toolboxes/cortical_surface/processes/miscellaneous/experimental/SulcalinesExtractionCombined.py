# -*- coding: utf-8 -*-
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
from __future__ import absolute_import
from brainvisa.processes import *
import numpy as np
from soma import aims

try:
    from brainvisa.cortical_surface.surface_tools import readSulcusLabelTranslationFile as rSLT
except:
  pass


name = 'HESCHL Sulcal Lines Extraction'

userLevel = 2

#def validation():
#    anatomist.validation()
    
signature = Signature(

    'white_sulcalines',ReadDiskItem( 'hemisphere Sulcal Lines texture', 'aims Texture formats'),
    'side', Choice('left', 'right'),
    'sulcalines_correspondance',ReadDiskItem( 'Graph Label Translation', 'Text File'),
    'white_sulci',ReadDiskItem( 'Sulci White Texture' ,'aims texture formats'),
    'sulci_correspondance',ReadDiskItem( 'Sulci To White Texture Translation', 'Text File'),
    'sulci_to_copy',ListOf( String() ),
    'combined_white_sulcalines',WriteDiskItem( 'hemisphere Sulcal Lines texture', 'aims Texture formats'),
#    'combined_sulcalines_correspondance',WriteDiskItem( 'Constraint coordinates values', 'Text File'),
    #'constraint_value', Choice('Basins Label','Lat/Lon'),
)

def initialization( self ):
    def linkSide( proc, dummy ):
        if proc.white_sulcalines is not None:
            side = proc.white_sulcalines.get( 'side' )
            return side
    self.linkParameters( 'sulcalines_correspondance','white_sulcalines' )
    self.linkParameters( 'side', 'white_sulcalines', linkSide )
    self.linkParameters( 'white_sulci','white_sulcalines' )
    self.linkParameters( 'sulci_correspondance','white_sulcalines' )
    self.linkParameters( 'combined_white_sulcalines','white_sulcalines' )
    self.sulci_to_copy = ['S.Heschl.','F.C.L.r.asc.','F.Cal.ant.-Sc.Cal.']
#    self.linkParameters( 'combined_sulcalines_correspondance','white_sulcalines' )

def execution( self, context ):

    tex_white_sulcalines = aims.read(self.white_sulcalines.fullPath())
    sulcalines_labels_dict = rSLT.readSulcusLabelTranslationFile(self.sulcalines_correspondance.fullPath(), invert=True)
    tex_white_sulci = aims.read(self.white_sulci.fullPath())
    sulci_labels_dict = rSLT.readSulcusLabelTranslationFile(self.sulci_correspondance.fullPath(), invert=True)
    sulci_to_copy = [sulc+'_'+self.side for sulc in self.sulci_to_copy]
    labels_to_copy = [sulcalines_labels_dict[sulc] for sulc in sulci_to_copy]
    labels_to_copy_sulci = [sulci_labels_dict[sulc] for sulc in sulci_to_copy]
    context.write('the labels of the sulci to copy is :')
    context.write(str(labels_to_copy)+' in the sulcal lines texture')
    context.write(str(labels_to_copy_sulci)+' in the sulci texture')
    if len(labels_to_copy) < 1 or len(labels_to_copy) != len(labels_to_copy_sulci):
        context.write('the sulcus cannot be found in both textures, aborting!')
    else:
        context.write('the sulcus was found in both textures, good!')
        atex_white_sulcalines = tex_white_sulcalines[0].arraydata()
        atex_white_sulci = tex_white_sulci[0].arraydata()
        for ind,lab in enumerate(labels_to_copy):
            inds_sulcalines = np.where(atex_white_sulcalines == lab)
            inds_sulci = np.where(atex_white_sulci == labels_to_copy_sulci[ind])
            "deleting the projection from sulcalines"
            atex_white_sulcalines[inds_sulcalines] = 0
            "copying the projection from sulci"
            atex_white_sulcalines[inds_sulci] = lab
        tex_out = aims.TimeTexture_S16()
        tex_out[0].assign(atex_white_sulcalines)
        aims.write(tex_out, self.combined_white_sulcalines.fullPath())
        context.write('Done')
