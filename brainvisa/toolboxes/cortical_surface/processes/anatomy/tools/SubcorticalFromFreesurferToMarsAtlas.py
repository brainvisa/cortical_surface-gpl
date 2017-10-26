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

import os, sys
from brainvisa.processes import *
#from freesurfer.brainvisaFreesurfer import launchFreesurferCommand
#from glob import glob
from brainvisa.cortical_surface.surface_tools import texture_tools as textureTls
import numpy as np
from freesurfer.brainvisaFreesurfer import launchFreesurferCommand
from freesurfer.brainvisaFreesurfer import testFreesurferCommand


name = 'Subcortical From Freesurfer To MarsAtlas Parcellation'
userlevel = 2

# def validation():
#     testFreesurferCommand()

signature = Signature(
    'left_Parcels_volume', ReadDiskItem( 'Left Gyri Volume', 'Aims writable volume formats' ), #parcellation volume', 'Aims writable volume formats' ),
    'right_Parcels_volume', ReadDiskItem( 'Right Gyri Volume', 'Aims writable volume formats' ), #parcellation volume', 'Aims writable volume formats' ),
    'freesurfer_database', ReadDiskItem( 'Directory', 'Directory' ),
    'subject', String(),
    'complete_Parcels_volume',WriteDiskItem( 'Parcellation volume',  'Aims writable volume formats'),
    )

def initialization(self):
    def linkDB( self, dummy ):
        if self.left_Parcels_volume is not None:
            return self.left_Parcels_volume.get('_database')
    def linkSubject( self, dummy ):
        if self.left_Parcels_volume is not None:
             return self.left_Parcels_volume.get('subject')
#    def linkside( self, dummy ):
#        if self.white_mesh is not None:
#            return self.white_mesh.get( 'side' )
    self.linkParameters('complete_Parcels_volume', 'left_Parcels_volume')
    self.linkParameters('right_Parcels_volume', 'left_Parcels_volume')
    self.linkParameters('subject', 'left_Parcels_volume', linkSubject )
    self.linkParameters('freesurfer_database', 'left_Parcels_volume', linkDB )

def FS_convert(vol_in, vol_out, context):
    try:
        launchFreesurferCommand( context, '',
                           'mri_convert',
                           vol_in,
                           vol_out )
        #context.system('mri_convert',
        #            '-i', vol_in,
        #            '-o', vol_out )
        #os.system('mri_convert -i '+vol_in+' -o '+vol_out)
    except IOError:
        context.write('problem with mri_convert, file format conversion failed')

def copy_ROI(vol_hiphop_L, vol_hiphop_R, vol_aseg, labels, context):
    # fusion L+R
    vol_hiphop_R = vol_hiphop_R + vol_hiphop_L

    asegData = vol_aseg.arraydata()
    hiphopData = np.array(vol_hiphop_R , copy=False )#vol_hiphop_R.arraydata()

    for lab in labels:
        context.write( 'freesurfer aseg label '+str(lab) )
        inds = np.where(asegData==lab)
        #inds = asegData==lab
        #print inds
        #print hiphopData.shape
        hiphopData[inds[3],inds[2],inds[1],0] = 200+lab
        #hiphopData[0,inds[0],inds[1],inds[2]] = lab
    return hiphopData



def execution(self, context):


    context.write( 'run freesurfer_setup before runing this code' )
    labels = [49,50,51,52,53,54,58,10,11,12,13,17,18,26]
    vol_aseg_mgz = os.path.join(self.freesurfer_database.fullPath(),self.subject,'mri','aseg.mgz')
    vol_aseg_nii = os.path.join(self.freesurfer_database.fullPath(),self.subject,'mri','aseg.nii')
    FS_convert(vol_aseg_mgz, vol_aseg_nii, context)
    context.write( 'file format conversion done, doing the fusion')
    vol_aseg = aims.read(vol_aseg_nii)
    vol_parcels_L = aims.read(self.left_Parcels_volume.fullPath())
    vol_parcels_R = aims.read(self.right_Parcels_volume.fullPath())
    complete_vol = copy_ROI(vol_parcels_L, vol_parcels_R, vol_aseg, labels, context)

    outVol = aims.Volume_S16(complete_vol)

    aims.write( outVol, self.complete_Parcels_volume.fullPath() )
    context.write('... Done')


