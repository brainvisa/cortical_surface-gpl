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
include( 'base' )

insert( '{protocol}/{subject}',
  'surface', SetContent(
    '<subject>_Lhippo', SetType( 'Left hippocampus pole texture' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rhippo', SetType( 'Right hippocampus pole texture' ), SetWeakAttr( 'side', 'right' ),
    '<subject>_L_poles', SetType( 'Left poles texture' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_R_poles', SetType( 'Right poles texture' ), SetWeakAttr( 'side', 'right' ),
    '<subject>_L_lat_cleaned', SetType( 'Left hemisphere latitude cleaned constraints texture'), SetWeakAttr( 'side', 'left' ),
    '<subject>_R_lat_cleaned', SetType( 'Right hemisphere latitude cleaned constraints texture'), SetWeakAttr( 'side', 'right' ),
    '<subject>_L_lon_cleaned', SetType( 'Left hemisphere longitude cleaned constraints texture'), SetWeakAttr( 'side', 'left' ),
    '<subject>_R_lon_cleaned', SetType( 'Right hemisphere longitude cleaned constraints texture'), SetWeakAttr( 'side', 'right' ),
    '<subject>_L_lat_roots', SetType( 'Left hemisphere latitude constraints texture'), SetWeakAttr( 'side', 'left' ),
    '<subject>_R_lat_roots', SetType( 'Right hemisphere latitude constraints texture'), SetWeakAttr( 'side', 'right' ),
    '<subject>_L_lon_roots', SetType( 'Left hemisphere longitude constraints texture'), SetWeakAttr( 'side', 'left' ),
    '<subject>_R_lon_roots', SetType( 'Right hemisphere longitude constraints texture'), SetWeakAttr( 'side', 'right' ),
    '<subject>_L_lat', SetType( 'Left hemisphere latitude texture'), SetWeakAttr( 'side', 'left' ),
    '<subject>_R_lat', SetType( 'Right hemisphere latitude texture'), SetWeakAttr( 'side', 'right' ),
    '<subject>_L_lon', SetType( 'Left hemisphere longitude texture'), SetWeakAttr( 'side', 'left' ),
    '<subject>_R_lon', SetType( 'Right hemisphere longitude texture'), SetWeakAttr( 'side', 'right' ),
    '<subject>_rootsValues', SetType( 'Constraint coordinates values'),
    '<subject>_L_grid', SetType( 'Left hemisphere coordinate grid'), SetWeakAttr( 'side', 'left' ),
    '<subject>_R_grid', SetType( 'Right hemisphere coordinate grid'), SetWeakAttr( 'side', 'right' ),
    '<subject>_L_gyri', SetType( 'Left hemisphere gyri parcellation texture'), SetWeakAttr( 'side', 'left' ),
    '<subject>_R_gyri', SetType( 'Right hemisphere gyri parcellation texture'), SetWeakAttr( 'side', 'right' ),
    '<subject>_L_gyri_regul', SetType( 'Left hemisphere regularized parcellation texture'), SetWeakAttr( 'side', 'left' ),
    '<subject>_R_gyri_regul', SetType( 'Right hemisphere regularized parcellation texture'), SetWeakAttr( 'side', 'right' ),
    '<subject>_Lhippo_Volume', SetType( 'Left Cingular Pole Template Subject' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rhippo_Volume', SetType( 'Right Cingular Pole Template Subject' ), SetWeakAttr( 'side', 'right' ),
    '<subject>_Talairach_To_Subject_Transformation', SetType( 'Talairach To Subject Transformation'),
    '<subject>_Subject_To_Template_Transformation', SetType( 'Subject To Template Transformation' ),
    '<subject>_L_gyriGraph', SetType( 'Left Gyri Graph' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_R_gyriGraph', SetType( 'Right Gyri Graph' ), SetWeakAttr( 'side', 'right' ),
    '<subject>_L_gyriVolume', SetType( 'Left Gyri Volume' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_R_gyriVolume', SetType( 'Right Gyri Volume' ), SetWeakAttr( 'side', 'right' ),
    '<subject>_Lwhite_KERNEL', SetType( 'Projection convolution kernels'), SetWeakAttr( 'side', 'left' ) ,
    '<subject>_Rwhite_KERNEL', SetType( 'Projection convolution kernels'), SetWeakAttr( 'side', 'right' ) ,
    '{volume}_<subject>_Lwhite_projection', SetType( 'Functional texture'), SetWeakAttr( 'side', 'left' ) ,
    '{volume}_<subject>_Rwhite_projection', SetType( 'Functional texture'), SetWeakAttr( 'side', 'right' ),
    '<subject>_Lwhite_thickness', SetType('Cortical thickness'), SetWeakAttr('side', 'left', 'mesh', 'white'),
    '<subject>_Rwhite_thickness', SetType('Cortical thickness'), SetWeakAttr('side', 'right', 'mesh', 'white'),
    '<subject>_Lhemi_thickness', SetType('Cortical thickness'), SetWeakAttr('side', 'left', 'mesh', 'hemi'),
    '<subject>_Rhemi_thickness', SetType('Cortical thickness'), SetWeakAttr('side', 'right', 'mesh', 'hemi'),
    "<subject>_Lwhite_curv_blobs", SetType( 'Blob White Curvature Texture' ), SetWeakAttr( 'side', 'left' ),
    "<subject>_Rwhite_curv_blobs", SetType( 'Blob White Curvature Texture' ), SetWeakAttr( 'side', 'right' ),
    "<subject>_Lwhite_curv_ss", SetType( 'Scale Space White Curvature Texture' ), SetWeakAttr( 'side', 'left' ),
    "<subject>_Rwhite_curv_ss", SetType( 'Scale Space White Curvature Texture' ), SetWeakAttr( 'side', 'right' ),
  ),
)

insert( '{protocol}/{subject}',
  'sulci', SetContent(
     '<subject>_L_{sulcus_name}', SetType( 'Sulcus mesh' ), SetWeakAttr( 'side', 'left' ),
     '<subject>_R_{sulcus_name}', SetType( 'Sulcus mesh' ), SetWeakAttr( 'side', 'right' ),
     '<subject>_L_{sulcus_name}-x', SetType( 'Sulcus x coordinate texture' ), SetWeakAttr( 'side', 'left' ),
     '<subject>_R_{sulcus_name}-x', SetType( 'Sulcus x coordinate texture' ), SetWeakAttr( 'side', 'right' ),
     '<subject>_L_{sulcus_name}-y', SetType( 'Sulcus y coordinate texture' ), SetWeakAttr( 'side', 'left' ),
     '<subject>_R_{sulcus_name}-y', SetType( 'Sulcus y coordinate texture' ), SetWeakAttr( 'side', 'right' ),
     '<subject>_L_{sulcus_name}-depth', SetType( 'Sulcus depth profile' ), SetWeakAttr( 'side', 'left' ),
     '<subject>_R_{sulcus_name}-depth', SetType( 'Sulcus depth profile' ), SetWeakAttr( 'side', 'right' ),
     '<subject>_L_{sulcus_name}-grid', SetType( 'Sulcus coordinate grid mesh' ), SetWeakAttr( 'side', 'left'),
     '<subject>_R_{sulcus_name}-grid', SetType( 'Sulcus coordinate grid mesh' ), SetWeakAttr( 'side', 'right')
  ),
)



