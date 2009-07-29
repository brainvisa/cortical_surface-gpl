include( 'base' )
include( 'sulci' )

insert( '{protocol}/{subject}',
  'surface', SetContent(
      ## que mettre dans les niveaux aquisition et analysis?
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
    '<subject>_Lhippo_Volume', SetType( 'Left Cingular Pole Template Subject' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rhippo_Volume', SetType( 'Right Cingular Pole Template Subject' ), SetWeakAttr( 'side', 'right' ),
    '<subject>_Talairach_To_Subject_Transformation', SetType( 'Talairach To Subject Transformation'),
    '<subject>_Subject_To_Template_Transformation', SetType( 'Subject To Template Transformation' ),
    #'<subject>_L_gyriGraph', SetType( 'Left Gyri Graph' ), SetWeakAttr( 'side', 'left' ),
    #'<subject>_R_gyriGraph', SetType( 'Right Gyri Graph' ), SetWeakAttr( 'side', 'right' ),
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
    ## Arnaud -> analyse de surface et structurelle -> a verifier
    ## utilise dans le traitement 2DGeodesicPrimalSketch de la toolbox cortical_surface
    "<subject>_Rwhite_primal",SetType( 'Primal Sketch' ),SetWeakAttr( 'side', 'right' ),
    "<subject>_Lwhite_primal",SetType( 'Primal Sketch' ),SetWeakAttr( 'side', 'left' ),
    "<subject>_Rwhite_GLB",SetType( 'Grey Level Blob Graph' ),SetWeakAttr( 'side', 'right' ),
    "<subject>_Lwhite_GLB",SetType( 'Grey Level Blob Graph' ),SetWeakAttr( 'side', 'left' ),
    "<subject>_Lwhite_curv_blob", SetType( 'Blob White Curvature Texture' ), SetWeakAttr( 'side', 'left' ),
    "<subject>_Rwhite_curv_blob", SetType( 'Blob White Curvature Texture' ), SetWeakAttr( 'side', 'right' ),
    "<subject>_Lwhite_curv_ss", SetType( 'Scale Space White Curvature Texture' ), SetWeakAttr( 'side', 'left' ),
    "<subject>_Rwhite_curv_ss", SetType( 'Scale Space White Curvature Texture' ), SetWeakAttr( 'side', 'right' ),
  ),
)

# Add To white texture translation and Gyri Graph under sulci_recognition_session
def _insertWhiteTextureGraph( recognition_type ):
  return ('{protocol}/{subject}/t1mri/{acquisition}/{analysis}/folds/{graph_version}/{sulci_recognition_session}_' + recognition_type,
      "<subject>_left_sulci_to_texture_<sulci_recognition_session>_" + recognition_type, SetType( 'Sulci To White Texture Translation' ),SetWeakAttr( 'side', 'left' ),
      "<subject>_right_sulci_to_texture_<sulci_recognition_session>_" + recognition_type, SetType( 'Sulci To White Texture Translation' ),SetWeakAttr( 'side', 'right' ),
      # GYRI - graphs, to white texture translation
      "<subject>_Rgyri_<sulci_recognition_session>_" + recognition_type, SetType( 'Gyri Graph' ),SetWeakAttr( 'side', 'right' ),
      "<subject>_Lgyri_<sulci_recognition_session>_" + recognition_type, SetType( 'Gyri Graph' ),SetWeakAttr( 'side', 'left' ),
      "<subject>_left_gyri_to_texture_<sulci_recognition_session>_" + recognition_type, SetType( 'Gyri To White Texture Translation' ),SetWeakAttr( 'side', 'left' ),
      "<subject>_right_gyri_to_texture_<sulci_recognition_session>_" + recognition_type, SetType( 'Gyri To White Texture Translation' ),SetWeakAttr( 'side', 'right' )
    )
 
# Add Sulci/Gyri White Texture and Volume in sulci_recognition_session/segmentation
def _insertWhiteTexture( recognition_type ):
  return ('<protocol>/<subject>/t1mri/<acquisition>/<analysis>/folds/<graph_version>/<sulci_recognition_session>_' + recognition_type + '/segmentation',
      "<subject>_Lwhite_sulci_<sulci_recognition_session>_" + recognition_type, SetType( 'Sulci White Texture' ), SetWeakAttr( 'side', 'left' ),
      "<subject>_Rwhite_sulci_<sulci_recognition_session>_" + recognition_type, SetType( 'Sulci White Texture' ), SetWeakAttr( 'side', 'right' ),
      "<subject>_Lwhite_gyri_<sulci_recognition_session>_" + recognition_type, SetType( 'Gyri White Texture' ), SetWeakAttr( 'side', 'left' ),
      "<subject>_Rwhite_gyri_<sulci_recognition_session>_" + recognition_type, SetType( 'Gyri White Texture' ), SetWeakAttr( 'side', 'right' ),
      "<subject>_Lwhite_gyri_<sulci_recognition_session>_" + recognition_type, SetType( 'Gyri White Volume' ), SetWeakAttr( 'side', 'left' ),
      "<subject>_Rwhite_gyri_<sulci_recognition_session>_" + recognition_type, SetType( 'Gyri White Volume' ), SetWeakAttr( 'side', 'right' )
      )

insert( *_insertWhiteTextureGraph( 'auto' ) )
insert( *_insertWhiteTexture( 'auto' ) )
insert( *_insertWhiteTextureGraph( 'manual' ) )
insert( *_insertWhiteTexture( 'manual' ) )
insert( *_insertWhiteTextureGraph( 'best' ) )
insert( *_insertWhiteTexture( 'best' ) )
del _insertWhiteTexture, _insertWhiteTextureGraph

insert( '{protocol}/{subject}',
  #'sulci', SetContent(
    #'{acquisition}', SetDefaultAttributeValue( 'acquisition', default_acquisition ), SetNonMandatoryKeyAttribute( 'acquisition' ),
  'sulci', SetContent(
          '<subject>_L_{sulcus_name}', SetType( 'Sulcus mesh' ), SetWeakAttr( 'side', 'left' ),
          '<subject>_R_{sulcus_name}', SetType( 'Sulcus mesh' ), SetWeakAttr( 'side', 'right' ),
          '<subject>_L_{sulcus_name}-x', SetType( 'Sulcus x coordinate texture' ), SetWeakAttr( 'side', 'left' ),
          '<subject>_R_{sulcus_name}-x', SetType( 'Sulcus x coordinate texture' ), SetWeakAttr( 'side', 'right' ),
          '<subject>_L_{sulcus_name}-y', SetType( 'Sulcus y coordinate texture' ), SetWeakAttr( 'side', 'left' ),
          '<subject>_R_{sulcus_name}-y', SetType( 'Sulcus y coordinate texture' ), SetWeakAttr( 'side', 'right' ),
          '<subject>_L_{sulcus_name}-depth', SetType( 'Sulcus depth profile' ), SetWeakAttr( 'side', 'left' ),
          '<subject>_R_{sulcus_name}-depth', SetType( 'Sulcus depth profile' ), SetWeakAttr( 'side', 'right' ),
          '<subject>_L_{sulcus_name}-grid.mesh', SetType( 'Sulcus coordinate grid mesh' ), SetWeakAttr( 'side', 'left'),
          '<subject>_R_{sulcus_name}-grid.mesh', SetType( 'Sulcus coordinate grid mesh' ), SetWeakAttr( 'side', 'right')
     ),
  #),
)
