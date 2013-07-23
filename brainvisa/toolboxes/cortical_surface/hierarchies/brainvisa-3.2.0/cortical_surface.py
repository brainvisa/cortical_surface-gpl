include( 'base' )
include( 'sulci' )

insert( '{center}/{subject}',
  'surface', SetContent(
      ## que mettre dans les niveaux aquisition et analysis?
    '<subject>_L_remeshed', SetType( 'Left remeshed mesh' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_R_remeshed', SetType( 'Right remeshed mesh' ), SetWeakAttr( 'side', 'right' ),
    '<subject>_rectangular_flat', SetType( 'Rectangular flat texture' ),
    '<subject>_Lwhite_sulcalines_rectangular_flat', SetType( 'Left hemisphere Sulcal Lines Rectangular Flat texture' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rwhite_sulcalines_rectangular_flat', SetType( 'Right hemisphere Sulcal Lines Rectangular Flat texture' ), SetWeakAttr( 'side', 'right' ),
    '<subject>_L_rectangular_flat', SetType( 'Left rectangular flat mesh' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_R_rectangular_flat', SetType( 'Right rectangular flat mesh' ), SetWeakAttr( 'side', 'right' ),
    '<subject>_L_rectangular_flat_cstr', SetType( 'Left rectangular flat cstr mesh' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_R_rectangular_flat_cstr', SetType( 'Right rectangular flat cstr mesh' ), SetWeakAttr( 'side', 'right' ),
    '<subject>_L_rectangular_flat_indices', SetType( 'Left rectangular flat indices texture' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_R_rectangular_flat_indices', SetType( 'Right rectangular flat indices texture' ), SetWeakAttr( 'side', 'right' ),
    '<subject>_L_rectangular_flat_boundary', SetType( 'Left rectangular boundary texture' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_R_rectangular_flat_boundary', SetType( 'Right rectangular boundary texture' ), SetWeakAttr( 'side', 'right' ),
    '<subject>_L_spherical', SetType( 'Left spherical mesh' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_R_spherical', SetType( 'Right spherical mesh' ), SetWeakAttr( 'side', 'right' ),
    '<subject>_Linsula', SetType( 'Left insula pole texture' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rinsula', SetType( 'Right insula pole texture' ), SetWeakAttr( 'side', 'right' ),
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
    '<subject>_Lconf_lat', SetType( 'Conformal latitude texture'), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rconf_lat', SetType( 'Conformal latitude texture'), SetWeakAttr( 'side', 'right' ),
    '<subject>_Lconf_lon', SetType( 'Conformal longitude texture'), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rconf_lon', SetType( 'Conformal longitude texture'), SetWeakAttr( 'side', 'right' ),
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
    
    #'<subject>_L_gyriGraph', SetType( 'Left Gyri Graph' ), SetWeakAttr( 'side', 'left' ),
    #'<subject>_R_gyriGraph', SetType( 'Right Gyri Graph' ), SetWeakAttr( 'side', 'right' ),
    '<subject>_L_gyriVolume', SetType( 'Left Gyri Volume' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_R_gyriVolume', SetType( 'Right Gyri Volume' ), SetWeakAttr( 'side', 'right' ),
    '<subject>_Lwhite_KERNEL', SetType( 'Projection convolution kernels'), SetWeakAttr( 'side', 'left' ) ,
    '<subject>_Rwhite_KERNEL', SetType( 'Projection convolution kernels'), SetWeakAttr( 'side', 'right' ) ,

    '<subject>_Lwhite_thickness', SetType('Cortical thickness'), SetWeakAttr('side', 'left', 'mesh', 'white'),
    '<subject>_Rwhite_thickness', SetType('Cortical thickness'), SetWeakAttr('side', 'right', 'mesh', 'white'),
    '<subject>_Lhemi_thickness', SetType('Cortical thickness'), SetWeakAttr('side', 'left', 'mesh', 'hemi'),
    '<subject>_Rhemi_thickness', SetType('Cortical thickness'), SetWeakAttr('side', 'right', 'mesh', 'hemi'), 
        
    ## Arnaud -> analyse de surface et structurelle -> a verifier
    ## utilise dans le traitement 2DGeodesicPrimalSketch de la toolbox cortical_surface
    "<subject>_Rwhite_primal",SetType( 'Primal Sketch' ),SetWeakAttr( 'side', 'right' ),
    "<subject>_Lwhite_primal",SetType( 'Primal Sketch' ),SetWeakAttr( 'side', 'left' ),
    "<subject>_Lwhite_primal_curv",SetType( 'Curvature Map Primal Sketch' ),SetWeakAttr( 'side', 'left' ),
    "<subject>_Rwhite_primal_curv",SetType( 'Curvature Map Primal Sketch' ),SetWeakAttr( 'side', 'right' ),
    "<subject>_Lwhite_primal_depth",SetType( 'Depth Map Primal Sketch' ),SetWeakAttr( 'side', 'left' ),
    "<subject>_Rwhite_primal_depth",SetType( 'Depth Map Primal Sketch' ),SetWeakAttr( 'side', 'right' ),
    "<subject>_Lwhite_curv_blobs", SetType( 'Blob White Curvature Texture' ), SetWeakAttr( 'side', 'left' ),
    "<subject>_Rwhite_curv_blobs", SetType( 'Blob White Curvature Texture' ), SetWeakAttr( 'side', 'right' ),
    "<subject>_Lwhite_curv_ss", SetType( 'Scale Space White Curvature Texture' ), SetWeakAttr( 'side', 'left' ),
    "<subject>_Rwhite_curv_ss", SetType( 'Scale Space White Curvature Texture' ), SetWeakAttr( 'side', 'right' ),
    "<subject>_Lwhite_depth_blobs", SetType( 'Blob White Depth Texture' ), SetWeakAttr( 'side', 'left' ),
    "<subject>_Rwhite_depth_blobs", SetType( 'Blob White Depth Texture' ), SetWeakAttr( 'side', 'right' ),
    "<subject>_Lwhite_depth_ss", SetType( 'Scale Space White Depth Texture' ), SetWeakAttr( 'side', 'left' ),
    "<subject>_Rwhite_depth_ss", SetType( 'Scale Space White Depth Texture' ), SetWeakAttr( 'side', 'right' ),  
    
    "<subject>_Lwhite_sulcalines", SetType( 'Left hemisphere Sulcal Lines texture' ), SetWeakAttr( 'side', 'left' ),
    "<subject>_Rwhite_sulcalines", SetType( 'Right hemisphere Sulcal Lines texture' ), SetWeakAttr( 'side', 'right' ),
    "<subject>_LgraphLabelBasins", SetType( 'Left Graph Label Translation' ), SetWeakAttr( 'side', 'left' ),
    "<subject>_RgraphLabelBasins", SetType( 'Right Graph Label Translation' ), SetWeakAttr( 'side', 'right' ),
  ),
)

# Cortical Surface Functional-related Types

insert( '{center}/{subject}/surface',
  'functional', SetDefaultAttributeValue( 'contrast', 'con' ), SetDefaultAttributeValue( 'volume', 'vol' ), SetContent(
      '<subject>_Anatomy_To_Mean_Function_Transformation', SetType( 'Anatomy To Mean Functional Volume Transformation' ),
      '<subject>_Mean_Function_To_Anatomy_Transformation', SetType( 'Mean Functional Volume To Anatomy Transformation' ),
    "<subject>_L_labels", SetType( 'Labelled Functional Blobs Texture' ), SetWeakAttr( 'side', 'left' ),
    "<subject>_R_labels", SetType( 'Labelled Functional Blobs Texture' ), SetWeakAttr( 'side', 'right' ),
    '<subject>_{volume}_Lwhite_proj', SetType( 'Functional Texture'), SetWeakAttr( 'side', 'left' ) ,
    '<subject>_{volume}_Rwhite_proj', SetType( 'Functional Texture'), SetWeakAttr( 'side', 'right' ),
    '<subject>_{volume}_Lwhite', SetType( 'Functional Time Texture'), SetWeakAttr( 'side', 'left' ) ,
    '<subject>_{volume}_Rwhite', SetType( 'Functional Time Texture'), SetWeakAttr( 'side', 'right' ),
    '<subject>_Lwhite_spmT_{contrast}', SetType( 'Surface-Based SPMt Map'), SetWeakAttr( 'side', 'left' ) ,
    '<subject>_Rwhite_spmT_{contrast}', SetType( 'Surface-Based SPMt Map'), SetWeakAttr( 'side', 'right' ),
    '<subject>_Lwhite_beta_{contrast}', SetType( 'Surface-Based Beta Map'), SetWeakAttr( 'side', 'left' ) ,
    '<subject>_Rwhite_beta_{contrast}', SetType( 'Surface-Based Beta Map'), SetWeakAttr( 'side', 'right' ),
   ),
)

    
    
# Add To white texture translation and Gyri Graph under sulci_recognition_session
def _insertWhiteTextureGraph( recognition_type ):
  return ('{center}/{subject}/t1mri/{acquisition}/{analysis}/folds/{graph_version}/{sulci_recognition_session}_' + recognition_type,
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
  return ('{center}/{subject}/t1mri/{acquisition}/{analysis}/folds/{graph_version}/{sulci_recognition_session}_' + recognition_type,
    'segmentation', SetContent(
      "<subject>_Lwhite_sulci_<sulci_recognition_session>_" + recognition_type, SetType( 'Sulci White Texture' ), SetWeakAttr( 'side', 'left' ),
      "<subject>_Rwhite_sulci_<sulci_recognition_session>_" + recognition_type, SetType( 'Sulci White Texture' ), SetWeakAttr( 'side', 'right' ),
      "<subject>_Lwhite_gyri_<sulci_recognition_session>_" + recognition_type, SetType( 'Gyri White Texture' ), SetWeakAttr( 'side', 'left' ),
      "<subject>_Rwhite_gyri_<sulci_recognition_session>_" + recognition_type, SetType( 'Gyri White Texture' ), SetWeakAttr( 'side', 'right' ),
      "<subject>_Lwhite_gyri_<sulci_recognition_session>_" + recognition_type, SetType( 'Gyri White Volume' ), SetWeakAttr( 'side', 'left' ),
      "<subject>_Rwhite_gyri_<sulci_recognition_session>_" + recognition_type, SetType( 'Gyri White Volume' ), SetWeakAttr( 'side', 'right' )
      ),
  )


insert( *_insertWhiteTextureGraph( 'auto' ) )
insert( *_insertWhiteTexture( 'auto' ) )
insert( *_insertWhiteTextureGraph( 'manual' ) )
insert( *_insertWhiteTexture( 'manual' ) )
insert( *_insertWhiteTextureGraph( 'best' ) )
insert( *_insertWhiteTexture( 'best' ) )
del _insertWhiteTexture, _insertWhiteTextureGraph

insert( '{center}/{subject}',
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
