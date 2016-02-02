include( 'base' )
include( 'sulci' )

insert( '{center}/{subject}/t1mri/{acquisition}/{analysis}/segmentation/mesh',
    "<subject>_Lwhite_resampled", SetType("Resampled Hemisphere White Mesh"),
        SetWeakAttr('side', 'left'),
    "<subject>_Rwhite_resampled", SetType("Resampled Hemisphere White Mesh"),
        SetWeakAttr('side', 'right'),
    "<subject>_Bwhite_resampled", SetType("Resampled Hemisphere White Mesh"),
        SetWeakAttr('side', 'both'),

    "<subject>_Lwhite_sphere", SetType("Spherical Mesh"),
        SetWeakAttr('side', 'left'),
    "<subject>_Rwhite_sphere", SetType("Spherical Mesh"),
        SetWeakAttr('side', 'right'),
)


insert( '{center}/{subject}/t1mri/{acquisition}/{analysis}/segmentation/mesh',
  'surface_analysis', SetContent(
    ## que mettre dans les niveaux aquisition et analysis?

    ## sulci voronoi
    '<subject>_Lwhite_sulci_voronoi', SetType( 'Left sulci voronoi texture' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rwhite_sulci_voronoi', SetType( 'Right sulci voronoi texture' ), SetWeakAttr( 'side', 'right' ),

    ## Laplacian Eigen vectors
    '<subject>_Lwhite_eigen_vectors', SetType( 'Left eigen vectors texture' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rwhite_eigen_vectors', SetType( 'Right eigen vectors texture' ), SetWeakAttr( 'side', 'right' ),

    ## sulcal pits
    '<subject>_Lwhite_pits_labeled', SetType( 'Left labeled pits texture' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rwhite_pits_labeled', SetType( 'Right labeled pits texture' ), SetWeakAttr( 'side', 'right' ),
    '<subject>_Lwhite_basins_labeled', SetType( 'Left labeled basins texture' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rwhite_basins_labeled', SetType( 'Right labeled basins texture' ), SetWeakAttr( 'side', 'right' ),

    '<subject>_Lwhite_pits', SetType( 'Left pits texture' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rwhite_pits', SetType( 'Right pits texture' ), SetWeakAttr( 'side', 'right' ),
    '<subject>_Lwhite_noisy_pits', SetType( 'Left noisy pits texture' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rwhite_noisy_pits', SetType( 'Right noisy pits texture' ), SetWeakAttr( 'side', 'right' ),
    '<subject>_Lwhite_ridges', SetType( 'Left ridges texture' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rwhite_ridges', SetType( 'Right ridges texture' ), SetWeakAttr( 'side', 'right' ),
    '<subject>_Lwhite_basins', SetType( 'Left basins texture' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rwhite_basins', SetType( 'Right basins texture' ), SetWeakAttr( 'side', 'right' ),

    ## Depth Potential Function
    '<subject>_Lwhite_DPF', SetType( 'Left DPF texture' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rwhite_DPF', SetType( 'Right DPF texture' ), SetWeakAttr( 'side', 'right' ),

    ## HIP-HOP parameterization
    '<subject>_Lwhite_remeshed_hiphop', SetType( 'Left remeshed mesh' ),
      SetWeakAttr('side', 'left', 'averaged', 'No', 'inflated', 'No',
                  'vertex_corr', 'Yes', 'vertex_corr_method', 'hiphop'),
    '<subject>_Rwhite_remeshed_hiphop', SetType( 'Right remeshed mesh' ),
      SetWeakAttr('side', 'right', 'averaged', 'No', 'inflated', 'No',
                  'vertex_corr', 'Yes', 'vertex_corr_method', 'hiphop'),
    '<subject>_Bwhite_remeshed_hiphop', SetType( 'Remeshed mesh' ),
      SetWeakAttr('side', 'both', 'averaged', 'No', 'inflated', 'No',
                  'vertex_corr', 'Yes', 'vertex_corr_method', 'hiphop'),
    '<subject>_Lwhite_remeshed_inflated_hiphop',
      SetType( 'Left remeshed mesh' ), SetPriorityOffset(-1),
      SetWeakAttr('side', 'left', 'averaged', 'No', 'inflated', 'Yes',
                  'vertex_corr', 'Yes', 'vertex_corr_method', 'hiphop'),
    '<subject>_Rwhite_remeshed_inflated_hiphop',
      SetType( 'Right remeshed mesh' ), SetPriorityOffset(-1),
      SetWeakAttr('side', 'right', 'averaged', 'No', 'inflated', 'Yes',
                  'vertex_corr', 'Yes', 'vertex_corr_method', 'hiphop'),
    '<subject>_Bwhite_remeshed_inflated_hiphop',
      SetType( 'Remeshed mesh' ), SetPriorityOffset(-1),
      SetWeakAttr('side', 'both', 'averaged', 'No', 'inflated', 'Yes',
                  'vertex_corr', 'Yes', 'vertex_corr_method', 'hiphop'),
    '<subject>_rectangular_flat', SetType( 'Rectangular flat texture' ),
    '<subject>_Lwhite_sulcalines_rectangular_flat', SetType( 'Left hemisphere Sulcal Lines Rectangular Flat texture' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rwhite_sulcalines_rectangular_flat', SetType( 'Right hemisphere Sulcal Lines Rectangular Flat texture' ), SetWeakAttr( 'side', 'right' ),
    '<subject>_Lwhite_rectangular_flat', SetType( 'Left rectangular flat mesh' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rwhite_rectangular_flat', SetType( 'Right rectangular flat mesh' ), SetWeakAttr( 'side', 'right' ),
    '<subject>_Lwhite_rectangular_flat_cstr', SetType( 'Left rectangular flat cstr mesh' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rwhite_rectangular_flat_cstr', SetType( 'Right rectangular flat cstr mesh' ), SetWeakAttr( 'side', 'right' ),
    '<subject>_Lwhite_rectangular_flat_indices', SetType( 'Left rectangular flat indices texture' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rwhite_rectangular_flat_indices', SetType( 'Right rectangular flat indices texture' ), SetWeakAttr( 'side', 'right' ),
    '<subject>_Lwhite_rectangular_flat_boundary', SetType( 'Left rectangular boundary texture' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rwhite_rectangular_flat_boundary', SetType( 'Right rectangular boundary texture' ), SetWeakAttr( 'side', 'right' ),
    '<subject>_Lwhite_parts', SetType( 'Left white mesh parts' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rwhite_parts', SetType( 'Right white mesh parts' ), SetWeakAttr( 'side', 'right' ),
    '<subject>_Lwhite_spherical', SetType( 'Left hemisphere spherical mesh' ),
      SetWeakAttr('side', 'left', 'averaged', 'No', 'inflated', 'Yes',
                  'vertex_corr', 'Yes', 'vertex_corr_method', 'hiphop'),
    '<subject>_Rwhite_spherical', SetType( 'Right hemisphere spherical mesh' ),
      SetWeakAttr('side', 'right', 'averaged', 'No', 'inflated', 'Yes',
                  'vertex_corr', 'Yes', 'vertex_corr_method', 'hiphop'),
    '<subject>_Lwhite_pole_insula', SetType( 'Left insula pole texture' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rwhite_pole_insula', SetType( 'Right insula pole texture' ), SetWeakAttr( 'side', 'right' ),
    '<subject>_Lwhite_pole_cingular', SetType( 'Left cingular pole texture' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rwhite_pole_cingular', SetType( 'Right cingular pole texture' ), SetWeakAttr( 'side', 'right' ),
    '<subject>_Lwhite_poles', SetType( 'Left poles texture' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rwhite_poles', SetType( 'Right poles texture' ), SetWeakAttr( 'side', 'right' ),
    '<subject>_Lwhite_lat_cleaned', SetType( 'Left hemisphere latitude cleaned constraints texture'), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rwhite_lat_cleaned', SetType( 'Right hemisphere latitude cleaned constraints texture'), SetWeakAttr( 'side', 'right' ),
    '<subject>_Lwhite_lon_cleaned', SetType( 'Left hemisphere longitude cleaned constraints texture'), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rwhite_lon_cleaned', SetType( 'Right hemisphere longitude cleaned constraints texture'), SetWeakAttr( 'side', 'right' ),
    '<subject>_Lwhite_lat_roots', SetType( 'Left hemisphere latitude constraints texture'), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rwhite_lat_roots', SetType( 'Right hemisphere latitude constraints texture'), SetWeakAttr( 'side', 'right' ),
    '<subject>_Lwhite_lon_roots', SetType( 'Left hemisphere longitude constraints texture'), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rwhite_lon_roots', SetType( 'Right hemisphere longitude constraints texture'), SetWeakAttr( 'side', 'right' ),
    '<subject>_Lwhite_lat', SetType( 'Left hemisphere latitude texture'),
      SetWeakAttr('side', 'left', 'averaged', 'No', 'vertex_corr', 'No'),
    '<subject>_Rwhite_lat', SetType( 'Right hemisphere latitude texture'),
      SetWeakAttr('side', 'right', 'averaged', 'No', 'vertex_corr', 'No'),
    '<subject>_Lwhite_lon', SetType( 'Left hemisphere longitude texture'),
      SetWeakAttr('side', 'left', 'averaged', 'No', 'vertex_corr', 'No'),
    '<subject>_Rwhite_lon', SetType( 'Right hemisphere longitude texture'),
      SetWeakAttr('side', 'right', 'averaged', 'No', 'vertex_corr', 'No'),
    '<subject>_Lwhite_lat_HIP', SetType( 'Left hemisphere latitude HIP texture'), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rwhite_lat_HIP', SetType( 'Right hemisphere latitude HIP texture'), SetWeakAttr( 'side', 'right' ),
    '<subject>_Lwhite_lon_HIP', SetType( 'Left hemisphere longitude HIP texture'), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rwhite_lon_HIP', SetType( 'Right hemisphere longitude HIP texture'), SetWeakAttr( 'side', 'right' ),
    '<subject>_Lconf_lat', SetType( 'Conformal latitude texture'), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rconf_lat', SetType( 'Conformal latitude texture'), SetWeakAttr( 'side', 'right' ),
    '<subject>_Lconf_lon', SetType( 'Conformal longitude texture'), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rconf_lon', SetType( 'Conformal longitude texture'), SetWeakAttr( 'side', 'right' ),
    '<subject>_rootsValues', SetType( 'Constraint coordinates values'),
    '<subject>_Lwhite_grid', SetType( 'Left hemisphere coordinate grid'), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rwhite_grid', SetType( 'Right hemisphere coordinate grid'), SetWeakAttr( 'side', 'right' ),
    '<subject>_Lhippo_Volume', SetType( 'Left Cingular Pole Template Subject' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rhippo_Volume', SetType( 'Right Cingular Pole Template Subject' ), SetWeakAttr( 'side', 'right' ),
    '<subject>_Talairach_To_Subject_Transformation', SetType( 'Talairach To Subject Transformation'),
    '<subject>_Subject_To_Template_Transformation', SetType( 'Subject To Template Transformation' ),

#    ## Parcels from HIP-HOP parameterization
    '<subject>_Lwhite_parcels_model', SetType( 'Left hemisphere model parcellation texture'),
      SetWeakAttr('side', 'left', 'averaged', 'No', 'vertex_corr', 'No',
                  'parcellation_type', 'model', 'regularized', 'false'),
    '<subject>_Rwhite_parcels_model', SetType( 'Right hemisphere model parcellation texture'),
      SetWeakAttr('side', 'right', 'averaged', 'No', 'vertex_corr', 'No',
                  'parcellation_type', 'model', 'regularized', 'false'),
    '<subject>_Lwhite_parcels_marsAtlas', SetType( 'Left hemisphere marsAtlas parcellation texture'),
      SetWeakAttr('side', 'left', 'averaged', 'No', 'vertex_corr', 'No',
                  'parcellation_type', 'marsAtlas', 'regularized', 'false'),
    '<subject>_Rwhite_parcels_marsAtlas', SetType( 'Right hemisphere marsAtlas parcellation texture'),
      SetWeakAttr('side', 'right', 'averaged', 'No', 'vertex_corr', 'No',
                  'parcellation_type', 'marsAtlas', 'regularized', 'false'),
  ),
)

insert( '{center}/{subject}/t1mri/{acquisition}/{analysis}/segmentation/mesh',
  'surface_analysis', SetContent(

    ## same with resampling to sphere coordinates
    '<subject>_Lwhite_parcels_model_remeshed_{vertex_corr_method}',
      SetType('Left hemisphere model parcellation texture'),
      SetWeakAttr('side', 'left', 'averaged', 'No', 'vertex_corr', 'Yes',
                  'parcellation_type', 'model'),
    '<subject>_Rwhite_parcels_model_remeshed_{vertex_corr_method}',
      SetType('Right hemisphere model parcellation texture'),
      SetWeakAttr('side', 'right', 'averaged', 'No', 'vertex_corr', 'Yes',
                  'parcellation_type', 'model'),
    '<subject>_Bwhite_parcels_model_remeshed_{vertex_corr_method}',
      SetType('Hemisphere model parcellation texture'),
      SetWeakAttr('side', 'both', 'averaged', 'No', 'vertex_corr', 'Yes',
                  'parcellation_type', 'model'),
    '<subject>_Lwhite_parcels_marsAtlas_remeshed_{vertex_corr_method}',
      SetType('Left hemisphere marsAtlas parcellation texture'),
      SetWeakAttr('side', 'left', 'averaged', 'No', 'vertex_corr', 'Yes',
                  'parcellation_type', 'marsAtlas'),
    '<subject>_Rwhite_parcels_marsAtlas_remeshed_{vertex_corr_method}',
      SetType('Right hemisphere marsAtlas parcellation texture'),
      SetWeakAttr('side', 'right', 'averaged', 'No', 'vertex_corr', 'Yes',
                  'parcellation_type', 'marsAtlas'),
    '<subject>_Bwhite_parcels_marsAtlas_remeshed_{vertex_corr_method}',
      SetType('Hemisphere marsAtlas parcellation texture'),
      SetWeakAttr('side', 'both', 'averaged', 'No', 'vertex_corr', 'Yes',
                  'parcellation_type', 'marsAtlas'),
#
## may be used to replace gyri regularized parcellation texture
##    '<subject>_Lwhite_model_parcels_regul', SetType( 'Left hemisphere regularized model parcellation texture'), SetWeakAttr( 'side', 'left' ), SetWeakAttr( 'regularized', 'true' ), SetWeakAttr( 'parcels_type', 'model' ),
##    '<subject>_Rwhite_model_parcels_regul', SetType( 'Right hemisphere regularized model parcellation texture'), SetWeakAttr( 'side', 'right' ), SetWeakAttr( 'regularized', 'true' ), SetWeakAttr( 'parcels_type', 'model' ),
#
#    '<subject>_L_model_parcelsVolume', SetType( 'Left model parcellation volume' ), SetWeakAttr( 'side', 'left' ),
#    '<subject>_R_model_parcelsVolume', SetType( 'Right model parcellation volume' ), SetWeakAttr( 'side', 'right' ),
#    '<subject>_L_coarse_parcelsVolume', SetType( 'Left coarse parcellation volume' ), SetWeakAttr( 'side', 'left' ),
#    '<subject>_R_coarse_parcelsVolume', SetType( 'Right coarse parcellation volume' ), SetWeakAttr( 'side', 'right' ),
# Parcels Graphs do not work....
#    '<subject>_L_parcelsGraph', SetType( 'Left Parcels Graph' ), SetWeakAttr( 'side', 'left' ),
#    '<subject>_R_parcelsGraph', SetType( 'Right Parcels Graph' ), SetWeakAttr( 'side', 'right' ),

    ## Old Gyri parcellation types, should be replaced by Parcels (above)
    '<subject>_Lwhite_gyri', SetType( 'Left hemisphere gyri parcellation texture'), SetWeakAttr( 'side', 'left' ), SetWeakAttr( 'regularized', 'false' ),
    '<subject>_Rwhite_gyri', SetType( 'Right hemisphere gyri parcellation texture'), SetWeakAttr( 'side', 'right' ), SetWeakAttr( 'regularized', 'false' ),
    '<subject>_Lwhite_gyri_regul', SetType( 'Left hemisphere regularized parcellation texture'), SetWeakAttr( 'side', 'left' ), SetWeakAttr( 'regularized', 'true' ),
    '<subject>_Rwhite_gyri_regul', SetType( 'Right hemisphere regularized parcellation texture'), SetWeakAttr( 'side', 'right' ), SetWeakAttr( 'regularized', 'true' ),
    
    ## volume parcellation s
    '<subject>_parcellation', SetType( 'Parcellation volume' ),
  ),
)


insert( '{center}/{subject}/t1mri/{acquisition}/{analysis}/segmentation/mesh/surface_analysis', 
    ## utilise dans le traitement 2DGeodesicPrimalSketch de la toolbox cortical_surface    
    '<subject>_L_gyriGraph', SetType( 'Left Gyri Graph' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_R_gyriGraph', SetType( 'Right Gyri Graph' ), SetWeakAttr( 'side', 'right' ),
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
)

# Cortical Surface Functional-related Types

insert( '{center}/{subject}/t1mri/{acquisition}/{analysis}/segmentation/mesh/surface_analysis',
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
