include( 'builtin' )
include( 'registration' )
include( 'anatomy' )


## surface coordinates
FileType( 'Hippocampus pole texture', 'Texture' )
FileType( 'Left hippocampus pole texture', 'Hippocampus pole texture' )
FileType( 'Right hippocampus pole texture', 'Hippocampus pole texture' )
FileType( 'Poles texture', 'Texture' )
FileType( 'Left poles texture', 'Poles texture' )
FileType( 'Right poles texture', 'Poles texture' )

FileType( 'Constraints texture', 'Label Texture' )
FileType( 'Latitude constraints texture', 'Constraints texture' )
FileType( 'Left hemisphere latitude constraints texture', 'Latitude constraints texture' )
FileType( 'Right hemisphere latitude constraints texture', 'Latitude constraints texture' )
FileType( 'Left hemisphere latitude cleaned constraints texture', 'Latitude constraints texture' )
FileType( 'Right hemisphere latitude cleaned constraints texture', 'Latitude constraints texture' )
FileType( 'Longitude constraints texture', 'Constraints texture' )
FileType( 'Left hemisphere longitude constraints texture', 'Longitude constraints texture' )
FileType( 'Right hemisphere longitude constraints texture', 'Longitude constraints texture' )
FileType( 'Left hemisphere longitude cleaned constraints texture', 'Longitude constraints texture' )
FileType( 'Right hemisphere longitude cleaned constraints texture', 'Longitude constraints texture' )

FileType( 'Constraint coordinates values', 'Text file')

FileType( 'Coordinate texture', 'Texture' )
FileType( 'Latitude coordinate texture', 'Coordinate texture' )
FileType( 'Longitude coordinate texture', 'Coordinate texture' )
FileType( 'Left hemisphere longitude texture', 'Longitude coordinate texture' )
FileType( 'Left hemisphere latitude texture', 'Latitude coordinate texture' )
FileType( 'Right hemisphere longitude texture', 'Longitude coordinate texture' )
FileType( 'Right hemisphere latitude texture', 'Latitude coordinate texture' )

FileType( 'Hemisphere gyri parcellation texture', 'Gyri White Texture' )
FileType( 'Left hemisphere gyri parcellation texture', 'Hemisphere gyri parcellation texture' )
FileType( 'Right hemisphere gyri parcellation texture', 'Hemisphere gyri parcellation texture' )

FileType( 'Coordinate grid', 'Mesh' )
FileType( 'Left hemisphere coordinate grid', 'Coordinate grid')
FileType( 'Right hemisphere coordinate grid', 'Coordinate grid')

## template for poles
FileType( 'Left Cingular Pole Template', '3D Volume' )
FileType( 'Right Cingular Pole Template', '3D Volume' )

FileType( 'Left Cingular Pole Template Subject', '3D Volume' )
FileType( 'Right Cingular Pole Template Subject', '3D Volume' )

## model files

FileType( 'Surface Label Translation', 'Label Translation')
FileType( 'Latitude Constraint Gyri Model', 'Gyri Model' )
FileType( 'Longitude Constraint Gyri Model', 'Gyri Model' )

## Transformation matrices for cingular pole template registration

FileType( 'Talairach To Subject Transformation', 'Transformation matrix' )
FileType( 'Subject To Template Transformation', 'Transformation matrix' )

## Gyri files

FileType( 'Left Gyri Graph', 'Data graph' )
FileType( 'Right Gyri Graph', 'Data graph' )

FileType( 'Left Gyri Volume', '4D Volume' )
FileType( 'Right Gyri Volume', '4D Volume' )

## Projection using convolution kernels

FileType( 'Projection convolution kernels', '4D Volume' )
FileType( 'Functional texture', 'Texture' )

# Cortical thickness

FileType( 'Cortical thickness', 'Texture' )

FileType( 'Surface Label Translation', 'Label Translation' )

# Sulci parameterizations

FileType( 'Sulcus mesh' , 'Mesh' )
FileType( 'Sulcus x coordinate texture', 'Texture' )
FileType( 'Sulcus y coordinate texture', 'Texture' ) 
FileType( 'Sulcus coordinate grid mesh', 'Mesh' )
FileType( 'Sulcus depth profile', 'Text file' )
