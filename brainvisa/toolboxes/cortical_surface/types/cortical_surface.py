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
include( 'builtin' )
include( 'registration' )
include( 'anatomy' )

## sulcal pits
FileType( 'pits texture', 'Texture' )
FileType( 'Left pits texture', 'pits texture' )
FileType( 'Right pits texture', 'pits texture' )
FileType( 'noisy pits texture', 'Texture' )
FileType( 'Left noisy pits texture', 'noisy pits texture' )
FileType( 'Right noisy pits texture', 'noisy pits texture' )
FileType( 'ridges texture', 'Texture' )
FileType( 'Left ridges texture', 'ridges texture' )
FileType( 'Right ridges texture', 'ridges texture' )
FileType( 'basins texture', 'Texture' )
FileType( 'Left basins texture', 'basins texture' )
FileType( 'Right basins texture', 'basins texture' )

## Depth Potential Function
FileType( 'DPF texture', 'Texture' )
FileType( 'Left DPF texture', 'DPF texture' )
FileType( 'Right DPF texture', 'DPF texture' )

## surface coordinates
FileType( 'Poles texture', 'Texture' )
FileType( 'Left poles texture', 'Poles texture' )
FileType( 'Right poles texture', 'Poles texture' )
FileType( 'Insula pole texture', 'Poles texture')
FileType( 'Left insula pole texture', 'Insula pole texture' )
FileType( 'Right insula pole texture', 'Insula pole texture' )
FileType( 'Hippocampus pole texture', 'Poles texture' )
FileType( 'Left hippocampus pole texture', 'Hippocampus pole texture' )
FileType( 'Right hippocampus pole texture', 'Hippocampus pole texture' )

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
FileType( 'Latitude coordinate HIP texture', 'Coordinate texture' )
FileType( 'Longitude coordinate HIP texture', 'Coordinate texture' )
FileType( 'Left hemisphere longitude HIP texture', 'Longitude coordinate HIP texture' )
FileType( 'Left hemisphere latitude HIP texture', 'Latitude coordinate HIP texture' )
FileType( 'Right hemisphere longitude HIP texture', 'Longitude coordinate HIP texture' )
FileType( 'Right hemisphere latitude HIP texture', 'Latitude coordinate HIP texture' )
# useless?
FileType( 'Conformal longitude texture', 'Coordinate texture' )
FileType( 'Conformal latitude texture', 'Coordinate texture' )

FileType( 'Coordinate grid', 'Mesh' )
FileType( 'Left hemisphere coordinate grid', 'Coordinate grid')
FileType( 'Right hemisphere coordinate grid', 'Coordinate grid')

FileType( 'Remeshed mesh', 'Mesh' )
FileType( 'Left remeshed mesh', 'Remeshed mesh')
FileType( 'Right remeshed mesh', 'Remeshed mesh')

FileType( 'Flat mesh', 'Mesh' )
FileType( 'Rectangular flat mesh', 'Flat mesh')
FileType( 'Left rectangular flat mesh', 'Rectangular flat mesh')
FileType( 'Right rectangular flat mesh', 'Rectangular flat mesh')
FileType( 'Rectangular flat cstr mesh', 'Flat mesh')
FileType( 'Left rectangular flat cstr mesh', 'Rectangular flat cstr mesh')
FileType( 'Right rectangular flat cstr mesh', 'Rectangular flat cstr mesh')
FileType( 'Disk flat mesh', 'Flat mesh')
FileType( 'Left disk flat mesh', 'Disk flat mesh')
FileType( 'Right disk flat mesh', 'Disk flat mesh')

## surface_tools for parameterization
FileType( 'Boundary texture', 'Texture' )
FileType( 'Rectangular boundary texture', 'Boundary texture' )
FileType( 'Left rectangular boundary texture', 'Rectangular boundary texture' )
FileType( 'Right rectangular boundary texture', 'Rectangular boundary texture' )
FileType( 'Indices texture', 'Texture' )
FileType( 'Rectangular flat indices texture', 'Indices texture' )
FileType( 'Left rectangular flat indices texture', 'Rectangular flat indices texture' )
FileType( 'Right rectangular flat indices texture', 'Rectangular flat indices texture' )
FileType( 'Rectangular flat texture','Texture' )
FileType( 'hemisphere Sulcal Lines Rectangular Flat texture', 'Texture' )
FileType( 'Left hemisphere Sulcal Lines Rectangular Flat texture', 'hemisphere Sulcal Lines Rectangular Flat texture' )
FileType( 'Right hemisphere Sulcal Lines Rectangular Flat texture', 'hemisphere Sulcal Lines Rectangular Flat texture' )
FileType( 'White mesh parts', 'Mesh' )
FileType( 'Left white mesh parts', 'White mesh parts')
FileType( 'Right white mesh parts', 'White mesh parts')

## template for poles
FileType( 'Cingular Pole Template', '3D Volume' )
FileType( 'Left Cingular Pole Template', 'Cingular Pole Template' )
FileType( 'Right Cingular Pole Template', 'Cingular Pole Template' )

FileType( 'Cingular Pole Template Subject', '3D Volume' )
FileType( 'Left Cingular Pole Template Subject', 'Cingular Pole Template Subject' )
FileType( 'Right Cingular Pole Template Subject', 'Cingular Pole Template Subject' )

## model files
FileType( 'Latitude Constraint Gyri Model', 'Gyri Model' )
FileType( 'Longitude Constraint Gyri Model', 'Gyri Model' )

## HIipHop model files
FileType( 'HipHop Model', 'Text File' )

## Transformation matrices for cingular pole template registration
FileType( 'Talairach To Subject Transformation', 'Transformation matrix' )
FileType( 'Subject To Template Transformation', 'Transformation matrix' )

## Tranformation matrices for registration between function and anatomy
FileType( 'Anatomy To Mean Functional Volume Transformation', 'Transformation matrix' )
FileType( 'Mean Functional Volume To Anatomy Transformation', 'Transformation matrix' )

## Gyri files !! different from Parcels files
FileType( 'Hemisphere gyri parcellation texture', 'Gyri White Texture' )
FileType( 'Left hemisphere gyri parcellation texture', 'Hemisphere gyri parcellation texture' )
FileType( 'Right hemisphere gyri parcellation texture', 'Hemisphere gyri parcellation texture' )
FileType( 'Left hemisphere regularized parcellation texture', 'Hemisphere gyri parcellation texture')
FileType( 'Right hemisphere regularized parcellation texture', 'Hemisphere gyri parcellation texture')

FileType( 'Left Gyri Graph', 'Data graph' )
FileType( 'Right Gyri Graph', 'Data graph' )

FileType( 'Left Gyri Volume', '4D Volume' )
FileType( 'Right Gyri Volume', '4D Volume' )

## Parcels files !! different from Gyri files
FileType( 'Hemisphere parcellation texture', 'Gyri White Texture' )
FileType( 'Left hemisphere parcellation texture', 'Hemisphere parcellation texture' )
FileType( 'Right hemisphere parcellation texture', 'Hemisphere parcellation texture' )

FileType( 'Left hemisphere model parcellation texture', 'Left hemisphere parcellation texture' )
FileType( 'Right hemisphere model parcellation texture', 'Right hemisphere parcellation texture' )

FileType( 'Left hemisphere marsAtlas parcellation texture', 'Left hemisphere parcellation texture')
FileType( 'Right hemisphere marsAtlas parcellation texture', 'Right hemisphere parcellation texture')

# may be used to replace hemisphere regularized parcellation texture
#FileType( 'Left hemisphere regularized model parcellation texture', 'Hemisphere gyri parcellation texture')
#FileType( 'Right hemisphere regularized model parcellation texture', 'Hemisphere gyri parcellation texture')

FileType( 'Parcels Graph', 'Data graph' )
FileType( 'Left parcels Graph', 'Parcels Graph' )
FileType( 'Right parcels Graph', 'Parcels Graph' )

FileType( 'Parcellation volume', '4D Volume' )
FileType( 'Left parcellation volume', 'Parcellation volume' )
FileType( 'Right parcellation volume', 'Parcellation volume' )
FileType( 'Left model parcellation volume', 'Left parcellation volume' )
FileType( 'Right model parcellation volume', 'Right parcellation volume' )
FileType( 'Left coarse parcellation volume', 'Left parcellation volume' )
FileType( 'Right coarse parcellation volume', 'Right parcellation volume' )

## Projection using convolution kernels
FileType( 'Labelled Functional Blobs Texture', 'Label Texture' )

## Types related to surface-based functional analysis
FileType( 'Projection convolution kernels', '4D Volume' )
FileType( 'Functional Texture', 'Texture' )
FileType( 'Functional Time Texture', 'Texture' )
FileType( 'Surface-Based SPMt Map', 'Texture' )
FileType( 'Surface-Based Beta Map', 'Texture' )

# Cortical thickness
FileType( 'Cortical thickness', 'Texture' )

# Sulci parameterizations
FileType( 'Sulcus mesh' , 'Mesh' )
FileType( 'Sulcus x coordinate texture', 'Texture' )
FileType( 'Sulcus y coordinate texture', 'Texture' ) 
FileType( 'Sulcus coordinate grid mesh', 'Mesh' )
FileType( 'Sulcus depth profile', 'Text file' )

# Sulcal lines extraction
FileType( 'hemisphere Sulcal Lines texture', 'Texture' )
FileType( 'Left hemisphere Sulcal Lines texture', 'hemisphere Sulcal Lines texture' )
FileType( 'Right hemisphere Sulcal Lines texture', 'hemisphere Sulcal Lines texture' )

FileType( 'Graph Label Translation', 'Text File' )
FileType( 'Left Graph Label Translation', 'Graph Label Translation' )
FileType( 'Right Graph Label Translation', 'Graph Label Translation' )



