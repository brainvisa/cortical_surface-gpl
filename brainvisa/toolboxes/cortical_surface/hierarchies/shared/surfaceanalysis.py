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
include( 'base' )
include( 'sulci' )

insert( 'nomenclature','surfaceanalysis',
  SetContent(
    'constraint_correspondance', SetType( 'Constraint coordinates values' ),SetWeakAttr( 'filename_variable', 'constraint_correspondance'),
    'constraint_correspondance_2012', SetType( 'Constraint coordinates values' ),SetWeakAttr( 'filename_variable', 'constraint_correspondance_2012'), SetPriorityOffset( +1 ),
    'surfaceReferential', SetType( 'Label Translation' ),SetWeakAttr( 'filename_variable', 'surfaceReferential'),SetWeakAttr( 'category', 'translation'),
    'surfaceRefModel_par', SetType( 'Latitude Constraint Gyri Model' ),
    'surfaceRefModel_mer', SetType( 'Longitude Constraint Gyri Model' ),
    'sulcal_coordinates_model_left', SetType( 'HipHop Model' ),
      SetWeakAttr( 'side', 'left' ),
    'sulcal_coordinates_model_right', SetType( 'HipHop Model' ),
      SetWeakAttr( 'side', 'right' ),
  )
)


insert( 'hemitemplate',
  "*Cingular", SetType( "Cingular Pole Template" ),
#  "*PoleLeft", SetType( "Left Cingular Pole Template" ), SetWeakAttr( 'side', 'left' ),
#  "*PoleRight", SetType( "Right Cingular Pole Template" ), SetWeakAttr( 'side', 'right' ),
)

insert( 'models/models_{sulci_database}/discriminative_models/{graph_version}/{model}',
    'gyri', SetType( 'Gyri Model' ),
)

