
include( 'base' )

insert( 'nomenclature','surfaceanalysis',
  SetContent(
    'constraint_correspondance', SetType( 'Constraint coordinates values' ),
    'surfaceReferential', SetType( 'Surface Label Translation' ),
    'surfaceRefModel_par', SetType( 'Latitude Constraint Gyri Model' ),
    'surfaceRefModel_mer', SetType( 'Longitude Constraint Gyri Model' ),
  )
)


insert( 'hemitemplate',
  "*PoleLeft", SetType( "Left Cingular Pole Template" ), 
  "*PoleRight", SetType( "Right Cingular Pole Template" ), 
)
