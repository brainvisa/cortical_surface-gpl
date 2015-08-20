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
from brainvisa.processes import *
import shfjGlobals   

name = 'Constraint Cleaner Left hemisphere'

userLevel = 2


signature = Signature(
    'left_white_mesh',ReadDiskItem( 'Left Hemisphere White Mesh' , shfjGlobals.aimsMeshFormats),
    'left_cingular_pole',ReadDiskItem( 'Left cingular pole texture'  , 'Texture',requiredAttributes={ 'side': 'left' } ),
    'left_white_sulci_mer',ReadDiskItem( 'Left hemisphere longitude constraints texture', 'Texture',requiredAttributes={ 'side': 'left' }  ),
    'left_white_sulci_par',ReadDiskItem( 'Left hemisphere latitude constraints texture', 'Texture',requiredAttributes={ 'side': 'left' }  ),
    'file_correspondance_constraint',ReadDiskItem( 'Constraint coordinates values', 'Text File' ),
    'left_sulci_label_to_sulci_name',ReadDiskItem( 'Sulci To White Texture Translation', 'Text File' ),
    'constraint_dist_param', Float(),
    'curvature_param', Float(),
    'elasticity_param', Float(),
    'left_poles_texture',WriteDiskItem( 'Left poles texture'  , 'Texture',requiredAttributes={ 'side': 'left' } ),
    'left_white_sulci_mer_cleaned',WriteDiskItem( 'Left hemisphere longitude cleaned constraints texture', 'Texture',requiredAttributes={ 'side': 'left' }  ),
    'left_white_sulci_par_cleaned',WriteDiskItem( 'Left hemisphere latitude cleaned constraints texture', 'Texture',requiredAttributes={ 'side': 'left' }  ),
    'process', Choice("Both","Longitude","Latitude" )
)

def initialization( self ):
    #self.linkParameters( 'left_cingular_pole','left_white_mesh')
    #self.linkParameters( 'left_white_sulci_mer','left_white_mesh')
    #self.linkParameters( 'left_white_sulci_par','left_white_mesh')
    #self.linkParameters( 'left_cingular_pole','left_white_sulci_mer')
    #self.linkParameters( 'left_cingular_pole','left_white_sulci_par')
    #self.linkParameters( 'left_white_sulci_par','left_white_sulci_mer')
    #self.linkParameters( 'left_white_sulci_mer','left_white_sulci_par')
    self.linkParameters( 'left_white_sulci_mer_cleaned','left_white_mesh')
    self.linkParameters( 'left_white_sulci_par_cleaned','left_white_mesh')
    self.linkParameters( 'left_poles_texture','left_white_mesh')
    #self.linkParameters( 'left_sulci_label_to_sulci_name', 'left_white_mesh' )
    self.setOptional( 'process', 'left_white_sulci_mer_cleaned', 'left_white_sulci_par_cleaned', 'constraint_dist_param', 'curvature_param', 'elasticity_param'  )
    self.findValue( 'file_correspondance_constraint', {} )

    self.constraint_dist_param = 20
    self.curvature_param = 500
    self.elasticity_param = 800
    
def execution( self, context ):
    if self.process == 'Longitude' :
        process_choice = 0
    else :
        if self.process == 'Latitude' :
            process_choice = 1
        else :
            process_choice = 2

    side = '_left'
    
    
    context.write('Left hemisphere')
    context.system('AimsConstraintCleaner', '-m', self.left_white_mesh.fullPath() ,  '-t', self.left_cingular_pole.fullPath() ,  '-p' , self.left_poles_texture.fullPath(), '-x' , self.left_white_sulci_mer.fullPath() , '-y', self.left_white_sulci_par.fullPath() , '-a' , self.left_white_sulci_mer_cleaned.fullPath(), '-b', self.left_white_sulci_par_cleaned.fullPath(), '-f' , self.file_correspondance_constraint.fullPath(), '-g', self.left_sulci_label_to_sulci_name.fullPath(), '-c',  process_choice, '-s', side, '-i', self.constraint_dist_param, '-j', self.curvature_param , '-k' , self.elasticity_param )
    context.write('Done')
