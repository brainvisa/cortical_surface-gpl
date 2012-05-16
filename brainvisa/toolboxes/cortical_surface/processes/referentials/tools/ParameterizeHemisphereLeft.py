
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

name = 'Parameterize Left hemisphere'

userLevel = 2


signature = Signature(
    'Side', Choice("Left"),
    'left_white_mesh',ReadDiskItem( 'Left Hemisphere White Mesh' , shfjGlobals.aimsMeshFormats),
    'left_white_sulci_par',ReadDiskItem( 'Left hemisphere latitude cleaned constraints texture', 'Texture',requiredAttributes={ 'side': 'left' }  ),
    'left_white_sulci_mer',ReadDiskItem( 'Left hemisphere longitude cleaned constraints texture', 'Texture',requiredAttributes={ 'side': 'left' }  ),
    'left_cingular_pole',ReadDiskItem( 'Left hippocampus pole texture'  , 'Texture',requiredAttributes={ 'side': 'left' } ),
    'left_poles_texture',ReadDiskItem( 'Left poles texture'  , 'Texture',requiredAttributes={ 'side': 'left' } ),
    'origin_mer', Integer(),
    'stop_criterium', Float(),
    'time_step', Float(),
    'process', Choice("Both","Longitude","Latitude","None"),
    'constraint_attach', Float(),
    'type_beta', Choice("Map","Value"),
    'left_latitude',WriteDiskItem( 'Left hemisphere latitude texture','Texture',requiredAttributes={ 'side': 'left' } ),
    'left_longitude',WriteDiskItem( 'Left hemisphere longitude texture','Texture',requiredAttributes={ 'side': 'left' })
)

def initialization( self ):
#    self.linkParameters( 'left_cingular_pole','left_white_mesh')
    self.linkParameters( 'left_latitude','left_white_mesh')
    self.linkParameters( 'left_longitude','left_white_mesh')
    self.setOptional( 'origin_mer', 'stop_criterium', 'time_step', 'process', 'constraint_attach', 'left_latitude', 'left_longitude', 'left_white_sulci_par', 'left_white_sulci_mer','type_beta' )
    self.stop_criterium = 1e-6
    self.time_step = 0.1
    self.constraint_attach = 0.2
    self.origin_mer = 0

def execution( self, context ):
    context.write(self.Side)
    if self.process == 'Longitude' :
        process_choice = 2
    else :
        if self.process == 'Latitude' :
            process_choice = 1
        else :
            if self.process == 'None' :
                process_choice = 3
            else :
                process_choice = 0

    if self.type_beta == 'Value' :
        t_beta = 0
    else :
        t_beta = 1

    context.write('Left hemisphere')
    context.system('AimsCorticalReferential', '-i', self.left_white_mesh.fullPath() ,  '-p' , self.left_white_sulci_par.fullPath(), '-m' , self.left_white_sulci_mer.fullPath() , '-r', self.left_poles_texture.fullPath() ,'-l', self.left_cingular_pole.fullPath() , '-c' , self.stop_criterium, '-d' , self.time_step, '-b' , self.origin_mer , '-f' , process_choice , '-t' , self.constraint_attach , '-a', t_beta , '-x' , self.left_latitude.fullPath() , '-y' , self.left_longitude.fullPath() )
    context.write('Done')
