
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
from brainvisa import anatomist

name = 'Right Cingular Pole Projection'

userLevel = 2

def validation():
    anatomist.validation()

signature = Signature(
    'Side', Choice("Right"),
    'right_white_mesh',ReadDiskItem( 'Right Hemisphere White Mesh' , shfjGlobals.aimsMeshFormats),
    'right_pole_template',ReadDiskItem( 'Right Cingular Pole Template Subject' , 'Aims readable volume formats' ),
    'right_pole',WriteDiskItem( 'Right hippocampus pole texture','Texture',requiredAttributes={ 'side': 'right' } )
)

def initialization( self ):
    self.linkParameters( 'right_pole', 'right_white_mesh' )
    #self.findValue( 'right_pole_template', {} )
    #self.setOptional('right_pole_template')
     
def execution( self, context ):
    a = anatomist.Anatomist()
    mesh = a.loadObject( self.right_white_mesh.fullPath() )
    vol = a.loadObject( self.right_pole_template.fullPath() )
    fusion = a.fusionObjects( [mesh, vol], method='Fusion3DMethod' )
    fusion.exportTexture(filename=self.right_pole.fullPath())

