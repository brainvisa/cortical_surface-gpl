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

from __future__ import absolute_import
from brainvisa.processes import *

name = '2D Parcellation to 3D parcellation'

userLevel = 0

signature = Signature(
    'Side', Choice("Both","Left","Right"),
    'left_gyri',ReadDiskItem( 'hemisphere parcellation texture','Aims texture formats',requiredAttributes={ 'side': 'left' } ),
    'right_gyri',ReadDiskItem( 'hemisphere parcellation texture','Aims texture formats',requiredAttributes={ 'side': 'right' } ),
    'left_white_mesh',ReadDiskItem( 'Left Hemisphere White Mesh' , 'aims Mesh Formats'),
    'right_white_mesh',ReadDiskItem( 'Right Hemisphere White Mesh' , 'aims Mesh Formats'),
    'texture_time', Integer(),
    'object_label', Integer(),
#    'translation',ReadDiskItem('Label Translation','Label Translation' ),
#    'translation_type', Choice("int_to_string","string_to_int"),
    'left_input_volume', ReadDiskItem( 'Left Grey White Mask', 'anatomist Volume Formats' ),
    'right_input_volume', ReadDiskItem( 'Right Grey White Mask', 'anatomist Volume Formats' ),
    'left_output_volume', WriteDiskItem( 'Left Gyri Volume', 'Aims writable volume formats' ), #parcellation volume', 'Aims writable volume formats' ),
    'right_output_volume', WriteDiskItem( 'Right Gyri Volume', 'Aims writable volume formats' ), #parcellation volume', 'Aims writable volume formats' ),
#    'left_output_graph', WriteDiskItem( 'Left Gyri Graph', 'Graph' ),#Parcels Graph', 'Graph' ),
#    'right_output_graph', WriteDiskItem( 'Right Gyri Graph', 'Graph' ),#Parcels Graph', 'Graph' ),

)

def initialization( self ):
    self.linkParameters( 'right_gyri', 'left_gyri' )
    self.linkParameters( 'left_gyri', 'right_gyri' )
    self.linkParameters( 'left_white_mesh', 'left_gyri' )
    self.linkParameters( 'right_white_mesh', 'right_gyri' )
    self.linkParameters( 'right_white_mesh', 'left_white_mesh' )
    self.linkParameters( 'left_white_mesh', 'right_white_mesh' )
    self.linkParameters( 'left_input_volume', 'left_white_mesh' )
    self.linkParameters( 'right_input_volume', 'right_white_mesh' )
    self.linkParameters( 'left_output_volume', 'left_white_mesh' )
    self.linkParameters( 'right_output_volume', 'right_white_mesh' )
#    self.linkParameters( 'left_output_graph', 'left_gyri' )
#    self.linkParameters( 'right_output_graph', 'right_gyri' )
    self.texture_time = 0
#    self.setOptional('left_white_mesh','right_white_mesh','left_gyri','right_gyri')#,'translation','object_label','left_output_graph','right_output_graph' )
    self.object_label = 100

def execution( self, context ): 
#    graph_version = '3.0'
    if self.Side in ('Left','Both'):    
        context.write('3D parcellation.')
        context.system('AimsMeshParcellation2VolumeParcellation',
                    '-m', self.left_white_mesh.fullPath(),
                    '-t', self.left_gyri.fullPath(),
                    '-o', self.left_output_volume.fullPath(),
                    '-v', self.left_input_volume.fullPath(),
                    '-T', self.texture_time,
                    '-l', self.object_label )
            
            
#        if self.left_output_graph is not None:
           #Build Roi graph from texture
#           call_list = ['AimsTex2Graph',
#                       '-m', self.left_white_mesh.fullPath(),
#                       '-t', self.left_gyri.fullPath(),
#                       '-o', self.left_output_graph.fullPath(),
#                       '-T', self.texture_time]
        
#           if self.translation is not None :
#               options = ['-c', self.translation.fullPath()]
#               if self.translation_type == "string_to_int":
#                   options += ['--reverse']
#               call_list += options
        
#           context.write('Convert mesh+texture to ROI graph.')
#           apply( context.system, call_list)
    
#           tmp = context.temporary( 'Graph' )
#           context.system('AimsGraphConvert',
#                           '-i', self.left_output_volume.fullPath(),
#                           '-o', tmp,
#                           '--bucket')
#           context.system('AimsGraphMerge',
#                           '-i', tmp,
#                           '-j', self.left_output_graph.fullPath(),
#                           '-k', 'roi_label',
#                           '-o', self.left_output_graph.fullPath())
    
#           context.system('AimsGraphComplete',
#                           '-i', self.left_output_graph.fullPath(),
#                           '-o', self.left_output_graph.fullPath(),
#                           '--dversion', graph_version,
#                           '--mversion', graph_version)
    
    if self.Side in ('Right','Both'):
        if self.right_input_volume is not None :    
            context.write('3D parcellation.')
            context.system('AimsMeshParcellation2VolumeParcellation',
                        '-m', self.right_white_mesh.fullPath(),
                        '-t', self.right_gyri.fullPath(),
                        '-o', self.right_output_volume.fullPath(),
                        '-v', self.right_input_volume.fullPath(),
                        '-T', self.texture_time,
                        '-l', self.object_label )
    

#        if self.right_output_graph is not None:
           #Build Roi graph from texture
#           call_list = ['AimsTex2Graph',
#                       '-m', self.right_white_mesh.fullPath(),
#                       '-t', self.right_gyri.fullPath(),
#                       '-o', self.right_output_graph.fullPath(),
#                       '-T', self.texture_time]
        
#           if self.translation is not None :
#               options = ['-c', self.translation.fullPath()]
#               if self.translation_type == "string_to_int":
#                   options += ['--reverse']
#               call_list += options
        
#           context.write('Convert mesh+texture to ROI graph.')
#           apply( context.system, call_list)
        
#           tmp = context.temporary( 'Graph' )
#           context.system('AimsGraphConvert',
#                           '-i', self.right_output_volume.fullPath(),
#                           '-o', tmp,
#                           '--bucket')
#           context.system('AimsGraphMerge',
#                           '-i', tmp,
#                           '-j', self.right_output_graph.fullPath(),
#                           '-k', 'roi_label',
#                           '-o', self.right_output_graph.fullPath())
#        
        
#           context.system('AimsGraphComplete',
#                           '-i', self.right_output_graph.fullPath(),
#                           '-o', self.right_output_graph.fullPath(),
#                           '--dversion', graph_version,
#                           '--mversion', graph_version)
    
