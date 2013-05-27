# -*- coding: utf-8 -*-
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

name = 'Sulcal Lines Extraction'

userLevel = 2

def validation():
    anatomist.validation()
    
signature = Signature(
                      
    'side', Choice('left', 'right'),
    'white_mesh',ReadDiskItem( 'Hemisphere White Mesh' , shfjGlobals.aimsMeshFormats),
    'graph', ReadDiskItem( 'Cortical folds graph', 'Graph'),
    'grey_white_input', ReadDiskItem( 'Right Grey White Mask', shfjGlobals.anatomistVolumeFormats),
    'sulcus_identification',Choice('name','label'),
    'labels_translation_map',ReadDiskItem( 'Label Translation' ,'Label Translation'),
    'graph_label_basins',WriteDiskItem( 'Right Graph Label Translation', 'Text File'),
    'file_correspondance_constraint',ReadDiskItem( 'Constraint coordinates values', 'Text File'),
    'mri', ReadDiskItem( "Raw T1 MRI", shfjGlobals.vipVolumeFormats ),
    'bucket_label_type', Choice('All', 'aims_junction', 'aims_bottom', 'aims_ss', 'aims_other'),  
    'white_sulcalines',WriteDiskItem( 'Right hemisphere Sulcal Lines texture', 'Texture'),
    'basin_min_size', Float(),
    'constraint_weight', Integer(),
    'constraint_value', Choice('Basins Label','Lat/Lon'),
)

def initialization( self ):
      
    self.linkParameters( 'mri', 'white_mesh' )
    self.linkParameters( 'graph','white_mesh')
    self.linkParameters( 'grey_white_input','white_mesh')
    self.linkParameters( 'white_sulcalines', 'white_mesh' )
    self.linkParameters( 'graph_label_basins','white_mesh')
    self.linkParameters( 'graph_label_basins','graph')
    self.findValue( 'labels_translation_map', {'filename_variable' : 'sulci_model_2008'} )
    self.findValue( 'file_correspondance_constraint', {'filename_variable' : 'constraint_correspondance_2012'} )
    self.bucket_label_type = 'All'
    self.sulcus_identification = 'label'
    self.basin_min_size = 50.0
    self.constraint_weight = 15
    self.constraint_value = 'Basins Label'
    
def execution( self, context ):
        
    curvatureIm=context.temporary(  'GIS image' )
    smoothIm=context.temporary( 'GIS image')

    context.write('computing right curvature texture')
    curv = [ 'AimsMeshCurvature',
                '-i', self.white_mesh.fullPath(),
                '-o', curvatureIm.fullPath(),
                '-m', 'barycenter' ]

    apply( context.system, curv )
    context.write('Done')

    context.write('smoothing right curvature texture')
    smooth = [ 'AimsTextureSmoothing',
                '-i', curvatureIm.fullPath(),
                '-o', smoothIm.fullPath(),
                '-m', self.white_mesh.fullPath(),
                '-s', 2,
                '-t', 0.01
                ]
    apply( context.system, smooth )
    context.write('Done')

    depthIm=context.temporary( 'GIS image')

    context.write('computing right depth texture')
    depth = [ 'AimsMeshGeodesicDepth',
                '-v', self.grey_white_input.fullPath(),
                '-i', self.white_mesh.fullPath(),
                '-c', 10,
                '-e', 8,
                '-o', depthIm.fullPath()
                ]
    apply( context.system, depth )
    context.write('Done')

    volumeGraphLabelBasins=context.temporary('NIFTI-1 image')

    context.write('computing Graph Label')

    if (self.bucket_label_type=='aims_junction'):
      graphBucketLabel = [ 'siGraph2Label','-g', self.graph.fullPath(),'-a', self.sulcus_identification,'-tv', self.mri.fullPath(),'-tr', self.labels_translation_map.fullPath(),'-o', volumeGraphLabelBasins.fullPath(),
      '-b', 'aims_junction', '-ot', self.graph_label_basins.fullPath()]
    elif (self.bucket_label_type=='aims_bottom'):
      graphBucketLabel = [ 'siGraph2Label','-g', self.graph.fullPath(),'-a', self.sulcus_identification,'-tv', self.mri.fullPath(),'-tr', self.labels_translation_map.fullPath(),'-o', volumeGraphLabelBasins.fullPath(),
      '-b', 'aims_bottom', '-ot', self.graph_label_basins.fullPath()]
    elif (self.bucket_label_type=='aims_ss'):
      graphBucketLabel = [ 'siGraph2Label','-g', self.graph.fullPath(),'-a', self.sulcus_identification,'-tv', self.mri.fullPath(),'-tr', self.labels_translation_map.fullPath(),'-o', volumeGraphLabelBasins.fullPath(),
      '-b', 'aims_ss', '-ot', self.graph_label_basins.fullPath()]
    elif (self.bucket_label_type=='aims_other'):
      graphBucketLabel = [ 'siGraph2Label','-g', self.graph.fullPath(),'-a', self.sulcus_identification,'-tv', self.mri.fullPath(),'-tr', self.labels_translation_map.fullPath(),'-o', volumeGraphLabelBasins.fullPath(),
      '-b', 'aims_other', '-ot', self.graph_label_basins.fullPath()]
    else : graphBucketLabel = [ 'siGraph2Label','-g', self.graph.fullPath(),'-a', self.sulcus_identification,'-tv', self.mri.fullPath(),'-tr', self.labels_translation_map.fullPath(),'-o', volumeGraphLabelBasins.fullPath(),
      '-b', 'aims_junction','-b', 'aims_bottom','-b', 'aims_ss', '-b', 'aims_other', '-ot', self.graph_label_basins.fullPath()]
#
    apply( context.system, graphBucketLabel )

#      context.runProcess
    context.write('Done')

    if self.constraint_value == 'Basins Label' :
      constraintValue = 1
    else :
      constraintValue = 2

    context.write('Sulcal Lines extraction')

    sulcalines = [ 'AimsSulcalLines',
                '-d', depthIm.fullPath(),
                '-c', smoothIm.fullPath(),
                '-i', self.white_mesh.fullPath(),
                '-b', volumeGraphLabelBasins.fullPath(),
                '-lb',self.graph_label_basins.fullPath(),
                '-ls',self.file_correspondance_constraint.fullPath(),
                '-m', 3,
                '-t', 2,
                '-st', self.constraint_weight,
                '-o', self.white_sulcalines.fullPath(),
                '-si', 'right', #self.side,
                '-sb', self.basin_min_size,
                '-cv', constraintValue
                ]
    apply( context.system, sulcalines )

    context.write('Done')

    