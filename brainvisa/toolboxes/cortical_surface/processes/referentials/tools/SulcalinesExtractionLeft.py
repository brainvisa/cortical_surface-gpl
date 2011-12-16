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
from neuroProcesses import *
import shfjGlobals   
from brainvisa import anatomist
import sigraph

name = 'Sulcal Lines Extraction Left'

userLevel = 2

def validation():
    anatomist.validation()
    
signature = Signature(
                      
    'Side', Choice('left'),
    'Lwhite_mesh',ReadDiskItem( 'Hemisphere White Mesh' , shfjGlobals.aimsMeshFormats),    
    'Lgraph', ReadDiskItem( 'Cortical folds graph', 'Graph',requiredAttributes={ 'side': 'left' } ),
    'Lgrey_white_input', ReadDiskItem( 'Left Grey White Mask', shfjGlobals.anatomistVolumeFormats,requiredAttributes={ 'side': 'left' } ),
    'sulcus_identification',Choice('name','label'),
    'labels_translation_map',ReadDiskItem( 'Label Translation' ,'Label Translation'),
    'Lgraph_label_basins',WriteDiskItem( 'Left Graph Label Translation', 'Text File',requiredAttributes={ 'side': 'left' }),
    'file_correspondance_constraint',ReadDiskItem( 'Constraint coordinates values', 'Text File'),
    'mri', ReadDiskItem( "Raw T1 MRI", shfjGlobals.vipVolumeFormats ),
    'bucket_label_type', Choice('All', 'aims_junction', 'aims_bottom', 'aims_ss', 'aims_other'),  
    'Lwhite_sulcalines',WriteDiskItem( 'Left hemisphere Sulcal Lines texture', 'Texture' ,requiredAttributes={ 'side': 'left' } ),
    'basin_min_size', Float(),
    'constraint_weight', Integer(),
    'constraint_value', Choice('Basins Label','Lat/Lon'),
)

def initialization( self ):
      
    self.linkParameters( 'mri', 'Lwhite_mesh' )
    self.linkParameters( 'Lgraph','Lwhite_mesh')
    self.linkParameters( 'Lgrey_white_input','Lwhite_mesh')
    self.linkParameters( 'Lwhite_sulcalines', 'Lwhite_mesh' )
    self.linkParameters( 'Lgraph_label_basins','Lwhite_mesh')
    self.linkParameters( 'Lgraph_label_basins','Lgraph')
    self.findValue( 'labels_translation_map', {'filename_variable' : 'sulci_model_2008'} )
    self.findValue( 'file_correspondance_constraint', {'filename_variable' : 'constraint_correspondance_2012'} )
    self.bucket_label_type = 'All'
    self.side = 'left'
    self.sulcus_identification = 'label'
    self.basin_min_size = 50.0
    self.constraint_weight = 15
    self.constraint_value = 'Basins Label'
    
def execution( self, context ):
        
    if self.side in ('left'):  
      curvatureIm=context.temporary(  'GIS image' )
      smoothIm=context.temporary( 'GIS image')
      
      context.write('computing left curvature texture')
      curv = [ 'AimsMeshCurvature',
                 '-i', self.Lwhite_mesh.fullPath(),
                 '-o', curvatureIm.fullPath(),
                 '-m', 'barycenter' ]

      apply( context.system, curv )     
      context.write('Done')
      
      context.write('smoothing left curvature texture')
      smooth = [ 'AimsTextureSmoothing',
                 '-i', curvatureIm.fullPath(),
                 '-o', smoothIm.fullPath(),
                 '-m', self.Lwhite_mesh.fullPath(),
                 '-s', 2,
                 '-t', 0.01
                 ]
      apply( context.system, smooth )
      context.write('Done')
      
      depthIm=context.temporary( 'GIS image')
      
      context.write('computing left depth texture')
      depth = [ 'AimsMeshGeodesicDepth',
                 '-v', self.Lgrey_white_input.fullPath(),
                 '-i', self.Lwhite_mesh.fullPath(),
                 '-c', 10,
                 '-e', 8,
                 '-o', depthIm.fullPath()
                 ]
      apply( context.system, depth )
      context.write('Done')
  
      volumeGraphLabelBasins=context.temporary('NIFTI-1 image')

      context.write('computing Graph Label')
      
      if (self.bucket_label_type=='aims_junction'):
        graphBucketLabel = [ 'siGraph2Label','-g', self.Lgraph.fullPath(),'-a', self.sulcus_identification,'-tv', self.mri.fullPath(),'-tr', self.labels_translation_map.fullPath(),'-o', volumeGraphLabelBasins.fullPath(),
        '-b', 'aims_junction', '-ot', self.Lgraph_label_basins.fullPath()]
      elif (self.bucket_label_type=='aims_bottom'):
        graphBucketLabel = [ 'siGraph2Label','-g', self.Lgraph.fullPath(),'-a', self.sulcus_identification,'-tv', self.mri.fullPath(),'-tr', self.labels_translation_map.fullPath(),'-o', volumeGraphLabelBasins.fullPath(),
        '-b', 'aims_bottom', '-ot', self.Lgraph_label_basins.fullPath()]
      elif (self.bucket_label_type=='aims_ss'):
        graphBucketLabel = [ 'siGraph2Label','-g', self.Lgraph.fullPath(),'-a', self.sulcus_identification,'-tv', self.mri.fullPath(),'-tr', self.labels_translation_map.fullPath(),'-o', volumeGraphLabelBasins.fullPath(),
        '-b', 'aims_ss', '-ot', self.Lgraph_label_basins.fullPath()]
      elif (self.bucket_label_type=='aims_other'):
        graphBucketLabel = [ 'siGraph2Label','-g', self.Lgraph.fullPath(),'-a', self.sulcus_identification,'-tv', self.mri.fullPath(),'-tr', self.labels_translation_map.fullPath(),'-o', volumeGraphLabelBasins.fullPath(),
        '-b', 'aims_other', '-ot', self.Lgraph_label_basins.fullPath()]
      else : graphBucketLabel = [ 'siGraph2Label','-g', self.Lgraph.fullPath(),'-a', self.sulcus_identification,'-tv', self.mri.fullPath(),'-tr', self.labels_translation_map.fullPath(),'-o', volumeGraphLabelBasins.fullPath(),
        '-b', 'aims_junction','-b', 'aims_bottom','-b', 'aims_ss', '-b', 'aims_other', '-ot', self.Lgraph_label_basins.fullPath()]
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
                 '-i', self.Lwhite_mesh.fullPath(),
                 '-b', volumeGraphLabelBasins.fullPath(),
                 '-lb',self.Lgraph_label_basins.fullPath(),
                 '-ls',self.file_correspondance_constraint.fullPath(),
                 '-m', 3,
                 '-t', 2,
                 '-st', self.constraint_weight,
                 '-o', self.Lwhite_sulcalines.fullPath(),
                 '-si', self.side,
                 '-sb', self.basin_min_size,
                 '-cv', constraintValue
                 ]
      apply( context.system, sulcalines )

      context.write('Done')
            
      