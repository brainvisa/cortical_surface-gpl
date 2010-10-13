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

from neuroProcesses import *
import shfjGlobals     

name = 'Create Surface-Based Functional Data'
userLevel = 0

signature = Signature(  'white_mesh', ReadDiskItem( 'Hemisphere White Mesh', 'BrainVISA mesh formats' ), 
        'fmri_data', ReadDiskItem('4D Volume', 'BrainVISA volume formats'),
        'timeserie_texture', WriteDiskItem('Functional Time Texture', 'Texture'),
  )

def getResolutionX(self,data):
  voxel_size  = data['voxel_size']
  return voxel_size[0]

def getResolutionY(self,data):

    voxel_size  = data['voxel_size']
    return voxel_size[1]

def getResolutionZ(self,data):
    voxel_size  = data['voxel_size']
    return voxel_size[2]

def getResolution( self, data):
  return [self.getResolutionX(data), self.getResolutionY(data), self.getResolutionZ(data)]


def getName(self,data):
  ##print len(self.fmri_data)
  ##if (len(self.functional_volumes) > 0 and self.intmesh is not None):
    #dial = defaultContext().dialog( 1, 'Enter a short name for the sequence to be projected - it will be included in the filenames :', Signature( 'param', String() ), _t_( 'OK' ), _t_( 'Cancel' ) )
    #dial.setValue( 'param', "" )
    #r = dial.call()
    #if r == 0:
      #v = dial.getValue( 'param' )
      #print v
      #result = None
    attributes = self.white_mesh.hierarchyAttributes()
    attributes[ 'volume' ] = 'fMRI'
    result = self.signature[ 'timeserie_texture' ].findValue( attributes )
    return result
    
    #return None




def initialization( self ):
  
  
  
    from neuroHierarchy import databases
    from copy import copy
    #self.signature['white_mesh']['volume'] = 'fMRI'

    eNode = SerialExecutionNode( self.name, parameterized=self )
    
    eNode.addChild( 'Average', ProcessExecutionNode( 'AverageVolumes', optional = 1 ) )
    eNode.addChild( 'Registration', ProcessExecutionNode( 'Register3DMutualInformation', optional = 1 ) )

    eNode.addChild( 'Kernels', ProcessExecutionNode( 'CreateKernels', optional = 1 ) )
    eNode.addChild( 'Projection', ProcessExecutionNode( 'ProjectionUsingKernels', optional = 1 ) )
    #eNode.addChild( 'FusionTextures', ProcessExecutionNode( 'TexturesAsTimeTexture', optional = 1 ) )

    eNode.addLink('Average.input','fmri_data')
    eNode.addLink('Registration.source_image','white_mesh')
    eNode.addLink('Registration.reference_image','Average.output')
    
    eNode.addLink('Kernels.intmesh','white_mesh')
    eNode.addLink('Kernels.resolution','fmri_data',self.getResolution)
    #eNode.addLink('Kernels.resolutionY','functional_volumes',self.getResolutionY)
    #eNode.addLink('Kernels.resolutionZ','functional_volumes',self.getResolutionZ)
    eNode.addLink('timeserie_texture', 'fmri_data', self.getName)
    eNode.addLink('timeserie_texture', 'white_mesh', self.getName)


    #for k in  eNode.childrenNames():
      #print k

    
    outAver = defaultContext().temporary('Nifti-1 image')

    eNode.Average.setValue('output',outAver )
    signat = copy( eNode.Registration.signature )
    signat[ 'source_image' ] = ReadDiskItem( 'T1 MRI Bias Corrected', 'BrainVISA volume formats' )
    signat[ 'source_to_reference' ] = WriteDiskItem( 'Mean Functional Volume To Anatomy Transformation', 'Transformation matrix' )
    signat[ 'reference_to_source' ] = WriteDiskItem( 'Anatomy To Mean Functional Volume Transformation', 'Transformation matrix')
    eNode.Registration.changeSignature( signat )
    eNode.addLink('Projection.registration_transformation','Registration.source_to_reference')
    eNode.addLink('Registration.source_to_reference', 'white_mesh')
    eNode.addLink('Registration.reference_to_source', 'white_mesh')
    eNode.addLink('Projection.white_mesh','white_mesh')
    eNode.addLink('Projection.fMRI_3D_data','fmri_data')
    eNode.Projection.removeLink('kernels', 'white_mesh')
    eNode.addLink('Projection.kernels','Kernels.output')
    eNode.addLink('Projection.fMRI_3D_data','fmri_data')
    eNode.addLink('Projection.fMRI_surface_data','timeserie_texture')


    #signat = copy( eNode.FusionTextures.signature )
    #signat[ 'output' ] = WriteDiskItem( 'Functional Time Texture', 'Texture' )
    #eNode.FusionTextures.changeSignature( signat )
    #eNode.addLink('FusionTextures.input','Projection.projection_textures')
    #eNode.addLink('FusionTextures.output','timetexture')


    self.setExecutionNode( eNode )


