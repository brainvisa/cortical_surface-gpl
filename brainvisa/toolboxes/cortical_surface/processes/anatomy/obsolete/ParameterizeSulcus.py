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
from numpy import *
from six.moves import range

name = 'Sulcus Parameterization'

userLevel = 0

signature = Signature(
    'graph', ReadDiskItem( 'Cortical folds graph', 'Graph' ),
    'mri', ReadDiskItem( 'T1 MRI Bias Corrected', 'Aims readable volume formats' ),
    'label_attributes', Choice( 'label', 'name' ),
    'label_values', String(),
    'orientation', Choice( 'Top->Bottom', 'Front->Back' ),
    'sulcus_mesh', WriteDiskItem( 'Sulcus mesh', 'Aims mesh formats' ),
    'texture_param1', WriteDiskItem( 'Sulcus x coordinate texture', 'Texture' ),
    'texture_param2', WriteDiskItem( 'Sulcus y coordinate texture', 'Texture' ),
    'coordinates_grid', WriteDiskItem( 'Sulcus coordinate grid mesh', 'Aims mesh formats' ),
    'depth_profile', WriteDiskItem( 'Sulcus depth profile', 'Text file' ),
    'dilation', Float(),
    'morpho_offset', Choice( '0.0', '0.5', '0.8', '1.0' ),
    'meshDecimation', Boolean(),
    'deltaT', Float(),
    'stop', Float(),
)

def initialization( self ):
     self.linkParameters( 'mri', 'graph' )
     self.linkParameters( 'sulcus_mesh', 'graph' )
     self.linkParameters( 'texture_param1', 'sulcus_mesh' )
     self.linkParameters( 'texture_param2', 'sulcus_mesh' )
     self.linkParameters( 'coordinates_grid', 'sulcus_mesh' )
     self.linkParameters( 'depth_profile', 'sulcus_mesh' )
     self.label_attributes = 'name'
     self.dilation = 1.0
     self.morpho_offset = '0.0'
     self.deltaT = 0.05
     self.stop = 0.000003
     self.meshDecimation=0


def execution( self, context ):
     sulcusIm=context.temporary( 'GIS image' )
     closedIm=context.temporary( 'GIS image' )
     hullIm=context.temporary(  'GIS image' )
#     dilHull=context.temporary(  'GIS image' )
     bottomIm=context.temporary(  'GIS image' )
#     dilBottom=context.temporary(  'GIS image' )
     simplesurf=context.temporary( 'GIS image' )
     dilatedIm=context.temporary(  'GIS image' )
     isoL=context.temporary( 'MESH mesh')
     meshNonDec=context.temporary( 'MESH mesh')
     
     distToPlan=context.temporary( 'Texture' )
     
     transform=''

     if (self.orientation=='Top->Bottom'):
          orient=0
     else :
          orient=1

     context.write('Extracting sulcus and buckets from graph')

     context.runProcess( 'Create Sulcus Label Volume',
                         graph = self.graph,
                         mri = self.mri,
                         sulci =  sulcusIm.fullPath(),
                         simple_surface = simplesurf.fullPath(),
                         compress= 'No',
                         bucket= 'Sulci',
                         label_attributes = self.label_attributes,
                         label_values = self.label_values,
                         node_edge_types='All' )
     context.runProcess( 'Create Sulcus Label Volume',
                         graph = self.graph,
                         mri = self.mri,
                         sulci =  sulcusIm.fullPath(),
                         simple_surface = simplesurf.fullPath(),
                         bottom = bottomIm.fullPath(),
                         compress= 'No',
                         bucket= 'Bottoms',
                         label_attributes = self.label_attributes,
                         label_values = self.label_values,
                         node_edge_types='All' )
     context.runProcess( 'Create Sulcus Label Volume',
                         graph = self.graph,
                         mri = self.mri,
                         sulci =  sulcusIm.fullPath(),
                         simple_surface = simplesurf.fullPath(),
                         hull_junction = hullIm.fullPath(),
                         compress= 'No',
                         bucket= 'Junctions with brain hull',
                         label_attributes = self.label_attributes,
                         label_values = self.label_values,
                         node_edge_types='All' )

     dilating = ['AimsMorphoMath', '-m', 'dil',
                 '-i', sulcusIm.fullPath(),
                 '-o', dilatedIm.fullPath(),
                 '-r', self.dilation ]

     context.system(*dilating)

#     context.write('Dilating ridges with 3.0')

#     dilatingH = [ 'AimsDilation',
#                 '-i', hullIm.fullPath(),
#                 '-o', dilHull.fullPath(),
#                 '-e', 3.0 ]

#     apply( context.system, dilatingH )

#     dilatingB = [ 'AimsDilation',
#                 '-i', bottomIm.fullPath(),
#                 '-o', dilBottom.fullPath(),
#                 '-e', 3.0 ]
#
#     apply( context.system, dilatingB )

     closing = [ 'AimsMorphoMath', '-m', 'clo',
                 '-i', dilatedIm.fullPath(),
                 '-o', closedIm.fullPath(),
                 '-r', 2 ]

     context.system(*closing)

     context.write('Remeshing sulcus')


     meshing = [ 'AimsMesh',
                 '-i', closedIm.fullPath(),
                 '-o', self.sulcus_mesh.fullPath(),
                 '-l', '32767',
                 '--smooth',
                 '--smoothIt', '20' ]
     context.system(*meshing)

     test=self.sulcus_mesh.fullName()
     sulcusMname=test + '_32767_0.mesh'

     if (self.meshDecimation==0) :
          shelltools.mv(sulcusMname, self.sulcus_mesh.fullPath())
          shelltools.mv(sulcusMname + '.minf', self.sulcus_mesh.fullPath() + '.minf')
     else :
          decimating = ['AimsMeshDecimation',
                    '-i', sulcusMname,
                    '-o', self.sulcus_mesh.fullPath() ]
          context.system(*decimating)


     context.write('Starting parameterisation')

     parameterising = [ 'AimsParameterizeSulcus',
                        '-i', self.sulcus_mesh.fullPath(),
                        '-b', bottomIm.fullPath(),
                        '-t', hullIm.fullPath(),
                        '-o', orient,
                        '-ox', self.texture_param1.fullPath(),
                        '-oy', self.texture_param2.fullPath(),
                        '-mo', self.morpho_offset,
                        '-d', self.deltaT,
                        '-s', self.stop ]

     context.system(*parameterising)

     context.write('Computing coordinate grid')

     i=0;
     iso = [ 'AimsMeshIsoLine',
              '-i', self.sulcus_mesh.fullPath(),
              '-t', self.texture_param2.fullPath(),
              '-o', isoL.fullPath(),
              '-v', i ]
     conc = [ 'AimsZCat',
              '-i', isoL.fullPath(),
              '-o', self.coordinates_grid.fullPath() ]
     context.system(*iso)
     context.system(*conc)
     iso = [ 'AimsMeshIsoLine',
              '-i', self.sulcus_mesh.fullPath(),
              '-t', self.texture_param1.fullPath(),
              '-o', isoL.fullPath(),
              '-v', i ]
     conc = [ 'AimsZCat',
              '-i', isoL.fullPath(),
              '-o', self.coordinates_grid.fullPath() ]
     context.system(*iso)
     context.system(*conc)
     i=i+10
     while (i<=100):
          iso = [ 'AimsMeshIsoLine',
                   '-i', self.sulcus_mesh.fullPath(),
                   '-t', self.texture_param2.fullPath(),
                   '-o', isoL.fullPath(),
                   '-v', i ]
          conc = [ 'AimsZCat',
                   '-i', isoL.fullPath(), self.coordinates_grid.fullPath(),
                   '-o', self.coordinates_grid.fullPath() ]
          context.system(*iso)
          context.system(*conc)
          i=i+10
     i=10
     while (i<=200):
          iso = [ 'AimsMeshIsoLine',
                   '-i', self.sulcus_mesh.fullPath(),
                   '-t', self.texture_param1.fullPath(),
                   '-o', isoL.fullPath(),
                   '-v', i ]
          conc = [ 'AimsZCat',
                   '-i', isoL.fullPath(), self.coordinates_grid.fullPath(),
                   '-o', self.coordinates_grid.fullPath() ]
          context.system(*iso)
          context.system(*conc)
          i=i+10
          
          
     read=aims.Reader()
     sulc=read.read(self.sulcus_mesh.fullPath())
     mesh=array(sulc.vertex())
     bary=mean(mesh, axis=0)
     mesh=mesh-bary
     tmesh=mesh.transpose()
     coord=dot(tmesh,mesh)
     val,vect=linalg.eig(coord)
     i=argmin(val)
     k=argmax(val)
     for t in range(3):
          if (t!=i) and (t!=k):
               j=t
        
     u1=vect[:,i]
     u2=vect[:,j]
     u3=vect[:,k]

     texturex=aims.TimeTexture_FLOAT() 
     nn= mesh.shape[0]
     for i in range(nn):
          texturex.push_back(dot(mesh[i], -u1))
    
     w=aims.Writer()
     w.write(texturex, distToPlan.fullPath())
          
     context.write('Computing depth profile')
     depth = [ 'AimsSulcusNormalizeDepthProfile',
               '-m', self.sulcus_mesh.fullPath(),
               '-x', self.texture_param1.fullPath(),
               '-y', self.texture_param2.fullPath(),
               '-d', distToPlan.fullPath(),
               '-o', self.depth_profile.fullPath() ]
     context.system(*depth)

     context.write('Finished')

