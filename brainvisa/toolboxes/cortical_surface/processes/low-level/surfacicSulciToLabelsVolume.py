# -*- coding: iso-8859-1 -*-

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
import string

name = 'Surfacic Sulci To Labels Volume'
userLevel = 2

def validation():
  ver = [ int(x) for x in neuroConfig.versionString().split('.') ]
  if ver[0] < 3 or ( ver[0] == 3 and ver[1] < 1 ):
    raise ValidationError( 'BrainVisa 3.1 or higher is needed' )
  from soma import aims
  from soma import aimsalgo

signature = Signature(
  'patches_sulci_graph', ReadDiskItem( 'Graph', 'Graph' ),
  'mri_volume', ReadDiskItem( '3D Volume', shfjGlobals.aimsVolumeFormats ),
  'labels_volume', WriteDiskItem( 'Label Volume',
    shfjGlobals.aimsWriteVolumeFormats ),
  'input_translation', ReadDiskItem( 'log file', 'text file' ),
  'output_translation', WriteDiskItem( 'log file', 'text file' ),
  'background_labels', ListOf( String() ),
  'dilation_size', Float(),
  'erosion_size', Float(),
  )


def initialization( self ):
  self.setOptional( 'input_translation' )
  self.setOptional( 'output_translation' )
  self.background_labels = [ 'Bkgr_l', 'Bkgr_r', 'brain', 'unknown' ]
  self.dilation_size = 2
  self.erosion_size = 1
  self.linkParameters( 'labels_volume', 'patches_sulci_graph' )
  self.linkParameters( 'mri_volume', 'patches_sulci_graph' )
  self.linkParameters( 'output_translation', 'input_translation' )


def execution( self, context ):
  import aims, aimsalgosip
  labelsmap = {}
  indexmap = {}
  usedindex = []
  if self.input_translation:
    if os.path.exists( self.input_translation.fullPath() ):
      f = open( self.input_translation.fullPath() )
      for line in f.xreadlines():
        el = line.split()
        if len( el ) >= 2:
          labelsmap[ string.join( el[:-1] ) ] = int( el[-1] )
          usedindex.append( int( el[-1] ) )
      f.close()
      del f
    else:
      context.warning( 'input_translation file', self.input_translation,
        'does not exist' )

  r = aims.Reader()
  graph = r.read( self.patches_sulci_graph.fullPath() )
  meshatt = 'aims_patch'
  if self.mri_volume:
    atts = shfjGlobals.aimsVolumeAttributes( self.mri_volume )
    dims = atts[ 'volume_dimension' ]
    vs = atts[ 'voxel_size' ]
  else:
    context.warning( 'mri_volume not specified: taking dimensions and '
      'voxel size from graph, which may be inaccurate or wrong' )
    bM = list( aims.PFObject( graph )[ 'boundingbox_max' ] )
    bm = list( aims.PFObject( graph )[ 'boundingbox_min' ] )
    dims = [ x-y+1 for x,y in zip( bM, bm ) ]
    vs = graph[ 'voxel_size' ]
    vs = [ x.getScalar() for x in vs ]
  context.write( 'dims:', dims )
  context.write( 'vs:', vs )
  vol = aims.Volume_S16( *dims[:3] )
  vol.getPropertySet()[ 'voxel_size' ] = vs
  arrvol = vol.arraydata()
  arrvol.fill( 0 )
  alvol = aims.AimsData_S16( *( dims[:3] + [ 1, 1 ] ) )
  lvol = alvol.volume()
  lvol.getPropertySet()[ 'voxel_size' ] = vs
  arrlvol = lvol.arraydata()

  for v in graph.vertices():
    try:
      label = v[ 'name' ].getString()
      if label not in self.background_labels:
        context.write( 'processing node', label, '...' )
        mesh = aims.AimsSurfaceTriangle.fromObject( v[meshatt].get() )
        if label in labelsmap:
          oindex = labelsmap[ label ]
        else:
          oindex = 1
          while oindex in usedindex:
            oindex += 1
          usedindex.append( oindex )
          labelsmap[ label ] = oindex
        arrlvol.fill( 0 )
        for p in mesh.vertex():
          lvol.setValue( 32767,
            *[ int( round(x/s) ) for x,s in zip(p,vs) ] )
        avol = aimsalgosip.AimsMorphoChamferDilation( alvol,
          self.dilation_size )
        avol = aimsalgosip.AimsMorphoErosion( avol, self.erosion_size )
        arrvol[ avol.volume().arraydata()[0:1,1:-1,1:-1,1:-1]==32767 ] = oindex
    except Exception, e:
      context.warning( 'node failed:', e )
      pass

  w = aims.Writer()
  w.write( vol, self.labels_volume.fullPath() )
  if self.output_translation:
    f = open( self.output_translation.fullPath(), 'w' )
    for x, y in labelsmap.items():
      print >> f, x, y
