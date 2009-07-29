# Copyright CEA and IFR 49 (2000-2005)
#
#  This software and supporting documentation were developed by
#      CEA/DSV/SHFJ and IFR 49
#      4 place du General Leclerc
#      91401 Orsay cedex
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

name = 'Simulate Activation Texture'
userLevel = 2

signature = Signature(
  'white_mesh', ReadDiskItem( 'Hemisphere White Mesh', 'MESH mesh' ),
  'output', WriteDiskItem('Texture', 'Texture'),
  'location', Float(),
  'intensity', Float(),
  'noise', Float(),
  'lissage', Float(),
  'focus1', Point3D(),
  'intensity1', Float(),
  'focus2', Point3D(),
  'intensity2', Float(),
  'focus3', Point3D(),
  'intensity3', Float(),
  'focus4', Point3D(),
  'intensity4', Float(),
  'focus5', Point3D(),
  'intensity5', Float()
  )
def initialization( self ):
  self.setOptional('focus1','focus2','focus3','focus4','focus5','intensity1', 'intensity2', 'intensity3', 'intensity4', 'intensity5')
  self.signature[ 'focus1' ].add3DLink( self, 'white_mesh' )
  self.signature[ 'focus2' ].add3DLink( self, 'white_mesh' )
  self.signature[ 'focus3' ].add3DLink( self, 'white_mesh' )
  self.signature[ 'focus4' ].add3DLink( self, 'white_mesh' )
  self.signature[ 'focus5' ].add3DLink( self, 'white_mesh' )
  
  self.noise = 6.0
  self.intensity = 0.1
  self.location = 10.0
  self.lissage = 16.0
  self.intensity1 = 5.0
  self.intensity2 = 5.0
  self.intensity3 = 5.0
  self.intensity4 = 5.0
  self.intensity5 = 5.0
  self.focus1 = [113.669, 69.4213, 29.1237]
  self.focus2 = [141.079, 81.536, 58.6003]
  self.focus3 = [129.447, 131.799, 40.0678]
  self.focus4 = [123.486, 49.0954, 62.2587]
  self.focus5 = [141.079, 81.536, 58.6003]

  
     

def execution( self, context ):
  from soma import aims
  import sys,os,string 
  f1 = []
  f1 = self.focus1
  f1.append(self.intensity1)
  f2 = []
  f2 = self.focus2
  f2.append(self.intensity2)
  f3 = []
  f3 = self.focus3
  f3.append(self.intensity3)
  f4 = []
  f4 = self.focus4
  f4.append(self.intensity4)
  f5 = []
  f5 = self.focus5
  f5.append(self.intensity5)
  process = ['surfTexActivationSimulation',
    '-o',  self.output,
    '-m', self.white_mesh,
    '-n', self.noise,
    '-i', self.intensity,
    '-l', self.location,
    '-s', self.lissage,
    '-f1x', f1[0], '-f1y', f1[1], '-f1z', f1[2], '-f1t', f1[3],
    '-f2x', f2[0], '-f2y', f2[1], '-f2z', f2[2], '-f2t', f2[3],
    '-f3x', f3[0], '-f3y', f3[1], '-f3z', f3[2], '-f3t', f3[3],
    '-f4x', f4[0], '-f4y', f4[1], '-f4z', f4[2], '-f4t', f4[3],
    '-f5x', f5[0], '-f5y', f5[1], '-f5z', f5[2], '-f5t', f5[3]]
  apply( context.system, process )
  reader = aims.Reader()
  object = reader.read( str(self.output) )
  print 'file:', object
  h = object.header()
  h['noise']=str(self.noise)
  h['location']=str(self.location)
  h['lissage']=str(self.lissage)
  h['intensity']=str(self.intensity)
  h['focus1']=str(f1[0]) + ' ' + str(f1[1]) + ' '  + str(f1[2]) + ' ' + str(f1[3])
  h['focus2']=str(f2[0]) + ' ' + str(f2[1]) + ' '  + str(f2[2]) + ' ' + str(f2[3])
  h['focus3']=str(f3[0]) + ' ' + str(f3[1]) + ' '  + str(f3[2]) + ' ' + str(f3[3])
  h['focus4']=str(f4[0]) + ' ' + str(f4[1]) + ' '  + str(f4[2]) + ' ' + str(f4[3])
  h['focus5']=str(f5[0]) + ' ' + str(f5[1]) + ' '  + str(f5[2]) + ' ' + str(f5[3])
  writer = aims.Writer()
  writer.write(object, str(self.output))
