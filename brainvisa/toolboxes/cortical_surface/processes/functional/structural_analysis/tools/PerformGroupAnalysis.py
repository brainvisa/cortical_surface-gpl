from brainvisa.processes import *
import shfjGlobals     

name = '3 - Perform Group Analysis'
userLevel = 2

signature = Signature(
  'primalsketches', ListOf(ReadDiskItem('Primal Sketch', 'Graph and data')),
  'labeled_primalsketches', ListOf(WriteDiskItem('Primal Sketch', 'Graph and data')),
  'ddx1', Float(),
  'ddx2', Float(),
  'simweight', Float(),
  'datadrivenweight', Float(),
  'intrapsweight', Float(),
  'lsweight', Float(),
  'ddh', Float(),
  'run', Boolean(),
  'save', Boolean(),
  'energypath', ReadDiskItem('Text File', 'Text File'),
  'recuitpath', ReadDiskItem('Text File', 'Text File'))


def initialization( self ):
  self.setOptional('energypath', 'recuitpath')
  self.ddx1 = 8
  self.ddx2 = 4
  self.simweight = 1.0
  self.datadrivenweight = 0.8
  self.intrapsweight = 4.0
  self.lsweight = 1.0
  self.ddh = 0.0001
  self.run = True
  self.save = True
  
     

def execution( self, context ):

  templist=''
  tempout=''
  for prims in self.primalsketches :
      templist += '|' + str(prims)
  for primsout in self.labeled_primalsketches :
      tempout += '|' + str(primsout)
  sketchlist=templist[1:]
  sketchout=tempout[1:]
  
  

  call_list = ['surfStructuralAnalysis',
  '-p', str(sketchlist),
  '-o', str(sketchout),
  '--ddx1', self.ddx1,
  '--ddx2', self.ddx2,
  '--simw', self.simweight,
  '--ddw', self.datadrivenweight,
  '--ipsw', self.intrapsweight,
  '--lsw', self.lsweight,
  '--ddh', self.ddh,  
  '--run', int(self.run),
  '--save', int(self.save)]


  if self.energypath != "":
    call_list += ['--energypath', self.energypath]

  if self.recuitpath != "":
    call_list += ['--recuitpath', self.recuitpath]
  
  context.write("Starting group analysis (see BrainVISA log for more details...)")
  apply(context.system, call_list)
  context.write("Finished")