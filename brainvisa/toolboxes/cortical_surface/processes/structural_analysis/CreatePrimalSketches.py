from neuroProcesses import *
import shfjGlobals     

name = 'Create Primal Sketches'
userLevel = 2

signature = Signature(
    'texture', ReadDiskItem( 'Texture','Texture'),
    'white', ReadDiskItem( 'Hemisphere White Mesh', 'MESH mesh'),
    'primalsketch', WriteDiskItem( 'Primal Sketch', 'Graph' ),
    'tMin', Float(),
    'tMax', Float(),
    'whiteAux', ReadDiskItem( 'Hemisphere White Mesh', 'MESH mesh'),
    'filterout', Float(),
    'intersectioncriterium', Integer(),
    'latitude', ReadDiskItem( 'Latitude coordinate texture','Texture'),
    'longitude', ReadDiskItem('Longitude coordinate texture','Texture')
    )

def initialization( self ):
     self.setOptional('whiteAux','filterout','intersectioncriterium','latitude','longitude')
     self.linkParameters( 'primalsketch', 'texture' )
     self.linkParameters( 'white', 'texture' )
     self.linkParameters( 'whiteAux', 'white' )
     self.linkParameters( 'latitude', 'texture' )
     self.linkParameters( 'longitude', 'texture' )
     self.tMin = 1.0
     self.tMax = 64.0
     self.filterout = 2.0
     self.intersectioncriterium = 10

def execution( self, context ):
     scales=context.temporary( 'Texture' )
     blobs=context.temporary( 'Texture' )

     call_list = [ 'AimsTexturePrimalSketch',
                   '-t', self.texture,
                   '-o', self.primalsketch,
       '-m', self.white,
                   '-os', scales,
                   '-ob', blobs,
                   '-t1', self.tMin,
                   '-t2', self.tMax,
       ]
     if (self.whiteAux is not None):
       call_list += ['-mX', self.whiteAux]
     s = self.white.get( 'subject')
     assert(s!=None)
     call_list += ['-sj', s]
         
     if (self.filterout is not ""):
       call_list += ['-f', self.filterout]
     if (self.latitude is not None):
       call_list += ['-l', self.latitude]
     if (self.longitude is not None):
       call_list += ['-L', self.longitude]
     if (self.intersectioncriterium is not ""):
       call_list += ['-iP', self.intersectioncriterium]


     context.write('Starting primal sketch computation')
     apply( context.system, call_list )
     context.write('Finished')

