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


def validation():
    try:
        from brainvisa.cortical_surface.surface_tools import PDE_tools as pdeTls
        from brainvisa.cortical_surface.surface_tools import mesh_watershed as watershed
    except:
        raise ValidationError( 'brainvisa.cortical_surface.surface_tools.mesh_watershed module can not be imported.' )
  

from brainvisa.processes import *
try:
    from brainvisa.cortical_surface.surface_tools import PDE_tools as pdeTls
    from brainvisa.cortical_surface.surface_tools import mesh_watershed as watershed
except:
    pass

from soma import aims
import numpy as np

name = 'Mesh Watershed'
userLevel = 0

# Argument declaration
signature = Signature(
    'input_mesh',ReadDiskItem( 'Hemisphere White Mesh' , 'Aims mesh formats' ),
    'DPF_texture',ReadDiskItem( 'DPF texture',  'Aims texture formats' ),
    'mask_texture',ReadDiskItem( 'Hippocampus pole texture','Aims Texture formats' ),
    'thresh_ridge', Float(),
    'thresh_dist', Float(),
    'thresh_area', Float(),
    'pits_texture',WriteDiskItem( 'pits texture',  'Aims texture formats' ),
    'noisypits_texture',WriteDiskItem( 'noisy pits texture',  'Aims texture formats' ),
    'ridges_texture',WriteDiskItem( 'ridges texture',  'Aims texture formats' ),
    'basins_texture',WriteDiskItem( 'basins texture',  'Aims texture formats' ),
#    'areas_texture',WriteDiskItem( ' texture',  'Aims texture formats' ),
)


# Default values
def initialization( self ):
    self.linkParameters( 'DPF_texture', 'input_mesh' )
    self.linkParameters( 'pits_texture', 'input_mesh' )
    self.linkParameters( 'noisypits_texture', 'input_mesh' )
    self.linkParameters( 'ridges_texture', 'input_mesh' )
    self.linkParameters( 'basins_texture', 'input_mesh' )
    self.thresh_dist=20
    self.thresh_ridge=1.5
    self.thresh_area=50
    self.setOptional('mask_texture')

def execution( self, context ):
    re = aims.Reader()
    ws = aims.Writer()


    #labelVoronoi=watershed.voronoiArea(mesh_file)
    mesh = re.read( self.input_mesh.fullPath() )
    vert_area = pdeTls.vertexVoronoi( mesh )
    #diff = np.max(np.array(labelVoronoi[0])-vert_area)
    #print diff
    #area=np.sum(labelVoronoi[0])

    depthTex = re.read( self.DPF_texture.fullPath() )
    depthArray = np.array( depthTex[0] )

    ##############################
    #depthArray = -depthArray
    ###################################

    if self.mask_texture is not None:
        maskTex = re.read( self.mask_texture.fullPath() )
        mask = np.array( maskTex[0] )
    else:
        mask = np.zeros( depthArray.shape )

    ## Watershed
    # first step : merging online
    labels_1, pitsKept_1, pitsRemoved_1, ridgePoints, parent_1=watershed.watershed(mesh, vert_area, depthArray, mask, self.thresh_dist, self.thresh_ridge)

    # second step : merging offline
    labels, infoBasins, pitsKept, pitsRemoved_2, parent=watershed.areaFiltering(mesh, vert_area, labels_1, pitsKept_1, parent_1, self.thresh_area)
    pitsRemoved=pitsRemoved_1+pitsRemoved_2

    ## Saving
    # texture of basins
    labelsTexture = aims.TimeTexture_FLOAT(1,len(labels))
    labelsTexture[0].assign(labels)
    ws.write( labelsTexture, self.basins_texture.fullPath() )
    # texture of pits
    PITS=np.zeros((len(labels),1))
    for pit in pitsKept:
        PITS[pit[0]]=1
    pitsTexture = aims.TimeTexture_FLOAT(1,len(labels))
    pitsTexture[0].assign(PITS)
    ws.write( pitsTexture, self.pits_texture.fullPath() )
    # texture of noisy pits
    NOISYPITS = np.zeros((len(labels),1))
    for pit in pitsRemoved:
        NOISYPITS[pit[0]]=1
    noisypitsTexture = aims.TimeTexture_FLOAT(1,len(labels))
    noisypitsTexture[0].assign(NOISYPITS)
    ws.write( noisypitsTexture, self.noisypits_texture.fullPath() )
    # texture of ridges
    RIDGES=np.zeros((len(labels),1))
    for ridge in ridgePoints:
        RIDGES[ridge[2]]=1
    ridgesTexture = aims.TimeTexture_FLOAT(1,len(labels))
    ridgesTexture[0].assign(RIDGES)
    ws.write( ridgesTexture, self.ridges_texture.fullPath() )
    # texture of basins areas
    ##AREAS=np.zeros((len(labels),1))
    ##for basin in infoBasins:
    ##    verts=np.where(labels==basin[0])[0]
    ##    for v in verts:
    ##        AREAS[v]=basin[3]
    ##areasTexture=aims.TimeTexture_FLOAT(1,len(labels))
    ##areasTexture[0].assign(AREAS)
    ##ws.write(areasTexture, areasTexture_file)
