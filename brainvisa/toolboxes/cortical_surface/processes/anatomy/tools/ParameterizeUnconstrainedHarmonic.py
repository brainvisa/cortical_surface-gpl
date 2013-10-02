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
    import brainvisa.cortical_surface.parameterization.mapping
  except:
    raise ValidationError( 'brainvisa.cortical_surface.parameterization.mapping module can not be imported.' )
  
from brainvisa.processes import *
from brainvisa.tools import aimsGlobals
from soma import aims
import numpy as np

try:
  from brainvisa.cortical_surface.parameterization import mapping as map#hipHop
  from brainvisa.cortical_surface.surface_tools import surface_tools as surfTls
  from brainvisa.cortical_surface.parameterization import sulcalLinesSet as slSet
  from brainvisa.cortical_surface.parameterization import model as md
except:
  pass

#from brainvisa import anatomist

name = 'Harmonic Intrinsic Parameterization (HIP)'

userLevel = 1

# def validation():
#     anatomist.validation()
    
signature = Signature(
                      
    'white_mesh',ReadDiskItem( 'Hemisphere White Mesh', aimsGlobals.aimsMeshFormats),
    'side', Choice('left', 'right'),
    'cingular_pole_texture',ReadDiskItem( 'Hippocampus pole texture', 'Texture'),
    'insular_pole_texture',ReadDiskItem( 'Insula pole texture', 'Texture'),
    'white_sulcalines',ReadDiskItem( 'hemisphere Sulcal Lines texture', 'Texture' ),                     
    'sulcus_labels',ReadDiskItem( 'Graph Label Translation', 'Text File'),
    'model_file',ReadDiskItem( 'HipHop Model', 'Text File'),
    'unfold_reversed_triangles', Choice('yes','no'),
    'nb_it_local_smoothing_for_unfolding', Integer(),
    'rectangular_mesh',WriteDiskItem( 'Rectangular flat mesh', aimsGlobals.aimsMeshFormats),
    'rectangular_white_sulcalines',WriteDiskItem( 'hemisphere Sulcal Lines Rectangular Flat texture', 'Texture' ),
    'boundary_texture',WriteDiskItem( 'Rectangular boundary texture', 'Texture'),
    'corresp_indices_texture',WriteDiskItem( 'Rectangular flat indices texture', 'Texture'),
    'white_mesh_parts',WriteDiskItem( 'White Mesh Parts', aimsGlobals.aimsMeshFormats)
)

def initialization( self ):
    def linkSide( proc, dummy ):
        if proc.white_mesh is not None:
            return proc.white_mesh.get( 'side' )
    self.linkParameters( 'side', 'white_mesh', linkSide )
    self.linkParameters( 'cingular_pole_texture', 'white_mesh')
    self.linkParameters( 'insular_pole_texture', 'white_mesh')
    self.linkParameters( 'rectangular_mesh','white_mesh' )
    self.linkParameters( 'boundary_texture','white_mesh')
    self.linkParameters( 'corresp_indices_texture','white_mesh')
    self.linkParameters( 'white_mesh_parts','white_mesh')
    self.linkParameters( 'white_sulcalines', 'white_mesh')
    self.linkParameters( 'rectangular_white_sulcalines', 'white_mesh')
    self.linkParameters( 'sulcus_labels', 'white_mesh')
    self.setOptional( 'model_file' )
    self.unfold_reversed_triangles = 'yes'
    self.nb_it_local_smoothing_for_unfolding = 100
    
def execution( self, context ):
  
    re = aims.Reader()
    ws = aims.Writer()
    context.write('Reading textures and mesh')
    cing_pole = re.read(self.cingular_pole_texture.fullPath())
    insula_pole = re.read(self.insular_pole_texture.fullPath())
    sulcal_lines = re.read(self.white_sulcalines.fullPath())
    mesh = re.read(self.white_mesh.fullPath())
    context.write('Reading sulcus-label correspondences file')
    sulc_labels_dict = surfTls.readSulcusLabelTranslationFile(self.sulcus_labels.fullPath())

    if self.model_file is not None:
        context.write('Reading model')
        model = md.Model().read(self.model_file.fullPath())
        for line in model.printArgs().splitlines():
            context.write(line)
    else:
        model = md.Model()
    length = model.right - model.left
    width = model.top - model.bottom

    context.write('HIP')
    (neoCortex_square, neoCortex_open_boundary, neocortex_indices, insula_indices, cingular_indices, insula_mesh, cingular_mesh, neoCortex_mesh) = map.hip(mesh, insula_pole[0].arraydata(), cing_pole[0].arraydata(), length, width)
    (nb_inward, inward) = map.invertedPolygon(neoCortex_square)
    vert = np.array(neoCortex_square.vertex())
    context.write('------------------number of vertices on folded triangles : '+str(nb_inward)+' => '+str(100.0 * nb_inward / vert.shape[0])+' %')
    if self.unfold_reversed_triangles == 'yes':
        poly = np.array(neoCortex_square.polygon())
        context.write('------------------unfolding reversed triangles')
        (neoCortex_square, nb_inward_evol, inward_evol) = map.solveInvertedPolygon(neoCortex_square, neoCortex_open_boundary, self.nb_it_local_smoothing_for_unfolding)
        vert = np.array(neoCortex_square.vertex())
        context.write('------------------evolution of the iterative unfolding : '+str(nb_inward_evol))
#         inward_tex = 'tmp.tex'
#         context.write('------------------writing inward tex in : '+inward_tex)
#         tmp_tex = np.zeros(len(neoCortex_square.vertex()))
# #        print np.unique(poly[inward, :])
#         tmp_tex[np.unique(poly[inward_evol[-1], :])] = 1
#         tex_unfold = aims.TimeTexture_S16()
#         tex_unfold[0].assign(tmp_tex)
#         ws.write(tex_unfold, inward_tex)
    context.write('Translating the barycenter of S.C. to 0')
    output_SL_tex = aims.TimeTexture(sulcal_lines)
    output_tex_tmp = sulcal_lines[0].arraydata()[neocortex_indices]
    tex_square_sulci = np.zeros(vert.shape[0], sulcal_lines[0].arraydata().dtype )
    tex_square_sulci[range( len(neocortex_indices) )] = output_tex_tmp
    for b in neoCortex_open_boundary:
        tex_square_sulci[b] = 0
    output_SL_tex[0].assign(tex_square_sulci)
    
    
    labels = np.unique(tex_square_sulci)
    labels = labels[labels!=0]
    context.write('found the following sulci in the texture :')
    context.write([sulc_labels_dict[lab] for lab in labels])
    context.write('associated to the following labels :')
    context.write(labels)
    full_sulci = slSet.SulcalLinesSet() 
    full_sulci.extractFromTexture(tex_square_sulci, neoCortex_square, sulc_labels_dict)
    SC_ind = full_sulci.names.index(('S.C._'+self.side))   
    SC_label = full_sulci.labels[SC_ind]
#    full_sulci.sulcalLines[SC_ind].printArgs()
    translation = -full_sulci.sulcalLines[SC_ind].barycenter[0]
    vert[:, 0] = vert[:, 0] + translation # * np.ones(vert.shape[0])
    neoCortex_square.vertex().assign([aims.Point3df(x) for x in vert])
    
    context.write('Writing meshes and textures')
    if self.side == 'right':
        poly = np.array(neoCortex_square.polygon())
        poly_tmp = poly.copy()
#        context.write(poly_tmp[0,:])
        poly_tmp[:,0] = poly[:,1]
        poly_tmp[:,1] = poly[:,0]
        pp = aims.vector_AimsVector_U32_3()
        for i in poly_tmp:
            pp.append(i)
#        context.write(np.array(pp)[0,:])
        neoCortex_square.polygon().assign(pp)
        neoCortex_square.updateNormals()
    ws.write( neoCortex_square, self.rectangular_mesh.fullPath() )
    ws.write(output_SL_tex, self.rectangular_white_sulcalines.fullPath())
    mesh_parts = aims.AimsTimeSurface_3()
    '''
    mesh_parts[0] = neoCortex
    mesh_parts[1] = insula
    mesh_parts[2] = cingular pole
    '''
    mesh_parts.vertex( 0 ).assign( neoCortex_mesh.vertex() )
    mesh_parts.normal( 0 ).assign( neoCortex_mesh.normal() )
    mesh_parts.polygon( 0 ).assign( neoCortex_mesh.polygon() )
    mesh_parts.vertex( 1 ).assign( insula_mesh.vertex() )
    mesh_parts.normal( 1 ).assign( insula_mesh.normal() )
    mesh_parts.polygon( 1 ).assign( insula_mesh.polygon() )
    mesh_parts.vertex( 2 ).assign( cingular_mesh.vertex() )
    mesh_parts.normal( 2 ).assign( cingular_mesh.normal() )
    mesh_parts.polygon( 2 ).assign( cingular_mesh.polygon() )
    ws.write(mesh_parts, self.white_mesh_parts.fullPath())

    '''boundaries (see mapping.path2Boundary for details:"
    boundary[0] == insula_boundary
    boundary[1] == neocortex_poles_path always from insula to cingular pole
    boundary[2] == cingular_boundary
    boundary[3] == new vertices always from insula to cingular pole
    '''
    tex_boundary = aims.TimeTexture_S16()
    for ind,bound in enumerate(neoCortex_open_boundary):
        tmp_tex = np.zeros(len(neoCortex_square.vertex()))
        tmp_tex[bound] = 1
        tex_boundary[ind].assign(tmp_tex)
    ws.write(tex_boundary, self.boundary_texture.fullPath())
    '''
    tex_corresp_indices contains the indices of the vertices in white_mesh for:
        neoCortex_square in time 0
        insula_indices in time 1
        cingular_indices in time 2
    '''
    tex_corresp_indices = aims.TimeTexture_S16()
    tmp_tex = np.zeros(len(mesh.vertex()))
    tmp_tex[neocortex_indices] = 1
    tex_corresp_indices[0].assign(tmp_tex)
    tmp_tex = np.zeros(len(mesh.vertex()))
    tmp_tex[insula_indices] = 1
    tex_corresp_indices[1].assign(tmp_tex)
    tmp_tex = np.zeros(len(mesh.vertex()))
    tmp_tex[cingular_indices] = 1
    tex_corresp_indices[2].assign(tmp_tex)
    ws.write(tex_corresp_indices, self.corresp_indices_texture.fullPath())

#     re = aims.Reader()
#     ws = aims.Writer()
#     context.write('Reading textures and mesh')
#     cing_pole = re.read(self.cingular_pole_texture.fullPath())
#     insula_pole = re.read(self.insular_pole_texture.fullPath())
#     texture_sulci = re.read(self.white_sulcalines.fullPath())
#     mesh = re.read(self.white_mesh.fullPath())
#     context.write('HIP-HOP')
#     context.write(np.unique(texture_sulci[0].arraydata()))
# #     from brainvisa.cortical_surface.parameterization import sulcalLinesSet as slSet
# #     full_sulci = slSet.SulcalLinesSet()
# #     full_sulci.extractFromTexture(texture_sulci[0].arraydata(), mesh)
# #     full_sulci.printArgs()
# 
#     
#     lon, lat = hipHop(mesh, insula_pole[0].arraydata(), cing_pole[0].arraydata(), texture_sulci[0].arraydata(), self.side)
#     context.write('Writing textures')
#     tex_lon = aims.TimeTexture_FLOAT()
#     tex_lon[0].assign(lon)
#     ws.write(tex_lon, self.longitude.fullPath())
#     tex_lat = aims.TimeTexture_FLOAT()
#     tex_lat[0].assign(lat)
#     ws.write(tex_lat, self.latitude.fullPath())

#    spherical_verts = sphericalMeshFromCoords(lat, lon, 50):
    
    context.write('Done')
            
      
