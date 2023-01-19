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

from __future__ import print_function

from __future__ import absolute_import
import numpy as np
from scipy import sparse
from soma import aims
from soma import aimsalgo
from brainvisa.cortical_surface.surface_tools import basic_tools as basicTls
from brainvisa.cortical_surface.surface_tools import PDE_tools as pdeTls
from six.moves import range
####################################################################
#
# ensure the texture corresponding to the value tex_val has only one connex component with simple boundary
#
####################################################################
def textureTopologicalCorrection(mesh, atex, tex_val, background_val=0, neigh=None):
    tex_val_indices = np.where(atex == tex_val)[0]
    if not tex_val_indices.size:
        print('no value ' + str(tex_val) + ' in the input texture!!')
        return tuple()
    else:
        print(str(tex_val_indices.size) + ' vertices have the texture value ' + str(tex_val))
        if neigh is None:
            neigh = aims.SurfaceManip.surfaceNeighbours(mesh)
        labels = np.unique(atex)
        labels = labels.tolist()
        #----------------begin:: ensure a single connex component
        tex2 = aims.TimeTexture_S16()
        tex2[0].reserve(atex.size)
        for i in range(atex.size):
            if i in tex_val_indices:
                tex2[0].push_back(tex_val)
            else:
                tex2[0].push_back(0)
        conn_tex = aimsalgo.AimsMeshLabelConnectedComponent(mesh, tex2, 0, 0)  # aimsalgo.AimsMeshLabelConnectedComponent(mesh, tex2,1,0)
#        aconn_tex = conn_tex[0].arraydata()
        conn_labels = set(conn_tex[0].arraydata())
        conn_labels = conn_labels.difference(set([-1]))  # value -1 corresponds to background
        if len(conn_labels) > 1:
            print(len(conn_labels), " connex component(s) for this value in the texture, keeping only the largest component")
            lab_val_indices_nb = [len(np.where(conn_tex[0].arraydata() == lab)[0]) for lab in conn_labels]#.difference(set([-1]))]
            #print(lab_val_indices_nb#conn_labels#.difference(set([-1])),lab_val_indices_nb)
            index_val_largest_comp = lab_val_indices_nb.index(max(lab_val_indices_nb))
          #  index_val_largest_comp = index_val_largest_comp.tolist()
            adeleted_conn_comp_indices = np.where(conn_tex[0].arraydata() == -1)[0]
            deleted_conn_comp_indices = adeleted_conn_comp_indices.tolist()
            for conn_labels_index in set(range(len(conn_labels))) - set([index_val_largest_comp]):
                    deleted_conn_comp_indices.extend(np.where(conn_tex[0].arraydata() == list(conn_labels)[conn_labels_index])[0])
            atex[deleted_conn_comp_indices] = background_val
#            conn_tex.erase()
#            for i in range(atex.size):
#                if i in deleted_conn_comp_indices:
#                    conn_tex[0].push_back(0)
#                else:
#                    conn_tex[0].push_back(atex[i])
        #----------------end:: ensure a single connex component
        #----------------begin:: fill any hole in this single connex component
        boundary = basicTls.textureBoundary(mesh, atex, tex_val, neigh)
        if len(boundary) > 1:
            print("filling holes in the largest connex component")
#            ws=aims.Writer()
#            tex_out = aims.TimeTexture_S16()
#            i=0
            while len(boundary) > 1:
                print('length of boundaries : ',[len(bound) for bound in boundary])
#                i+=1
#                nm= 's12158_Rhippo_cleaned'+str(i)+'.tex'
#                tex_out[0].assign(atex)
#                ws.write(tex_out, nm)
                for l in range(len(boundary) - 1):
                    atex = dilateTexture(atex, neigh, boundary[l])
                boundary = basicTls.textureBoundary(mesh, atex, tex_val, neigh)
        print('length of boundaries : ',[len(bound) for bound in boundary])
        #----------------end:: fill any hole in this single connex component
        #----------------begin:: ensure that the boundary do not contain the 3 vertices of a triangle
        print('cleaning the longest boundary of the texture')
        atex, boundary = cleanTextureBoundary(mesh, atex, tex_val, boundary[-1], background_val, neigh)
        #----------------end:: ensure that the boundary do not contain the 3 vertices of a triangle

    return (atex, boundary)

####################################################################
# ensure that the boundary do not contain the 3 vertices of a triangle
####################################################################
def cleanTextureBoundary(mesh, tex, tex_val, bound, background_val=0, neigh=None):
    if neigh is None:
        neigh = aims.SurfaceManip.surfaceNeighbours(mesh)
    '''check if the 3 vertices of any polygon is on the boundary
    if it is the case, the vertex that can be removed is selected and removed'''    
    poly = np.array(mesh.polygon())
    I = basicTls.ismember(poly, bound)
    poly_set = poly[I[:, 0] & I[:, 1] & I[:, 2], :]
    if poly_set.shape[0]==0:
        return (tex, basicTls.textureBoundary(mesh, tex, tex_val, neigh))
    else:
        while poly_set.shape[0]>0:
            print('nb triangles on the boundary ',poly_set.shape[0])
            pts_to_remove = []
            for pb_poly in poly_set:
                ind_V1 = bound.index(pb_poly[0])
                ind_V2 = bound.index(pb_poly[1])
                ind_V3 = bound.index(pb_poly[2])
                list_inds = [ind_V1, ind_V2, ind_V3]
                m = min(list_inds)
                M = max(list_inds)
                pts_to_remove.append( bound[list(set(list_inds).difference(set([m,M])))[0]] )
#             print('pts_to_remove ', np.array(pts_to_remove))
#             print('poly_set ',poly_set)
#             print('bound ', bound)
            if len(pts_to_remove) > 0:
                tex[np.array(pts_to_remove)] = background_val
            boundary = basicTls.textureBoundary(mesh, tex, tex_val, neigh)
            bound = boundary[-1]
            I = basicTls.ismember(poly, bound)
            poly_set = poly[I[:, 0] & I[:, 1] & I[:, 2], :]
        if len(boundary) > 1:
            print('complex boundary after cleanTextureBoundary')
        return (tex, boundary)


def dilateTexture(tex, neigh, inds):
    for i in inds:
        tex[np.array(neigh[i].list())] = tex[i]
    return tex

def TextureExtrema(mesh, atex, neigh=None):
    if neigh is None:
        neigh = aims.SurfaceManip.surfaceNeighbours(mesh)
    extrema = np.zeros_like(atex)
    for v,ne_v in enumerate(neigh):
        mi = np.min(atex[ne_v.list()])
        ma = np.max(atex[ne_v.list()])
        if atex[v] < mi:
            extrema[v] = -1
        elif atex[v] > ma:
            extrema[v] = 1
    return extrema
