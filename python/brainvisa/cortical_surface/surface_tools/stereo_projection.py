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

from __future__ import absolute_import
import numpy as np
import os
from soma import aims


def stereo_projection(vertices, h=None):
    """
    compute the stereographic projection from the unit sphere (center = 0, radius = 1) onto the horizontal plane which 3rd coordinate is h of the vertices given
    :param vertices: Nx3 coordinates of the vertices to be projected onto the plane
    :param h: 3rd coordinate of the projection plane
    :return: Nx3 coordinates of the vertices projected onto the plane, their 3rd coordinate is thus equal to h
    """
    if h is None:
        h = -1
    for ind, vert in enumerate(vertices):
        vertices[ind, 0] = (-h + 1) * vert[0] / (1 - vert[2])
        vertices[ind, 1] = (-h + 1) * vert[1] / (1 - vert[2])
        vertices[ind, 2] = h
    return vertices


def inverse_stereo_projection(vertices, h=None):
    """
    compute the inverse stereograhic projection from an horizontal plane onto the unit sphere (center = 0, radius = 1)
    :param vertices: Nx3 vertices to be inverse  projected onto the sphere
    :param h: 3rd coordinate of the projection plane
    :return: Nx3 coordinates of the vertices onto the unit sphere
    """
    if h is None:
        h = vertices[0, 2]
    for ind, vert in enumerate(vertices):
        denom = ((1 - h) ** 2 + vert[0] ** 2 + vert[1] ** 2)
        vertices[ind, 2] = (-(1 - h) ** 2 + vert[0] ** 2 + vert[1] ** 2) / denom
        vertices[ind, 1] = 2 * (1 - h) * vert[1] / denom
        vertices[ind, 0] = 2 * (1 - h) * vert[0] / denom
    return vertices


def rotation(src, targ):
    """
    compute the rotation between the two points src and targ onto a sphere centered in 0
    :param src: 3x1 coordinates of the source point onto the sphere
    :param targ: 3x1 coordinates of the target point onto the sphere
    :return:  rot = 3x3 rotation matrix in 3D
    """
    # compute the angle between src and targ
    tet = np.arccos(np.sum((src * targ)) / (np.linalg.norm(src) * np.linalg.norm(targ)))
    # compute the axis arround which the sphere will be rotated = cross product of the vector src and src+targ
    u = np.cross((src + targ), src)
    u = u / np.linalg.norm(u)
    # the rotation is given by the following formula, see e.g. wikipedia rotation matrix ;-)
    c = np.cos(tet)
    s = np.sin(tet)
    rot = np.array([[u[0] ** 2 + (1 - u[0] ** 2) * c, u[0] * u[1] * (1 - c) - u[2] * s, u[0] * u[2] * (1 - c) + u[1] * s], [u[0] * u[1] * (1 - c) + u[2] * s, u[1] ** 2 + (1 - u[1] ** 2) * c, u[1] * u[2] * (1 - c) - u[0] * s], [u[0] * u[2] * (1 - c) - u[1] * s, u[1] * u[2] * (1 - c) + u[0] * s, u[2] ** 2 + (1 - u[2] ** 2) * c]])
    return rot
