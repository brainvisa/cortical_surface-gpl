############################################################################
#  This software and supporting documentation are distributed by
#      CEA/NeuroSpin, Batiment 145,
#      91191 Gif-sur-Yvette cedex
#      France
# This software is governed by the CeCILL license version 2 under
# French law and abiding by the rules of distribution of free software.
# You can  use, modify and/or redistribute the software under the
# terms of the CeCILL license version 2 as circulated by CEA, CNRS
# and INRIA at the following URL "http://www.cecill.info".
############################################################################

# BrainVisa modules
from __future__ import absolute_import
from brainvisa.processes import *
from soma import aims
from soma.aimsalgo import mesh_coordinates_sphere_resampling
import numpy as np
from brainvisa import registration

name = "Resample Mesh From Coordinates Fields"
userLevel = 2


# Argument declaration
signature = Signature(
    "sphere", ReadDiskItem("Ico Mesh", "aims Texture formats"),
    "mesh", ReadDiskItem("Hemisphere White Mesh", "Aims mesh formats"),
    "latitude", ReadDiskItem(
        "Latitude coordinate texture", "aims Texture formats"),
    "longitude", ReadDiskItem(
        "Longitude coordinate texture", "aims Texture formats"),
    "resampled_mesh", WriteDiskItem(
        "Remeshed mesh", "Aims mesh formats"),
    "invert_sphere", Boolean(),
)


def linkInversion(self, proc, dummy):
    if self.mesh is not None:
        side = self.mesh.get("side")
        if side == "right":
            return True
    return False


def initialization(self):
    self.invert_sphere = False
    self.linkParameters("latitude", "mesh")
    self.linkParameters("longitude", "mesh")
    self.linkParameters("resampled_mesh", "mesh")
    self.linkParameters("invert_sphere", "mesh", self.linkInversion)
    self.sphere = self.signature['sphere'].findValue({})


def execution(self, context):
    sphere = aims.read(self.sphere.fullPath())
    mesh = aims.read(self.mesh.fullPath())
    lon = aims.read(self.longitude.fullPath())
    lat = aims.read(self.latitude.fullPath())
    if self.invert_sphere:
        sphere.vertex().assign(
            [aims.Point3df(x)
             for x in np.asarray(sphere.vertex()) * [-1, 1, 1]])
    resampled_mesh \
        = mesh_coordinates_sphere_resampling.resample_mesh_to_sphere(
            mesh, sphere, lon, lat) #, inversion=self.invert_longitude)
    aims.write(resampled_mesh, self.resampled_mesh.fullPath())
    tm = registration.getTransformationManager()
    tm.copyReferential(self.mesh, self.resampled_mesh)

