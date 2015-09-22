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
from brainvisa.processes import *
from soma import aims
from soma.aimsalgo import mesh_coordinates_sphere_resampling

name = "Resampled Mesh"
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
)


def initialization(self):
    self.linkParameters("latitude", "mesh")
    self.linkParameters("longitude", "mesh")
    self.linkParameters("resampled_mesh", "mesh")

def execution(self, context):
    sphere = aims.read(self.sphere.fullPath())
    mesh = aims.read(self.mesh.fullPath())
    lon = aims.read(self.longitude.fullPath())
    lat = aims.read(self.latitude.fullPath())
    resampled_mesh = mesh_coordinates_sphere_resampling.resample_mesh_to_sphere(
        mesh,sphere,lon,lat)
    aims.write(resampled_mesh, self.resampled_mesh.fullPath())
        