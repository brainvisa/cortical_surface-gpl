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

name = "Resample Texture From Coordinates Fields"
userLevel = 2


# Argument declaration
signature = Signature(
    "texture", ReadDiskItem("Texture", "aims texture formats"),
    "sphere", ReadDiskItem("Ico Mesh", "aims Texture formats"),
    "mesh", ReadDiskItem("Hemisphere White Mesh", "Aims mesh formats"),
    "latitude", ReadDiskItem(
        "Latitude coordinate texture", "aims Texture formats"),
    "longitude", ReadDiskItem(
        "Longitude coordinate texture", "aims Texture formats"),
    "interpolation", Choice("linear", "nearest_neighbour"),
    "resampled_texture", WriteDiskItem(
        "Texture", "Aims texture formats",
        requiredAttributes={"vertex_corr": "Yes",
                            "vertex_corr_method": "hiphop"}),
)


def initialization(self):
    self.linkParameters("mesh", "texture")
    self.linkParameters("latitude", "mesh")
    self.linkParameters("longitude", "mesh")
    self.linkParameters("resampled_texture", "texture")
    self.interpolation = "nearest_neighbour"
    self.sphere = self.signature['sphere'].findValue({})

def execution(self, context):
    sphere = aims.read(self.sphere.fullPath())
    mesh = aims.read(self.mesh.fullPath())
    lon = aims.read(self.longitude.fullPath())
    lat = aims.read(self.latitude.fullPath())
    texture = aims.read(self.texture.fullPath())
    resampled_texture \
        = mesh_coordinates_sphere_resampling.resample_texture_to_sphere(
            mesh, sphere, lon, lat, texture, self.interpolation)
    aims.write(resampled_texture, self.resampled_texture.fullPath())
