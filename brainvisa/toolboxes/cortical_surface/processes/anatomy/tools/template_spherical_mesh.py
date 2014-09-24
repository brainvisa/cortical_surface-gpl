############################################################################
# This software and supporting documentation are distributed by
# CEA/NeuroSpin, Batiment 145, 91191 Gif-sur-Yvette cedex, France.
# This software is governed by the CeCILL license version 2 under
# French law and abiding by the rules of distribution of free software.
# You can  use, modify and/or redistribute the software under the
# terms of the CeCILL license version 2 as circulated by CEA, CNRS
# and INRIA at the following URL "http://www.cecill.info".
############################################################################

# brainvisa import
from brainvisa.processes import *

# soma import
from soma.path import find_in_path

name = "Make Template Spherical Mesh"
userlevel = 2

signature = Signature(
    "sphere", ReadDiskItem("Ico Mesh", "Aims mesh formats"),
    "mesh", ListOf(
        ReadDiskItem("Hemisphere White Mesh", "Aims mesh formats")),
    "latitude", ListOf(
        ReadDiskItem("Latitude coordinate texture", "Aims texture formats")),
    "longitude", ListOf(
        ReadDiskItem("Longitude coordinate texture", "Aims texture formats")),
    "distance", Float(),
    "refined_mesh", WriteDiskItem("Spherical Mesh", "Aims mesh formats"),
    "inversion", Boolean())


def initialization(self):
    self.distance = 4.0
    self.inversion = False
    self.linkParameters("latitude", "mesh")
    self.linkParameters("longitude", "latitude")


def execution(self, context):
    cmd_args = []
    context.write("inversion: ", self.inversion)
    if self.inversion:
        cmd_args += ["-t", self.inversion]
    else:
        cmd_args += ["-f", self.inversion]
    for m in self.mesh:
        cmd_args += ["-m", m]
    for lat in self.latitude:
        cmd_args += ["-l", lat]
    for lon in self.longitude:
        cmd_args += ["-g", lon]
    cmd_args += ["-s", self.sphere, "-d", str(self.distance),
                 "-o", self.refined_mesh]
    context.system(
        "python", find_in_path("make_spherical_mesh.py"), *cmd_args)