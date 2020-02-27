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
from __future__ import absolute_import
from brainvisa.processes import *

# soma import
from soma.path import find_in_path

name = "Draw a sphere"
userlevel = 2

signature = Signature(
    "white_mesh", ReadDiskItem("Hemisphere White Mesh", "Aims mesh formats"),
    "latitude", ReadDiskItem(
        "Latitude coordinate texture", "aims Texture formats"),
    "longitude", ReadDiskItem(
        "Longitude coordinate texture", "aims Texture formats"),
    "sphere_mesh", WriteDiskItem(
        "Spherical Mesh", "aims Texture formats"),
    "inversion", Boolean())


def initialization(self):
    self.inversion = False
    self.linkParameters("latitude", "white_mesh")
    self.linkParameters("longitude", "latitude")


def execution(self, context):
    cmd_args = []
    
    if self.inversion:
        cmd_args += ["-t", self.inversion]
    else:
        cmd_args += ["-f", self.inversion]  
        
    cmd_args += ["-m", self.white_mesh, "-l", self.latitude, 
                 "-g", self.longitude, "-o", self.sphere_mesh]
                 
    context.system(
        sys.executable, find_in_path("mesh_to_sphere.py"), *cmd_args)
  