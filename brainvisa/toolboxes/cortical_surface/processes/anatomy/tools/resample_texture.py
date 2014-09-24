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

name = "Resample Texture"
userlevel = 2

signature = Signature(
    "white_mesh", ReadDiskItem("AimsWhite", "Aims mesh formats"),
    "gyri_segmentation", ReadDiskItem("ResampledGyri", "Aims texture formats"),
    "resampled_white_mesh", ReadDiskItem(
        "Resampled Hemisphere White Mesh", "Aims mesh formats"),
    "resampled_gyri_segmentation", WriteDiskItem(
        "Resampled Hemisphere Gyri Texture", "Aims texture formats"),)


def initialization(self):
    self.linkParameters("gyri_segmentation", "white_mesh")
    self.linkParameters("white_mesh", "resampled_white_mesh")
    self.linkParameters("gyri_segmentation", "resampled_gyri_segmentation")


def execution(self, context):
    white_mesh = aims.read(self.white_mesh.fullPath())
    gyri_segmentation = aims.read(self.gyri_segmentation.fullPath())
    resampled_white_mesh = aims.read(self.resampled_white_mesh.fullPath())

    interpoler = aims.MeshInterpoler(white_mesh, resampled_white_mesh)
    interpoler.project()

    texture = interpoler.resampleTexture(
        gyri_segmentation, aims.MeshInterpoler.NearestNeighbour)

    aims.write(texture, self.resampled_gyri_segmentation.fullPath())