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
from soma import aims
from soma.aimsalgo import mesh_coordinates_sphere_resampling

name = 'Make Template Spherical Mesh'
userlevel = 2

signature = Signature(
    'sphere', ReadDiskItem('Ico Mesh', 'GIFTI File'),
    'mesh', ListOf(
        ReadDiskItem('Hemisphere White Mesh', 'Aims mesh formats')),
    'latitude', ListOf(
        ReadDiskItem('Latitude coordinate texture', 'aims Texture formats')),
    'longitude', ListOf(
        ReadDiskItem('Longitude coordinate texture', 'aims Texture formats')),
    'distance', Float(),
    'refined_mesh', WriteDiskItem('Ico Mesh', 'aims Texture formats'))

def initialization(self):
    self.distance = 2.2

def execution(self, context):
    sphere = aims.read(self.sphere.fullPath())
    for index, mesh in enumerate(self.mesh):
        # load object
        lat = aims.read(self.latitude[index].fullPath())
        lon = aims.read(self.longitude[index].fullPath())
        m = aims.read(self.mesh[index].fullPath())
        resampled_mesh = mesh_coordinates_sphere_resampling.resample_mesh_to_sphere(
            m, sphere, lon, lat)
        distance_texture = aims.SurfaceManip.meshEdgeLengthRatioTexture(
            resampled_mesh, sphere)
        if index == 0:
            averaged_texture = distance_texture[0].arraydata()
        else:
            averaged_texture = averaged_texture + distance_texture[0].arraydata()
        

    averaged_texture = averaged_texture / len(self.mesh)
    averaged_texture = aims.TimeTexture(averaged_texture)

    # Builds a sphere mesh with vertices density driven by an average 
    # distance map
    resampled_mesh = mesh_coordinates_sphere_resampling.spere_mesh_from_distance_map(
        sphere, averaged_texture, self.distance)
    
    aims.write(resampled_mesh, self.refined_mesh.fullPath())
    