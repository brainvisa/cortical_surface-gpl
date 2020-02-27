#! /usr/bin/env python
############################################################################
# This software and supporting documentation are distributed by
# CEA/NeuroSpin, Batiment 145, 91191 Gif-sur-Yvette cedex, France.
# This software is governed by the CeCILL license version 2 under
# French law and abiding by the rules of distribution of free software.
# You can  use, modify and/or redistribute the software under the
# terms of the CeCILL license version 2 as circulated by CEA, CNRS
# and INRIA at the following URL "http://www.cecill.info".
############################################################################

from __future__ import print_function

# python system modules
from __future__ import absolute_import
import sys
import numpy
import optparse

# soma import
from soma import aims
import soma.aimsalgo.mesh_coordinates_sphere_resampling as resampling


def parseOpts(argv):
    desc = """usage: %prog [options] filename"
    Make template spherical mesh.
    
    Builds a sphere mesh with vertices density driven by an average distance
    map, coming with another initial sphere mesh and a target length.
    """

    parser = optparse.OptionParser(desc)

    parser.add_option("-s", "--icosphere",
                      dest="ico_sphere",
                      help="An initial sphere mesh, typically an icosphere")
    parser.add_option("-m", "--mesh",
                      dest="mesh",
                      action="append",
                      help="List of individual hemisphere white mesh")
    parser.add_option("-l", "--latitude",
                      dest="latitude",
                      action="append",
                      help="List of individual latitude coordinate texture")
    parser.add_option("-g", "--longitude",
                      dest="longitude",
                      action="append",
                      help="List of individual longitude coordinate texture")
    parser.add_option("-d", "--distance",
                      dest="distance",
                      help="Distance in millimeters")
    parser.add_option("-o", "--output",
                      dest="refined_mesh",
                      help="Refined mesh")
    parser.add_option("-t", "--trueinversion",
                      dest="inversion", default=False, action="store_true",
                      help="Check in the case of righ hemisphere")
    parser.add_option("-f", "--falseinversion",
                      dest="inversion", default=False, action="store_false",
                      help="Check in the case of righ hemisphere")

    return parser, parser.parse_args(argv)


def main():
    parser, (options, args) = parseOpts(sys.argv)

    # load sphere
    ico_sphere = aims.read(options.ico_sphere)
    for index, mesh in enumerate(options.mesh):
        # load object (subject by subject)
        lat = aims.read(options.latitude[index])
        lon = aims.read(options.longitude[index])
        m = aims.read(options.mesh[index])

        # resample a mesh to the sphere.
        print("inversion dans la commande", options.inversion)
        resampled_mesh = resampling.resample_mesh_to_sphere(
            m, ico_sphere, lon, lat, inversion=options.inversion)
        distance_texture = aims.SurfaceManip.meshEdgeLengthRatioTexture(
            resampled_mesh, ico_sphere)
        # add all textures
        if index == 0:
            print()
            averaged_texture = numpy.asarray(
                distance_texture[0], dtype=numpy.float32)
        else:
            averaged_texture = averaged_texture + \
                numpy.asarray(distance_texture[0], dtype=numpy.float32)

    # average textures
    averaged_texture = averaged_texture / len(options.mesh)
    averaged_texture = aims.TimeTexture(averaged_texture)

    # Builds a sphere mesh with vertices density driven by an average
    # distance map
    resampled_mesh = resampling.sphere_mesh_from_distance_map(
        ico_sphere, averaged_texture, float(options.distance), 
        inversion=options.inversion)

    # write the final mesh (resampled)
    aims.write(resampled_mesh, options.refined_mesh)

if __name__ == "__main__":
    main()