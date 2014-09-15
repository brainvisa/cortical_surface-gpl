#!/usr/bin/env python

# python system modules
import sys
import optparse
import numpy

# soma import
from soma import aims
import soma.aimsalgo.mesh_coordinates_sphere_resampling as resampling


def parseOpts(argv):
    desc = """usage: %prog [options] filename"
    Make template spherical mesh.
    """

    parser = optparse.OptionParser(desc)

    parser.add_option("-s", "--sphere",
                      dest="sphere",
                      help="A sphere mesh with center 0")
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

    return parser, parser.parse_args(argv)


def main():
    parser, (options, args) = parseOpts(sys.argv)

    # load sphere
    sphere = aims.read(options.sphere)
    for index, mesh in enumerate(options.mesh):
        # load object (subject by subject)
        lat = aims.read(options.latitude[index])
        lon = aims.read(options.longitude[index])
        m = aims.read(options.mesh[index])

        # resample a mesh to the sphere.
        resampled_mesh = resampling.resample_mesh_to_sphere(
            m, sphere, lon, lat)
        distance_texture = aims.SurfaceManip.meshEdgeLengthRatioTexture(
            resampled_mesh, sphere)

        # add all textures
        if index == 0:
            averaged_texture = numpy.asarray(
                distance_texture[0], dtype=numpy.float32)
        else:
            averaged_texture = averaged_texture + \
                distance_texture[0].arraydata()

    # average textures
    averaged_texture = averaged_texture / len(options.mesh)
    averaged_texture = aims.TimeTexture(averaged_texture)

    # Builds a sphere mesh with vertices density driven by an average
    # distance map
    resampled_mesh = resampling.spere_mesh_from_distance_map(
        sphere, averaged_texture, float(options.distance))

    # write the final mesh (resampled)
    aims.write(resampled_mesh, options.refined_mesh)

if __name__ == "__main__":
    main()