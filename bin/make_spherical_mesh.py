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
    Make template sphere mesh.
    """

    parser = optparse.OptionParser(desc)

    parser.add_option("-s", "--sphere", dest="sphere")
    parser.add_option("-m", "--mesh",
                      dest="mesh",
                      action="append",
                      help="list of individual hemisphere white mesh")
    parser.add_option("-l", "--latitude",
                      dest="latitude",
                      action="append",
                      help="Latitude coordinate texture")
    parser.add_option("-g", "--longitude",
                      dest="longitude",
                      action="append",
                      help="Longitude coordinate texture")
    parser.add_option("-d", "--distance",
                      dest="distance",
                      help="distance in millimeters")
    parser.add_option("-o", "--output",
                      dest="refined_mesh",
                      help="Refined mesh")

    return parser, parser.parse_args(argv)


def main():
    parser, (options, args) = parseOpts(sys.argv)

    sphere = aims.read(options.sphere)
    for index, mesh in enumerate(options.mesh):
        # load object
        lat = aims.read(options.latitude[index])
        lon = aims.read(options.longitude[index])
        m = aims.read(options.mesh[index])
        resampled_mesh = resampling.resample_mesh_to_sphere(
            m, sphere, lon, lat)
        aims.write(resampled_mesh, '/tmp/resampled_mesh' + str(index) + '.gii')
        distance_texture = aims.SurfaceManip.meshEdgeLengthRatioTexture(
            resampled_mesh, sphere)
        if index == 0:
            averaged_texture = numpy.asarray(
                distance_texture[0], dtype=numpy.float32)
        else:
            averaged_texture = averaged_texture + \
                distance_texture[0].arraydata()

    averaged_texture = averaged_texture / len(options.mesh)
    averaged_texture = aims.TimeTexture(averaged_texture)

    # Builds a sphere mesh with vertices density driven by an average
    # distance map
    aims.write(averaged_texture, '/tmp/averaged_texture.gii')
    resampled_mesh = resampling.spere_mesh_from_distance_map(
        sphere, averaged_texture, float(options.distance))

    aims.write(resampled_mesh, options.refined_mesh)

if __name__ == "__main__":
    main()