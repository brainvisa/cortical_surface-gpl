#!/usr/bin/env python

# python system modules
import sys
import optparse

# soma import
from soma import aims
from soma.aimsalgo import mesh_coordinates_sphere_resampling


def parseOpts(argv):
    desc = """usage: %prog [options] filename"
    Draw a sphere from a triangulation of cortical hemisphere of the subject.
    """

    parser = optparse.OptionParser(desc)

    parser.add_option("-m", "--mesh",
                      dest="mesh",
                      metavar="FILE",
                      help="Individual hemisphere white mesh")
    parser.add_option("-l", "--latitude",
                      dest="latitude",
                      metavar="FILE",
                      help="Latitude coordinate texture")
    parser.add_option("-g", "--longitude",
                      dest="longitude",
                      metavar="FILE",
                      help="Longitude coordinate texture")
    parser.add_option("-o", "--output",
                      dest="sphere_mesh",
                      help="Spherical triangulation on a sphere")

    return parser, parser.parse_args(argv)


def main():
    parser, (options, args) = parseOpts(sys.argv)

    # load objects
    mesh = aims.read(options.mesh)
    latitude = aims.read(options.latitude)
    longitude = aims.read(options.longitude)

    # a spherical triangulation of the subject of its cortical hemisphere,
    # projected on a sphere
    sphere_mesh = mesh_coordinates_sphere_resampling.draw_sphere(
        mesh, longitude, latitude)

    # write the new mesh
    aims.write(sphere_mesh, options.sphere_mesh)

if __name__ == "__main__":
    main()