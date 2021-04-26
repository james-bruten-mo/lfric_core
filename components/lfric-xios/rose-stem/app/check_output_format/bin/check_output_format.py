#!/usr/bin/env python
##############################################################################
# (c) Crown copyright 2021 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
'''
Check that diagnostic output from apps matches what is expected for UGRID
output with XIOS.
'''
import sys
import glob
import iris
import numpy as np


def load_cube_by_varname(filename, var):
    '''
    Loads an Iris cube according to the varname
        filename:   a string holding the path to the file to be tested
        var:        the varname id to be loaded
    '''
    variable_constraint = iris.Constraint(cube_func=(lambda c:
                                          c.var_name == var))
    return iris.load_cube(filename, constraint=variable_constraint)


def identify_mesh(filename):
    '''
    Identifies the mesh type and size based on the header of the diagnostic
    file
        filename:   a string holding the path to the file to be tested
    '''
    mesh_dict = {}
    links_dict = {}

    mesh_face_links = load_cube_by_varname(filename, "Mesh2d_face_face_links")
    n_mesh_faces = mesh_face_links.shape[0]

    if np.sqrt((n_mesh_faces/6)).is_integer():
        mesh_dict['type'] = "cubedsphere"
        mesh_dict['panel_length'] = int(np.sqrt((n_mesh_faces/6)))
        mesh_dict['name'] = "C{0}".format(mesh_dict['panel_length'])

        links_dict['face_face'] = n_mesh_faces
        links_dict['face_edge'] = n_mesh_faces * 2
        links_dict['edge_edge'] = n_mesh_faces * 2
        links_dict['node'] = n_mesh_faces + 2

        mesh_dict['links'] = links_dict

    else:
        sys.exit("Unable to identify mesh in {0}".format(filename))

    return mesh_dict


def test_header(filename, mesh_dict):
    '''
    Tests the mesh connectivity information from the NetCDF header against the
    expected values for the identified mesh
        filename:   a string holding the path to the file to be tested
        mesh_dict:  a dictionary holding information about the mesh
    '''
    cubes = iris.load(filename)
    header_errors = []

    mesh_links = mesh_dict['links']

    for link in mesh_links:
        for cube in cubes:
            if "Mesh2d" in cube.var_name:
                if "_{0}_".format(link) in cube.var_name:
                    if cube.shape[0] != mesh_links[link]:
                        header_errors.append([cube.var_name, cube.shape[0],
                                              mesh_links[link]])

    return header_errors


if __name__ == "__main__":

    try:
        output_path = sys.argv[1]
    except ValueError:
        sys.exit("Usage: {0} <output_path>".format(sys.argv[0]))

    TOTAL_ERRORS = 0

    diagnostic_paths = glob.glob(output_path + "*/*/*/*/*/*diag.nc")
    for data_path in diagnostic_paths:

        mesh = identify_mesh(data_path)
        errors = test_header(data_path, mesh)

        TOTAL_ERRORS = TOTAL_ERRORS + len(errors)

        if errors:
            print("\n{0} errors found in {1}:".format(len(errors), data_path))
            print("{0} {1} mesh identified...".format(mesh['name'],
                                                      mesh['type']))

            for err in errors:
                print("{0}: expected {1}, got {2}".format(err[0], err[2],
                                                          err[1]))

    if TOTAL_ERRORS > 0:
        sys.exit("Errors found in output file mesh links - see std.out")
