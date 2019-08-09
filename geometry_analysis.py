# -*- coding: utf-8 -*-
"""
geometry.analysis.py
This module contains the geometry analysis project.
"""

import numpy
import os
import sys
# sys library allows us to use command line arguments 

def open_xyz(filename): 
    """This function reads in and processes an xyz file.
    Input: filename
    Returns: symbols and coordinates
    """
    xyz_file = numpy.genfromtxt(fname=filename, skip_header=2, dtype='unicode')
    symbols = xyz_file[:,0]
    coordinates = xyz_file[:,1:]
    coordinates = coordinates.astype(numpy.float)
    return symbols, coordinates  


def calculate_distance(atom1_coord, atom2_coord):
    """This function takesdi the coordinates of two atoms and calculates the distances between them.
    Input: atom1_coord, atom2_coord
    Returns: distance
    """
    x_distance = atom1_coord[0] - atom2_coord[0]
    y_distance = atom1_coord[1] - atom2_coord[1]
    z_distance = atom1_coord[2] - atom2_coord[2]
    distance = numpy.sqrt(x_distance**2 + y_distance**2 + z_distance**2)
    return distance

def bond_check(distance, minimum_length=0, maximum_length=1.5):
    """
    This function checks if the distance is between the specified minimum (default: 0 Angstroms) and maximum values 
    (default: 1.5).
    Input: distance 
    Returns: True or False
    
    """
    
    if distance < 0:
        raise ValueError("A negative distance has been detected. Please check your input! ")
    
    if distance > minimum_length and distance < maximum_length:
        return True
    else:
        return False

#coordinate_file = os.path.join('data', 'water.xyz')
if __name__ == "__main__":
    print(F"Running {sys.argv[0]}")

    if len(sys.argv) < 2:
        raise NameError("Incorrect input! Please specifiy an xyz file to be analyzed!")


    file_location=sys.argv[1]
    
    atom, data = open_xyz(file_location)


    data = data.astype(numpy.float)


    num_atoms = len(atom)


    for atom1 in range(0,num_atoms):
        for atom2 in range(0, num_atoms):
            if atom1 < atom2:
                bondlength_12 = calculate_distance(data[atom1], data[atom2])
                if bond_check(bondlength_12) is True:
                    print(F'{atom[atom1]} to {atom[atom2]} : {bondlength_12 : .3F}')

