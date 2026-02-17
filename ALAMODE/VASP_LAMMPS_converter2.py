#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  VASP_LAMMPS_converter2.py
#
#  Copyright 2026 Diego Martinez Gutierrez <diego.martinez@ehu.eus>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#
"""
Este script python utiliza ASE para convertir los datos del formato usado por
LAMMPS al formato usado por VASP.
"""
from ase.io import read, write
from ase import Atoms
import glob

def fromVASPtoLAMMPS(arg):
    atoms = read(arg, format='lammps-data')
    write(f'{arg}.lammps', atoms, format='lammps-data', style='atomic')
    pass

def fname(arg):
    # Example: Assign atom types
    atoms = Atoms('Si2O', positions=[[0, 0, 0], [1, 1, 1], [2, 2, 2]])
    atoms.set_tags([1, 1, 2])  # Assign atom types (e.g., 1 for Si, 2 for O)
    atoms.set_initial_charges([0.0, 0.0, -1.0])  # Set charges for each atom

    write(f'{poscar}.lammps', atoms, format='vasp')
    pass

if __name__ == '__main__':
#    for poscar in glob.glob('POSCAR'):
    atoms = read('relax.dat', format='lammps-data')
    write(f'AA222.SPOSCAR', atoms, format='vasp')
    print(atoms.get_chemical_symbols())  # Check atom types
    print(atoms.get_positions())         # Check positions
    print(atoms.get_initial_charges())   # Check charges, if applicable

def image(atoms):
    from ase.visualize import view
    write('image.png', atoms)
