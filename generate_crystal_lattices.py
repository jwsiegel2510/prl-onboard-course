# Author: Jonathan Siegel
#
# Contains code which generates each of the lattices for assignment 3 of the course.

import pymatgen as mg
import numpy as np

from pymatgen.core.structure import IStructure
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.ext.matproj import MPRester

from user_api_key import USER_API_KEY

def generate_Po_cubic(lattice_param: float):
  lattice_vectors = lattice_param * np.identity(3)
  lattice = Lattice(lattice_vectors)
  return Structure(lattice, ["Po"], [[0.0, 0.0, 0.0]])

def generate_Fe_BCC(lattice_param: float, cell_type: str):
  if cell_type == "primitive":
    lattice_vectors = np.identity(3)
    lattice_vectors[2,:] = [0.5, 0.5, 0.5]
    lattice_vectors *= lattice_param
    lattice = Lattice(lattice_vectors)
    return Structure(lattice, ["Fe"], [[0.0, 0.0, 0.0]])
  elif cell_type == "conventional":
    lattice_vectors = lattice_param * np.identity(3)
    lattice = Lattice(lattice_vectors)
    return Structure(lattice, ["Fe", "Fe"], [[0.0, 0.0, 0.0],
                                             [0.5, 0.5, 0.5]])
  else:
    raise SystemExit("The cell type needs to be either primitive or conventional.")

def generate_Al_FCC(lattice_param: float, cell_type: str):
  if cell_type == "primitive":
    lattice_vectors = np.array([[0.5, 0.5, 0.0],
                                   [0.5, 0.0, 0.5],
                                   [0.0, 0.5, 0.5]])
    lattice_vectors *= lattice_param
    lattice = Lattice(lattice_vectors)
    return Structure(lattice, ["Al"], [[0.0, 0.0, 0.0]])
  elif cell_type == "conventional":
    lattice_vectors = lattice_param * np.identity(3)
    lattice = Lattice(lattice_vectors)
    return Structure(lattice, ["Al", "Al", "Al", "Al"], [[0.0, 0.0, 0.0],
                                                         [0.5, 0.5, 0.0],
                                                         [0.5, 0.0, 0.5],
                                                         [0.0, 0.5, 0.5]])
  else:
    raise SystemExit("The cell type needs to be either primitive or conventional.")

def generate_NaCl(lattice_param: float):
  return Structure(Lattice.cubic(lattice_param), ["Na", "Cl"], [[0.0, 0.0, 0.0],
                                                                [0.5, 0.5, 0.5]])
def generate_Ti_hex(lattice_param_a: float, lattice_param_c: float):
  return Structure(Lattice.hexagonal(lattice_param_a, lattice_param_c), 
                    ["Ti", "Ti", "Ti"], [[0.0, 0.0, 0.0],
                                         [0.66, 0.33, 0.5],
                                         [0.33, 0.66, 0.5]])

# Generate all of the crystal structures and write to output files.
def main():
  Po_cubic = generate_Po_cubic(lattice_param=3.35)
  Po_cubic.to(filename="polonium_cubic.cif")
  Fe_FCC = generate_Fe_BCC(lattice_param=2.25, cell_type="conventional")
  Fe_FCC.to(filename="iron-bcc.cif")
  # A lattice parameter of 4.038 produces a cell volume of 65.89.
  Al_FCC = generate_Al_FCC(lattice_param=4.038, cell_type="conventional")
  Al_FCC.to(filename="aluminum-fcc.cif")
  NaCl = generate_NaCl(lattice_param=5.6402)
  NaCl.to(filename="nacl.cif")
  Ti_hex = generate_Ti_hex(lattice_param_a=4.577, lattice_param_c=2.829)
  Ti_hex.to(filename="ti-hex.cif")

  with MPRester(USER_API_KEY) as m:
    aluminum = m.get_structure_by_material_id("mp-134")
    aluminum.to(filename="al.cif")

if __name__ == "__main__":
  main()
