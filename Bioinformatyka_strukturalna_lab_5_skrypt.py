# -*- coding: utf-8 -*-
from Bio import PDB
import numpy as np
from Bio.PDB import Superimposer
import os

def calculate_rmsd(coords1, coords2):
    diff = coords1 - coords2
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    return rmsd

# Wczytaj strukturê referencyjn¹
reference_structure = PDB.PDBParser(QUIET=True).get_structure("reference", "R1107_reference.pdb")

# Wczytaj modele
models = []
model_filenames = ["struktura1.pdb", "struktura2.pdb", "struktura3.pdb", "struktura4.pdb"]
for filename in model_filenames:
    model = PDB.PDBParser(QUIET=True).get_structure(filename[:-4], filename)
    models.append(model)

# Utwórz obiekt Superimposer
superimposer = Superimposer()

# Iteruj przez modele
for model, filename in zip(models, model_filenames):
    # Pobierz atomy do superpozycji
    atoms_ref = [atom for atom in reference_structure.get_atoms()]
    atoms_model = [atom for atom in model.get_atoms()]

    # PrzeprowadŸ superpozycjê
    superimposer = Superimposer()
    superimposer.set_atoms(atoms_ref, atoms_model)
    superimposer.apply(model)

    # Pobierz wspó³rzêdne atomów do obliczenia RMSD
    atom_coords_ref = np.array([atom.get_coord() for atom in atoms_ref])
    atom_coords_model = np.array([atom.get_coord() for atom in atoms_model])

    rmsd_value = calculate_rmsd(atom_coords_ref, atom_coords_model)

    print(f"RMSD dla {filename} wynosi: {rmsd_value}")

    # Skrypt PyMOLa
    pymol_script = f"""
    load {filename}
    super {filename}, reference
    color red, reference
    color blue, {filename}
    show cartoon
    zoom
    """

    with open(f"pymol_script_{filename}.pml", "w") as f:
        f.write(pymol_script)

    # Uruchom PyMOLa z skryptem
    pymol_path = "D:\\PyMOL\\PyMOLWin.exe"
    os.system(f"{pymol_path} -c pymol_script_{filename}.pml")