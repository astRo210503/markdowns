import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import RWMol, SanitizeMol, rdchem
import pandas as pd
import selfies as sf

# Function to check valence
def is_valence_ok(atom_element, existing_bonds, new_bond_order=1):
    max_valences = {'H': 1, 'O': 2, 'C': 4, 'Cl': 1}
    return (existing_bonds + new_bond_order) <= max_valences.get(atom_element, 4)

# Function to add atom and bond
def add_atom_and_bond(molecule, atom_index, element):
    atom = molecule.GetAtomWithIdx(atom_index)
    existing_bonds = sum([bond.GetBondTypeAsDouble() for bond in atom.GetBonds()])
    if not is_valence_ok(atom.GetSymbol(), existing_bonds):
        return None

    new_atom_index = molecule.AddAtom(Chem.Atom(element))
    molecule.AddBond(atom_index, new_atom_index, rdchem.BondType.SINGLE)
    try:
        Chem.SanitizeMol(molecule)
        return new_atom_index
    except:
        molecule.RemoveAtom(new_atom_index)
        return None

# Function to generate molecules
def generate_molecules(max_depth, current_molecule, current_depth=0, molecules=None):
    if molecules is None:
        molecules = []

    smiles = Chem.MolToSmiles(current_molecule)
    molecules.append(smiles)

    if current_depth <= max_depth:
        for atom_index in range(current_molecule.GetNumAtoms()):
            for element in ['C', 'O', 'H', 'Cl']:
                existing_bonds = sum([bond.GetBondTypeAsDouble() for bond in current_molecule.GetAtomWithIdx(atom_index).GetBonds()])
                if is_valence_ok(element, existing_bonds):
                    new_molecule = Chem.RWMol(current_molecule)
                    if add_atom_and_bond(new_molecule, atom_index, element) is not None:
                        generate_molecules(max_depth, new_molecule, current_depth + 1, molecules)

    return list(set(molecules))

# Updated generate_valid_smiles function
def generate_valid_smiles(elements, max_depth):
    init_molecule = Chem.RWMol()
    init_molecule.AddAtom(Chem.Atom(elements[0]))
    generated_molecules = generate_molecules(max_depth, init_molecule)
    return generated_molecules

# Streamlit app
st.title("SMILES Generator and Visualizer")

# Define the tabs
tab1, tab2 = st.tabs(["Generate SMILES", "Display SMILES"])

# Tab 1: Generate SMILES
with tab1:
    st.header("Generate SMILES")
    elements = st.text_input("Enter elements (comma separated, e.g., C,O,N)")
    max_depth = st.number_input("Enter max depth", min_value=1, value=2)
    if st.button("Generate"):
        elements_list = elements.split(",")
        smiles_list = generate_valid_smiles(elements_list, max_depth)
        st.write("Generated SMILES strings:")
        smiles_df = pd.DataFrame(smiles_list, columns=["SMILES"])
        st.dataframe(smiles_df)

# Tab 2: Display SMILES
with tab2:
    st.header("Display SMILES")
    smiles_string = st.text_input("Enter SMILES string")
    if st.button("Display"):
        if smiles_string:
            mol = Chem.MolFromSmiles(smiles_string)
            if mol:
                img = Draw.MolToImage(mol)
                st.image(img, caption=f"SMILES: {smiles_string}")
            else:
                st.error("Invalid SMILES string")
        else:
            st.error("Please enter a SMILES string")
