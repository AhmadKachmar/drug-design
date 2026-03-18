" drug design project"
import os
import networkx as nx
import argparse
import rdkit
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import rdmolops

def edge_density(graph):
    """Computes the edge density of graph """
    num_nodes = graph.number_of_nodes()
    num_edges = graph.number_of_edges()
    
    if num_nodes > 1:
        return (2 * num_edges) / (num_nodes * (num_nodes - 1)) 
    else: 
        return 0
    
    
def wiener_index(graph):
    """Computes the Wiener Index (sum of shortest path distances)."""
    wiener_index = 0
    for node1 in graph.nodes:
        for node2 in graph.nodes:
            if node1 < node2:
                wiener_index += nx.shortest_path_length(graph, node1, node2)
    return wiener_index


def petitjean_index(graph):
    """Computes the PetitJean Index."""
    shortest_paths = nx.all_pairs_shortest_path_length(graph)
    eccentricities = {}
    for node, paths in shortest_paths:
        eccentricities[node] = max(paths.values())
    D = max(eccentricities.values())
    R = min(eccentricities.values())
    if D != 0:
        return (D - R) / R
    return 0


def mol_to_graph(mol):
    """Converts an RDKit molecule to a NetworkX graph."""
    graph = nx.Graph()
    for atom in mol.GetAtoms():
        graph.add_node(atom.GetIdx())
    for bond in mol.GetBonds():
        graph.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    return graph


def visualize_graph(graph, title):
    """Displays the molecular graph."""
    plt.figure(figsize=(6, 6))
    pos = nx.spring_layout(graph)
    nx.draw(graph, pos, with_labels=True, node_color='lightblue', edge_color='gray', node_size=500, font_size=10)
    plt.title(title)
    plt.show()

def process_file(file_path):
    """Reads a .mol or .sdf file and returns a list of molecule names and structures."""
    molecules = []
    if file_path.endswith(".mol"):
        mol = Chem.MolFromMolFile(file_path)
        if mol is not None:
            molecules.append((os.path.basename(file_path), mol))
    elif file_path.endswith(".sdf"):
        suppl = Chem.SDMolSupplier(file_path)
        for i, mol in enumerate(suppl, start=1):
            if mol is not None:
                molecules.append((f"Molecule {i}", mol))
    return molecules

def compute_indices(molecules, indices):
    """Calculates the chosen indices for molecules."""
    results = []
    for name, mol in molecules:
        graph = mol_to_graph(mol)
        visualize_graph(graph, name)
        print("Molecule:", name)
        if "edge_density" in indices:
            print("Edge Density:", edge_density(graph))
        if "wiener_index" in indices:
            print("Wiener Index:", wiener_index(graph))
        if "petitjean_index" in indices:
            print("PetitJean Index:", petitjean_index(graph))
    return results

    
def main():
    """Handles user input and computes indices for molecules."""
    input_path = input("Enter the file or folder path: ")
    indices_choice = input("Enter indices to compute (edge_density, wiener_index, petitjean_index, all): ")
    
    if indices_choice == "all":
        indices = ["edge_density", "wiener_index", "petitjean_index"]
    else:
        indices = indices_choice.split()
    
    molecules = []
    if os.path.isdir(input_path):
        for filename in os.listdir(input_path):
            if filename.endswith(".mol"):
                molecules.extend(process_file(os.path.join(input_path, filename)))
    else:
        molecules = process_file(input_path)
    
    compute_indices(molecules, indices)

if __name__ == "__main__":
    main()