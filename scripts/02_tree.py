# Install necessary packages if not already installed
# !pip install ete4 pandas matplotlib seaborn

import pandas as pd
from ete4 import Tree
from ete4.treeview import TreeStyle, NodeStyle, AttrFace, TextFace, ImgFace, faces
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
import seaborn as sns
import re

# Function to preprocess Newick string
def preprocess_newick(newick_str):
    # Replace non-numeric bootstrap values with numeric ones
    newick_str = re.sub(r'(\d+)/(\d+)', r'\1', newick_str)
    return newick_str

# Load the Newick tree file and preprocess it
with open("/playground/dataset_salmonella/bactopia-runs/snippy-20240321-155218/trimal_out_seq90_res07.aln.treefile") as f:
    newick_str = f.read()

newick_str = preprocess_newick(newick_str)

# Load the phylogenetic tree
tree = Tree(newick_str)

tree.write(format=1, outfile="/playground/dataset_salmonella/bactopia-runs/snippy-20240321-155218/trimal_out_seq90_res07.aln.treefile.nw")
# Load the CSV file
df = pd.read_csv("/playground/dataset_salmonella/bactopia-runs/snippy-20240321-155218/metadata.csv")

# Set silhouette sources
img_path = "/home/jlvillanueva/Documents/SALMONELLA/scripts/"
humanFace = faces.ImgFace(img_path+"human.svg")
pidgeonFace = faces.ImgFace(img_path+"columba_livia.png")
larus1Face = faces.ImgFace(img_path+"larus1.png")
larus2Face = faces.ImgFace(img_path+"larus2.png")

# Ensure the column names are appropriate for your dataset
leaf_name_column = "accession"  # Replace with the column name matching the leaf names in your tree
isolation_source_column = "isolation_source_fixed"

# Create a dictionary mapping leaf names to isolation sources
isolation_source_dict = pd.Series(df[isolation_source_column].values, index=df[leaf_name_column]).to_dict()

# Generate a color palette based on 'set1' from ggplot2
palette = sns.color_palette("deep", len(df[isolation_source_column].unique()))
color_dict = dict(zip(df[isolation_source_column].unique(), palette))

# Function to get color for a node
def get_node_color(node_name):
    return to_hex(color_dict.get(isolation_source_dict.get(node_name, "white")))

# Function to apply custom style to tree nodes
def apply_node_styles(tree):
    for node in tree.traverse():
        if node.is_leaf:
            sample_name = node.name
            if sample_name in isolation_source_dict:
                node.img_style["bgcolor"] = get_node_color(sample_name)
                node.name=sample_name + ': ' + isolation_source_dict[sample_name]
            else:
                node.img_style["bgcolor"] = "white"  # Default color for other leaf nodes
            node.img_style["size"] = 0
            if sample_name in human_samples:
                node.add_face(humanFace, column=1, position='aligned')
            elif sample_name in audouin_samples:
                node.add_face(larus1Face, column=1, position='aligned')
        else:
                node.img_style["bgcolor"] = "white"  # Default color for other leaf nodes

# Define subset of samples
samples_of_interest = [
    "SAL1-14_12598", "SAL7-16_18792", "SAL2-13_13995", "REF-GSJ_2017-Sal-008",
    "SAL3-14_11275", "SAL5-14_28823", "SAL4-14_28806", "SAL6-16_18784", "Reference"
]

human_samples = ["SAL5-14_28823", "SAL6-16_18784","SAL4-14_28806", "SAL1-14_12598", "SAL7-16_18792"
]

audouin_samples = ["SAL3-14_11275", "SAL2-13_13995"
]

# Apply styles to the full tree
apply_node_styles(tree)

michahellis_samples = []

'''
pruned_tree = tree.copy()
for leaf in pruned_tree.get_leaves():
    if leaf.name not in samples_of_interest:
        leaf.detach()

# Apply styles to the pruned tree
apply_node_styles(pruned_tree)
'''
# Add a legend to the plot
legend = [TextFace(source, fsize=21, fgcolor=to_hex(color)) for source, color in color_dict.items()]

# Define tree style
def get_tree_style():
    ts = TreeStyle()
    ts.mode = "c"
    ts.show_leaf_name = True
    ts.legend.add_face(TextFace("Isolation Source:", fsize=24, bold=True), column=0)
    for text_face in legend:
        ts.legend.add_face(text_face, column=0)
    ts.legend_position = 3
    return ts

# Define tree style
def get_tree_style2():
    ts = TreeStyle()
    ts.mode = "r"
    ts.show_leaf_name = True
    #ts.legend.add_face(TextFace("Isolation Source:", fsize=24, bold=True), column=0)
    for text_face in legend:
        ts.legend.add_face(text_face, column=0)
    ts.legend_position = 3
    return ts

# Render the full tree
tree.render("full_circular_tree.png", w=1832, units="px", tree_style=get_tree_style())

# Load the Newick tree file and preprocess it
with open(img_path+ "subtree.tree") as f:
    newick_str = f.read()

newick_str = preprocess_newick(newick_str)

# Load the phylogenetic tree
subtree = Tree(newick_str)

# Apply styles to the full tree
apply_node_styles(subtree)


# Render the pruned tree
pruned_tree.render("pruned_circular_tree.png", w=1832, units="px", tree_style=get_tree_style())

# Optionally display the trees in interactive windows (useful during development)
tree.explore(name="thetree")
#tree.show(tree_style=get_tree_style())
#subtree.show(tree_style=apply_node_styles())