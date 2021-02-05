import os.path
import xml.etree.ElementTree as ET
import sys
import networkx as nx
import numpy as np
"""
	Taken from graphfiles.py
	Reads an xml or cxl file listing all the graph files in a collection.
	Uses the function loadGXL to load each graph file
	parameter:
	- filename: path and file name of the xml or cxl file containing the listing of the graphs in the collection.
	- dataset_path: if the graph files are not in the same folder that the collection file, this variable specifies the directory where the graph files are to be read from. If not included, the files are searched for in the same folder as the collection file.
	
	returns:
	- data: list of nx graph objects
	- y: list corresponding to the class of each graph
	- struct: the node and edge structures. See loadGXL for description
"""
def loadFromXML(filename, dataset_path=None):
	if not dataset_path is None:
		dirname_dataset = dataset_path
	else:
		dirname_dataset = os.path.dirname(filename)
	tree = ET.parse(filename)
	root = tree.getroot()
	data = []
	y = []
	struct = {}
	struct_aux={}
	node_labels = {}
	edge_labels = {}
	
	for graph in root.iter('graph'):
		mol_filename = graph.attrib['file']
		mol_class = graph.attrib['class']
		g, struct_aux = loadGXL(dirname_dataset + '/' + mol_filename)
		g.name = mol_filename
		for a1 in struct_aux['nodes']:
			if not a1 in node_labels: node_labels[a1] = []
			attrs = nx.get_node_attributes(g,a1)			
			node_labels[a1] += list(attrs.values())

		for a2 in struct_aux['edges']:
			if not a2 in edge_labels: edge_labels[a2] = []
			attrs = nx.get_edge_attributes(g,a2)			
			edge_labels[a2] += list(attrs.values())
		
		struct.update(struct_aux)
		data.append(g)
		y.append(mol_class)
	labels = {'nodes': node_labels, 'edges': edge_labels}		
	return data, y, struct, labels
	
"""
	Taken from graphfiles.py
	Reads a gxl file describing a graph
	- filename: The path and file containing the gxl description of the file
	returns:
	- g: a nx graph object
	- struct: the node and edge structures. It is a dictionary with the folowwing fields:
		nodes: a dictionary having the name and tag of each attribute of a node. For example, for the LETTER dataset, the structure is
			<node><attr name="x"><float>1.2</float></attr><attr name="y"><float>2.2</float></attr></node>
		In this case, struct['nodes'] = {"x": "float", "y": "float"}
		edges: a dictionary having the name and tag of each attribute of an edge.
	- labels: a dictionary with the same structure as struct but having as values a list of all the values found for each attribute
"""
def loadGXL(filename):
	from os.path import basename
	import networkx as nx
	import xml.etree.ElementTree as ET

	tree = ET.parse(filename)
	root = tree.getroot()
	index = 0
	g = nx.Graph(**root[0].attrib)
	dic = {}  # used to retrieve incident nodes of edges
	node_attr = {}
	edge_attr = {}
	done = False
	for node in root.iter('node'):
		dic[node.attrib['id']] = index
		labels = {}
		for attr in node.iter('attr'):
			labels[attr.attrib['name']] = attr[0].text
			if not done:
				for a in attr:
					node_attr[attr.attrib['name']] = a.tag
		done = True
		g.add_node(node.attrib['id'], **labels)
		index += 1
	done = False
	for edge in root.iter('edge'):
		labels = {}
		for attr in edge.iter('attr'):
			labels[attr.attrib['name']] = attr[0].text
			if not done:
				for a in attr:
					edge_attr[attr.attrib['name']] = a.tag
		done = True

		g.add_edge(edge.attrib['from'], edge.attrib['to'], **labels)
	struct = {'nodes': node_attr, 'edges': edge_attr}
	return g, struct
	

def all_match(dict1, dict2, verbose = True):
	for k in dict1.keys():
		if dict1[k].strip()=="default":
			continue
		if k not in dict2:
			return False
		else:
			if dict1[k].strip() != dict2[k].strip():
				return False
	for k in dict2.keys():
		if dict2[k].strip()=="default":		
			continue
		if k not in dict1:
			return False
		else:
			if dict1[k].strip() != dict2[k].strip():
				return False
	return True
		
def test_iso(filename, original_path, decoded_path):
	data, y, struct, labels = loadFromXML(filename,original_path)
	data_decoded, y_decoded, struct, labels = loadFromXML(filename,decoded_path) 
	iso = False
	print("Sizes: Orig -> ", len(data), ", compressed -> ", len(data_decoded))

	
		
	good = True
	not_iso = 0
	for i, g in enumerate(data):
		G1 = g
		G2 = data_decoded[i]
		iso = nx.is_isomorphic(G1, G2, node_match=all_match, edge_match=all_match)
		
		if not iso: 
			print(G1.name + " vs " + G2.name + " = " + str(iso))
			not_iso +=1
			good = False
	
	print(filename)
	if good:
		print(" -- All good --")
	else:
		print(" -- PROBLEMS: check in the terminal to see graphs that are not isomorphic to their decompressed versions --")
		print(not_iso, " of ", len(data), " non-isomorphic")


def get_graph_dir(ds, base):
	if ds=="acyclic" or ds=="mao" or ds=="pah" or ds=="msts_int_w" or ds=="msts_float_w" or ds=="msts_no_w":
		return base + ds
	
	if ds=="AIDS" or ds=="Mutagenicity" or ds=="Protein":
		return base + ds + "/data"
	
	return base + ds + "/MED"
	

def get_collection_file(ds, base):
	return base + ds +".xml"
	
import argparse

parser = argparse.ArgumentParser(description="Check every pair of graphs with the same name in two collections are isomorphic.")
parser.add_argument("--base_xml", help="path to the folder where the original xml collection files", default="../../data/collections/")
parser.add_argument("--graph_dir_orig", help="path to the folder containig the original data sets with the graphs",default="../../data/datasets/")
parser.add_argument("--graph_dir_decoded", help="path to the folder containig the decoded data sets with the graphs",default="../data/output/")
parser.add_argument("--suffix_decoded", help="name of the folder where the decoded files are",default="/decoded_bin")
parser.add_argument("--datasets", help="names of the datasets. Can be many, separated by spaces",nargs="+")

args = parser.parse_args()

for d in args.datasets:
	test_iso(get_collection_file(d, args.base_xml), get_graph_dir(d, args.graph_dir_orig), args.graph_dir_decoded + d + args.suffix_decoded)

						
	

	
	
	

	
	
	
	

