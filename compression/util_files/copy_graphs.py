import shutil
from os import listdir
from os.path import isfile, join
import os.path
import xml.etree.ElementTree as ET
import sys
import networkx as nx
import numpy as np

def get_graph_dir(ds, base):
	if ds=="acyclic" or ds=="mao" or ds=="pah":
		return base + ds
	
	if ds=="AIDS" or ds=="Mutagenicity" or ds=="Protein":
		return base + ds + "/data"
	
	return base + ds + "/MED"
	

def get_collection_file(ds, base):
	return base + ds +".xml"


datasets = ["acyclic"] #, "AIDS", "Letter", "mao", "Mutagenicity", "pah", "Protein"]

destiny =  "../data/orig_datasets_to_tar"
base_xml = "../../data/collections/"
base_datasets = "../../data/datasets/"
for d in datasets:
	print(d)
	filename = get_collection_file(d, base_xml)
	original_path = get_graph_dir(d, base_datasets)
	# Read all graphs and copy them
	tree = ET.parse(filename)
	root = tree.getroot()
	for graph in root.iter('graph'):
		g_name = graph.attrib['file']
		orig_g_path = original_path + "/" + g_name
		new_path = destiny + "/" + d
		shutil.copy(orig_g_path, new_path)
	print("done")

