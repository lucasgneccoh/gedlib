import shutil
from os import listdir, mkdir
from os.path import isfile, join
import os.path
import os
import xml.etree.ElementTree as ET
import sys
import networkx as nx
import numpy as np

def get_graph_dir(ds, base):
	if ds=="acyclic" or ds=="mao" or ds=="pah" or "msts" in ds:
		return base + ds
	
	if ds=="AIDS" or ds=="Mutagenicity" or ds=="Protein":
		return base + ds + "/data"
	
	return base + ds + "/MED"
	

def get_collection_file(ds, base):
	return base + ds +".xml"

	
import argparse

parser = argparse.ArgumentParser(description="Check every pair of graphs with the same name in two collections are isomorphic.")
parser.add_argument("--target", help="target folder to create the copies", default="../data/orig_datasets_to_tar")
parser.add_argument("--base_xml", help="path to the folder containig the original xml collection files",default="../../data/collections/")
parser.add_argument("--base_datasets", help="path to the folder containig the original data sets with the graphs",default="../../data/datasets/")

parser.add_argument("--datasets", help="names of the datasets. Can be many, separated by spaces",nargs="+")

args = parser.parse_args()

for d in args.datasets:
	print(d, end="...")
	filename = get_collection_file(d, args.base_xml)
	original_path = get_graph_dir(d, args.base_datasets)
	new_path = args.target + "/" + d
	try:
		shutil.rmtree(new_path)
	except:
		pass
		
	try:
		mkdir(new_path)
	except:
		pass
	# Read all graphs and copy them
	tree = ET.parse(filename)
	root = tree.getroot()
	for graph in root.iter('graph'):
		g_name = graph.attrib['file']
		orig_g_path = original_path + "/" + g_name
		shutil.copy(orig_g_path, new_path)
	shutil.copy(filename, new_path)
	print("done")

