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
#		if 'chem' in labels:
#			labels['label'] = labels['chem']
#			node_attr['label'] = 'string'
#			labels['atom'] = labels['chem']
#			node_attr['atom'] = 'string'
		#g.add_node(index, **labels)
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
#		if 'valence' in labels:
#			labels['label'] = labels['valence']
#			edge_attr['label'] = 'int' 
#			labels['bond_type'] = labels['valence']
#			edge_attr['bond_type'] = 'int' 
		#g.add_edge(dic[edge.attrib['from']], dic[edge.attrib['to']], **labels)
		g.add_edge(edge.attrib['from'], edge.attrib['to'], **labels)
	struct = {'nodes': node_attr, 'edges': edge_attr}
	return g, struct
	
"""
	Taken from graphfiles.py
	Saves a gxl file describing a graph
	parameters:
	- g: a nx graph object
	- filename: path and file name to write the graph into.
	- struct: it is recommended to always specify the node and edge structure to write the file. See loadGXL for the descriptions of the struct dictionary. 
"""
def saveGXL(g,filename, struct = None):
	with open(filename, 'w') as gxl_file:
		gxl_file.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
		gxl_file.write("<!DOCTYPE gxl SYSTEM \"http://www.gupro.de/GXL/gxl-1.0.dtd\">\n")
		gxl_file.write("<gxl xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n")

		gxl_file.write("<graph ")
		for k in g.graph:
			gxl_file.write(k + "=\"" + str(g.graph[k]) + "\" ")
		gxl_file.write(">\n")
		if not struct is None:
			node_attr = struct['nodes']
			edge_attr = struct['edges']
			tag = True
		else:
			node_attr = None
			edge_attr = None
			tag = False
			
		for v, attrs in g.nodes(data=True):
			gxl_file.write("<node id=\"" + str(v) + "\">")
			for a in attrs:	
				gxl_file.write("<attr name=\"" + a + "\">")
				if tag: gxl_file.write("<" + node_attr[a] + ">")
				gxl_file.write(str(attrs[a]))
				if tag: gxl_file.write("</" + node_attr[a] + ">")
				gxl_file.write("</attr>")
			gxl_file.write("</node>\n")
			
		for v1, v2, attrs in g.edges(data=True):
			gxl_file.write("<edge from=\"" + str(v1) + "\" to=\"" + str(v2) + "\">")
			for a in attrs:	
				gxl_file.write("<attr name=\"" + a + "\">")
				if tag: gxl_file.write("<" + edge_attr[a] + ">")
				gxl_file.write(str(attrs[a]))
				if tag: gxl_file.write("</" + edge_attr[a] + ">")
				gxl_file.write("</attr>")
			gxl_file.write("</edge>\n")
		gxl_file.write("</graph>\n")
		gxl_file.write("</gxl>")
		gxl_file.close()
	

"""
	Taken from graphfiles.py
	Saves a gxl collection file and the corresponding graph files
	parameters:
	- graphs: a list of nx graph objects to save
	- y: a list with the classes corresponding to the graphs
	- filename: path and file name to write the collection file into.
	- graph_dir: If the graphs are not to be stored in the same folder where the collection file is, this variable specifies the folder to save the graph files into.
	- struct: it is recommended to always specify the node and edge structure to write the file. See loadGXL for the descriptions of the struct dictionary. 
"""
def saveToXML(graphs, y, filename='gfile', graph_dir=None, struct = None):
	import os
	dirname_ds = os.path.dirname(filename)
	if dirname_ds != '':
		dirname_ds += '/'
		if not os.path.exists(dirname_ds) :
			os.makedirs(dirname_ds)
				
	if graph_dir is not None:
		graph_dir = graph_dir + '/'
		if not os.path.exists(graph_dir):
			os.makedirs(graph_dir)
	else:
		graph_dir = dirname_ds 
		
	with open(filename, 'w') as fgroup:
		fgroup.write("<?xml version=\"1.0\"?>")
		fgroup.write("\n<!DOCTYPE GraphCollection SYSTEM \"http://www.inf.unibz.it/~blumenthal/dtd/GraphCollection.dtd\">")
		fgroup.write("\n<GraphCollection>")
		for idx, g in enumerate(graphs):
			fname_tmp = "graph" + str(idx) + ".gxl"
			saveGXL(g, graph_dir + fname_tmp, struct=struct)
			fgroup.write("\n\t<graph file=\"" + fname_tmp + "\" class=\"" + str(y[idx]) + "\"/>")
		fgroup.write("\n</GraphCollection>")
		fgroup.close()


"""
	Prints a representation of the graph including the ID and attributes of the graph itself, and the ID and attributes of all nodes and edges
	parameters:
	- graph: a nx graph object
"""
def visual(graph):
	print("GRAPH")
	print(graph.graph)
	print("NODES:")
	for n in graph.nodes:
		attr = graph.nodes[n]
		print(n, " -> ", attr)
	print("EDGES:")
	for e in graph.edges:
		attr = graph.edges[e]
		print(e, " -> ", attr)   


"""
	Draws the letter graph in the matplotlib axes object ax
	parameters:
	- graph: a nx graph object
	- ax: a matplotlib.Axes object
"""
def draw_Letter(graph, ax):
	import numpy as np
	import networkx as nx
	import matplotlib.pyplot as plt
	pos = {}
	for n in graph.nodes:
		pos[n] = np.array([float(graph.nodes[n]['x']),float(graph.nodes[n]['y'])])
	nx.draw_networkx(graph,pos, ax = ax)
	
def draw_lucas(graph, ax):
	import numpy as np
	import networkx as nx
	import matplotlib.pyplot as plt
	pos = {}
	labels = {}
	edge_color = []
	for n in graph.nodes:
		labels[n] = graph.nodes[n]['letra']
	for e in graph.edges:
		edge_color.append(graph.edges[e]['color'])
		
	nx.draw_networkx(graph,pos = nx.circular_layout(graph), ax = ax, labels = labels, edge_color = edge_color)
		
def test_iso(filename, original_path, decoded_path, verbose=False):
	data, y, struct, labels = loadFromXML(filename,original_path)
	data_decoded, y_decoded, struct, labels = loadFromXML(filename,decoded_path) 
	iso = False
	print("Sizes: Orig -> ", len(data), ", compressed -> ", len(data_decoded))

	def all_match(dict1, dict2, verbose = True):
		for k in dict1.keys():
			if dict1[k].strip()=="default":
				continue
			if k not in dict2:
				#if(verbose): print(k, " not in dict2")
				return False
			else:
				if dict1[k].strip() != dict2[k].strip():
					#if(verbose): print(k, " : ", dict1[k], " != ", dict2[k])
					return False
		for k in dict2.keys():
			if dict2[k].strip()=="default":		
				continue
			if k not in dict1:
				#if(verbose): print(k, " not in dict1")
				return False
			else:
				if dict1[k].strip() != dict2[k].strip():
					#if(verbose): print(k, " : ", dict1[k], " != ", dict2[k])
					return False
		return True
		
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
	
		
if __name__ == "__main__":

	import matplotlib.pyplot as plt
	# "/home/lucas/Documents" -> "/home/lamsade/lgnecco"
	base = "/home/lamsade/lgnecco"
	#base = "/home/lucas/Documents"
	
	suffix = "AIDS"
	verbose = True
	
	filename_input = base + "/stage/gedlib/compression/data/test_collections/collections_list_" + suffix + ".txt"
	original_path_input = base +  "/stage/gedlib/compression/data/test_collections/graphs_dir_list_"+ suffix +".txt"
	dataset_input= base + "/stage/gedlib/compression/data/test_collections/dataset_"+ suffix +".txt"

	output_root_input= base + "/stage/gedlib/compression/data/test_collections/output_root.txt"

	
	with open(filename_input, "r") as col:
		with open(original_path_input, "r") as g_dir:
			with open(output_root_input, "r") as out_root:
				with open(dataset_input, "r") as dataset:
					filename = col.readline().rstrip('\n')
					while filename:
						original_path = g_dir.readline().rstrip('\n')
						output_root = out_root.readline().rstrip('\n')
						dataset_name = dataset.readline().rstrip('\n')
					
						decoded_path = output_root + "/" + dataset_name + "/decoded"
						test_iso(filename, original_path, decoded_path, verbose)
						filename = col.readline().rstrip('\n')
					col.close()
					g_dir.close()
					out_root.close()
					dataset.close()
						
	

	
	
	

	
	
	
	

