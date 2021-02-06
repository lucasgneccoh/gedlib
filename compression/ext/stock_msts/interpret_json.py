import json
import numpy as np
import networkx as nx
import graphs_gxl as gxl
from progress.bar import Bar
import os
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

def parseInputs():
	parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
	parser.add_argument("--json_file", help="json file returned by mstEvolve.py containing the msts information")
	parser.add_argument("--out_folder", help="path to store the collections of graphs")
	parser.add_argument("--suffix", help="suffix to append to the dataset name", default="")
	args = parser.parse_args()
	return args


# START

args = parseInputs()
preffix = args.out_folder
if preffix[-1] != '/': preffix += '/'

with open(args.json_file) as file:
	graphs = json.load(file)


nx_graphs = []
y = []
bar = Bar('Float weight', max=len(graphs))
for i, g in enumerate(graphs):
	edge_list = [(e[0], e[1], {'w':e[2]}) for e in g]
	g_nx = nx.Graph(edge_list)
	g_nx.graph = {'id': str(i)}
	for n in g_nx.nodes():
		g_nx.nodes[n]["stock"] = n
	nx_graphs.append(g_nx)
	y.append('0')
	bar.next()
	
bar.finish()


import matplotlib.pyplot as plt
plt.plot([g.number_of_nodes() for g in nx_graphs])
plt.show()

filename = preffix + "msts_float_w" + args.suffix + ".xml"
folder = preffix + "msts_float_w" + args.suffix
os.makedirs(folder, exist_ok=True)

struct = {'nodes':{'stock': 'string'}, 'edges':{'w':'float'}}
gxl.saveToXML(nx_graphs, y, filename=filename, graph_dir=folder, struct = struct)


nx_graphs = []
y = []
bar = Bar('Int weight', max=len(graphs))
for i,g in enumerate(graphs):
	edge_list = [(e[0], e[1], {'w':str(int(100*float(e[2])))}) for e in g]
	g_nx = nx.Graph(edge_list)
	g_nx.graph = {'id': str(i)}
	for n in g_nx.nodes():
		g_nx.nodes[n]["stock"] = n
	nx_graphs.append(g_nx)
	y.append('0')
	bar.next()
	
bar.finish()

filename = preffix + "msts_int_w" + args.suffix + ".xml"
folder = preffix + "msts_int_w" + args.suffix
os.makedirs(folder, exist_ok=True)
struct = {'nodes':{'stock': 'string'}, 'edges':{'w':'int'}}
gxl.saveToXML(nx_graphs, y, filename=filename, graph_dir=folder, struct = struct)



nx_graphs = []
y = []
bar = Bar('No weight', max=len(graphs))
for i,g in enumerate(graphs):
	edge_list = [(e[0], e[1]) for e in g]
	g_nx = nx.Graph(edge_list)
	g_nx.graph = {'id': str(i)}
	for n in g_nx.nodes():
		g_nx.nodes[n]["stock"] = n
	nx_graphs.append(g_nx)
	y.append('0')
	bar.next()
	
bar.finish()

filename = preffix + "msts_no_w" + args.suffix + ".xml"
folder = preffix + "msts_no_w" + args.suffix
os.makedirs(folder, exist_ok=True)
struct = {'nodes':{'stock': 'string'}, 'edges':{}}
gxl.saveToXML(nx_graphs, y, filename=filename, graph_dir=folder, struct = struct)

