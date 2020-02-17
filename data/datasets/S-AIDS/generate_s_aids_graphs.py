#//////////////////////////////////////////////////////////////////////////#
#                                                                          #
#   Copyright (C) 2020 by David B. Blumenthal                              #
#                                                                          #
#   This file is part of GEDLIB.                                           #
#                                                                          #
#   GEDLIB is free software: you can redistribute it and/or modify it      #
#   under the terms of the GNU Lesser General Public License as published  #
#   by the Free Software Foundation, either version 3 of the License, or   #
#   (at your option) any later version.                                    #
#                                                                          #
#   GEDLIB is distributed in the hope that it will be useful,              #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of         #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the           #
#   GNU Lesser General Public License for more details.                    #
#                                                                          #
#   You should have received a copy of the GNU Lesser General Public       #
#   License along with GEDLIB. If not, see <http://www.gnu.org/licenses/>. #
#                                                                          #
#//////////////////////////////////////////////////////////////////////////#

##
# @file generate_molecules.py
# @brief Python script for generating synthetic molecules.
# @details The synthetic molecules generated by this script were used for the experiments in the following paper:
# - D. B. Blumenthal, S. Bougleux, N. Boria, J. Gamper, L. Brun:
#   &ldquo;Comparing heuristics for graph edit distance computation&rdquo;,
#   Accepted for publication in VLDB J. 
# 
# Usage: 
# ```sh
# $ python generate molecules.py
# ```
#
# @warning Running this script overrides the molecules distributed with GEDLIB that were used for the experiments in the VLDB J. paper.
'''
Python script for generating synthetic graphs with the topology of AIDS graphs but fewer node labels.
'''

import random
import xml.etree.ElementTree as ET
import argparse
import os.path
import subprocess

def parse_graph(dir, gxl_file):
    num_nodes = 0
    graph = ET.parse(os.path.join(dir, gxl_file)).getroot()
    str_to_id = {}
    for node in graph.findall("graph/node"):
        str_to_id.update({node.attrib["id"] : num_nodes})    
        num_nodes = num_nodes + 1
    edge_list = []
    adj_list = {node_id : [] for node_id in range(num_nodes)}
    adj_matrix = [[0 for col in range(num_nodes)] for row in range(num_nodes)]
    for edge in graph.findall("graph/edge"):
        tail = str_to_id[edge.attrib["from"]]
        head = str_to_id[edge.attrib["to"]]
        if adj_matrix[tail][head] == 1:
            continue
        adj_matrix[tail][head] = 1
        adj_matrix[head][tail] = 1
        edge_list.append((tail,head))
        adj_list[tail].append(head)
        adj_list[head].append(tail)
    return edge_list, adj_list, adj_matrix, num_nodes, len(edge_list)

class Graph:
    
    def __init__(self, dir, filename):
        # Initialize member variables.
        self.filename = filename
        self.nodes = []
        self.edges = []
        self.node_alphabet = []
        self.num_node_labels = 0
        self.node_label_ids = {}
        
        # Parse the nodes.
        node_id = 0
        graph = ET.parse(os.path.join(dir, filename)).getroot()
        str_to_id = {}
        for node in graph.findall("graph/node"):
            str_to_id[node.attrib["id"]] = node_id
            node_label = ""
            for attr in node.findall("attr"):
                if attr.get("name") == "chem":
                    node_label = attr.find("int").text
                    break
            if node_label == "":
                raise Exception("node " + str(node_id) + " does not have attribute chem.")
            self.nodes.append((node_id, node_label))
            self.node_alphabet.append(node_label)
            node_id += 1
        
        # Parse the edges.
        adj_matrix = [[0 for col in range(len(self.nodes))] for row in range(len(self.nodes))]
        for edge in graph.findall("graph/edge"):
            tail = str_to_id[edge.attrib["from"]]
            head = str_to_id[edge.attrib["to"]]
            if adj_matrix[tail][head] == 1:
                continue
            adj_matrix[tail][head] = 1
            adj_matrix[head][tail] = 1
            self.edges.append((tail, head, edge.find("attr/int").text))
        
        
    def shrink_node_label_alphabet(self):
        self.num_node_labels -= 1
        if self.num_node_labels > 0:
            for node in range(len(self.nodes)):
                label = self.nodes[node][1]
                if self.node_label_ids[label] >= self.num_node_labels:
                    label_id = random.randint(0, self.num_node_labels - 1)
                    self.nodes[node] = (self.nodes[node][0], self.node_alphabet[label_id])
    
    def __repr__(self):
        string = "num_nodes =  " + str(len(self.nodes)) + "\n"
        string = string + "nodes =  " + str(self.nodes) + "\n"
        string = string + "edges =  " + str(self.edges)
        return string
            
    def to_gxl(self, directory):
        gxl_file_name = os.path.join(directory, self.filename)
        gxl_file = open(gxl_file_name, "w")
        gxl_file.write("<?xml version=\"1.0\"?>\n")
        gxl_file.write("<!DOCTYPE gxl SYSTEM \"http://www.gupro.de/GXL/gxl-1.0.dtd\">\n")
        gxl_file.write("<gxl>\n")
        gxl_file.write("<graph id=\"" + self.filename + "\" edgeids=\"false\" edgemode=\"undirected\">\n")
        for id, label in self.nodes:
            gxl_file.write("<node id=\"_" + str(id) + "\">\n")
            gxl_file.write("<attr name=\"chem\"><int>" + str(label) + "</int></attr>\n")
            gxl_file.write("</node>\n")
        for edge in self.edges:
            gxl_file.write("<edge from=\"_" + str(edge[0]) + "\" to=\"_" + str(edge[1]) + "\">\n")
            gxl_file.write("<attr name=\"valence\"><int>" + str(edge[2]) + "</int></attr>\n")
            gxl_file.write("</edge>\n")
        gxl_file.write("</graph>\n")
        gxl_file.write("</gxl>\n")
        gxl_file.close()
        

collection = ET.parse("../../collections/S-AIDS.xml").getroot()
graphs = [Graph("../AIDS/data/", entry.attrib["file"]) for entry in collection]
node_alphabet = []
for graph in graphs:
    node_alphabet += graph.node_alphabet
node_alphabet = [label for label in set(node_alphabet)]
node_alphabet.sort()
random.shuffle(node_alphabet)
num_node_labels = len(node_alphabet)
node_label_ids = {node_alphabet[i] : i for i in range(num_node_labels)}
for graph in graphs:
    graph.node_alphabet = node_alphabet
    graph.node_label_ids = node_label_ids
    graph.num_node_labels = num_node_labels
while num_node_labels > 0:
    dir = "NL" + str(graphs[0].num_node_labels)
    subprocess.call("mkdir -p " + dir, shell=True)
    for graph in graphs:
        graph.to_gxl(dir)
        graph.shrink_node_label_alphabet()
    num_node_labels -= 1