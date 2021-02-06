#!/usr/bin/env python3

import glob
import sys
import os 
import pandas as pd 
import numpy as np  
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
import time
import networkx as nx
import matplotlib.pyplot as plt
from pyvis.network import Network

from scipy.stats import zscore
# def replaceN(x):
# 	if x == "\\N":
# 		return "NaN"
# 	return x 

def checkDelimFromLine(line):
	delim = ","
	if "\t" in line:
		delim = "\t"
	elif "," in line:
		delim = ","
	elif ";" in line:
		delim = ";"
	elif " " in line:
		delim = " "
	return delim

def loadCsv(testFile,observedTickers=""):
	delim = ","
	with open(testFile) as ff:

		for line in ff:
	# 		line = line.replace('\\N',"NaN")
			if "\t" in line:
				delim = "\t"
			elif "," in line:
				delim = ","
			elif ";" in line:
				delim = ";"
			elif " " in line:
				delim = " "
			break


	# df = pd.read_csv(testFile, index_col = None, header=0, delimiter=delim)
	df = pd.read_csv(testFile, index_col = 0, header=0, delimiter=delim, dtype='unicode')
	# print(df)

	missing_tickers = ""
	if observedTickers != "" :
		observedTickersArr = set([x for x in observedTickers.split(",") if x.strip() != ""])
		original_nodenames = df.columns.values.tolist()
		missing_tickers = set(observedTickersArr) - set(original_nodenames)
		# observedTickersArr = list(observedTickersArr - missing_tickers)
		selected_tickers = [x for x in original_nodenames if x not in missing_tickers and x in observedTickersArr]
		# print(selected_tickers)
		df = df[selected_tickers]

	return df, ",".join(list(missing_tickers))

def formatDF(df):
	
	
	# print(df.columns.values)
	# print(observedTickersArr)
	# sys.exit()
	# print(missing_tickers)
	# print(df)
	# sys.exit()
	# print(df)
	# sys.exit()
	df = df.replace('\\N',np.NaN)
	df = df.replace('NaN',np.NaN)
	df = df.astype(float)
	# how to remove non-numeric values ?

	

	return df

def dataFrameFromFile(testFile,observedTickers=""):
	# content = []
	delim = ","
	with open(testFile) as ff:

		for line in ff:
	# 		line = line.replace('\\N',"NaN")
			if "\t" in line:
				delim = "\t"
			elif "," in line:
				delim = ","
			elif ";" in line:
				delim = ";"
			elif " " in line:
				delim = " "
			break


	# df = pd.read_csv(testFile, index_col = None, header=0, delimiter=delim)
	df = pd.read_csv(testFile, index_col = 0, header=0, delimiter=delim)

	# classes = None
	# if "classes" in df.index.values:
	# 	classes = df.loc["classes"]
	# 	df = df.drop("classes")


	observedTickersArr = set([x for x in observedTickers.split(",") if x.strip() != ""])
	missing_tickers = ""
	# print(df.columns.values)
	# print(observedTickersArr)
	# sys.exit()
	# print(missing_tickers)
	if observedTickers != "" :
		original_nodenames = df.columns.values.tolist()
		missing_tickers = set(observedTickersArr) - set(original_nodenames)
		# observedTickersArr = list(observedTickersArr - missing_tickers)
		selected_tickers = [x for x in original_nodenames if x not in missing_tickers and x in observedTickersArr]
		# print(selected_tickers)
		df = df[selected_tickers]

	# print(df)
	# sys.exit()
	df = df.replace('\\N',np.NaN)
	df = df.replace('NaN',np.NaN)
	df = df.astype(float)
	# how to remove non-numeric values ?

	return df, ",".join(list(missing_tickers))

def zNorm(df):
	for col in df.columns.values:
		df[col] = (df[col] - df[col].mean())/df[col].std(ddof=0)

	# df.apply(zscore)
	return df

def converVal(x, edgeThreshold):
	val = 2-x
	if val < edgeThreshold:
		val = 0.0
	return val


def corrToDistance(corr_input,edgeThreshold=0.0,ifIncludeNegativeEdges=True):

	corr = corr_input.copy()
	
	# distance formula
	corr = corr.apply(lambda x:2-abs(x)) if ifIncludeNegativeEdges else corr.apply(lambda x:2-x)
	# print(ifIncludeNegativeEdges,corr)
	# sys.exit()
	# make unqualified pairs unreachable
	corr[corr > 2 - edgeThreshold] = 0
	corr = corr.fillna(0)
	
	# buttom part to 0
	matrix_unsymmetric = corr.values.tolist()
	for idx,row in enumerate(matrix_unsymmetric):
		for x in range(0,idx):
			matrix_unsymmetric[idx][x] = 0.0

	# diagonal to 0
	for i in range(len(matrix_unsymmetric)):
		matrix_unsymmetric[i][i] = 0.0

	return matrix_unsymmetric

def addToGraphSimple(arr):
	
	G=nx.DiGraph()
	for x in range(len(arr)):
		nd = x
		G.add_node(nd)
	# show the MST result
	for x in range(len(arr)):
		for y in range(len(arr[x])):
			if arr[x][y] != 0:
				G.add_edge(x,y)

	g=Network(height=800,width=800,notebook=True)
	g.from_nx(G)
	return g

def arrToNX(arr):
	G=nx.Graph()
	edges = []
	for x in range(len(arr)):
		G.add_node(x)
		for y in range(len(arr[x])):
			if float(arr[x][y]) != 0.0:
				edges.append([x,y,round(2-arr[x][y],2)])

	
	for edge in edges:
		x = edge[0]
		y = edge[1]
		G.add_edge(x,y)
	return G

def toJson2(corr):
	print('in')
	return 1,2
def toJson(corr, arr, names , trim=None ):
	orig_corr = corr.values.tolist()
	edges = []

	for x in range(len(arr)):
		for y in range(len(arr[x])):
			# print(arr[x][y])
			if float(arr[x][y]) != 0.0: # 0.0 are unreachable pair points

				# edges.append([x,y,round(2-arr[x][y],2)])
				edges.append([x,y,round(orig_corr[x][y],2)]) 
				# print(x,y,orig_corr[x][y])
				# relevant_nodes.append(x)
				# relevant_nodes.append(y)



	edges.sort(key=lambda x: x[2],reverse=True)

	if trim:
		edges = edges[:trim]

	relevant_nodes = []
	for e in edges:
		relevant_nodes.append(e[0])
		relevant_nodes.append(e[1])

	if len(names) != 0:
		relevant_one = {x:names[x] for x in relevant_nodes}
	else:
		relevant_one = {x:x for x in relevant_nodes}



	return relevant_one,edges


def addToColorGraph(arr,size=800, tickerfile="demo_data/Tickers_Tong.csv"):
	names = []
	if tickerfile != "":
		df = pd.read_csv(tickerfile,header=-1)
		names = df[0].tolist()


	edges = []
	relevant_nodes = []
	for x in range(len(arr)):
		for y in range(len(arr[x])):
			if arr[x][y] != 0:
				edges.append([x,y,arr[x][y]])
				relevant_nodes.append(x)
				relevant_nodes.append(y)

	relevant_nodes = set(relevant_nodes)
	g=Network(height=size,width=size,notebook=True)
	for x in range(len(arr)):
		if x not in relevant_nodes:
			continue
		nd = x if len(names)==0 else names[x]
		nd_name = str(nd) if len(names)==0 else names[x]
		g.add_node(nd,title=nd_name)

    # add tags
	# g.add_nodes([1,2,3], value=[10, 100, 400], title=["I am node 1", "node 2 here", "and im node 3"], x=[21.4, 54.2, 11.2], y=[100.2, 23.54, 32.1], label=["NODE 1", "NODE 2", "NODE 3"], color=["#00ff1e", "#162347", "#dd4b39"])

	for edge in edges:
		val = 2-edge[2]
		x = edge[0]
		y = edge[1]
		if len(names)!=0:
			x = names[edge[0]]
			y = names[edge[1]]


		color='gray'
		weight = 5
		if val >= 0.5:
			color = 'red'
			weight = 20
		elif val > 0.3 and val < 0.5:
			color = 'orange'
			weight = 10
		g.add_edge(x,y,color=color,label = str(val), width=weight)

	# for x in range(len(arr)):
	# 	for y in range(len(arr[x])):
	# 		if arr[x][y] != 0:
	# 			val = 2-arr[x][y]
	# 			color='gray'
	# 			weight = 5
	# 			if val >= 0.5:
	# 				color = 'red'
	# 				weight = 20
	# 			elif val > 0.3 and val < 0.5:
	# 				color = 'orange'
	# 				weight = 10
	# 			g.add_edge(x,y,color=color,label = str(val), width=weight)


	g.toggle_hide_edges_on_drag(True)
	g.barnes_hut()
	return g



def subFrame(df,start,end):
	reduceDF = df.iloc[start:end,:]
	return reduceDF

def unsymmetricMatrixToMST(matrix_unsymmetric):
	# X = csr_matrix(matrix_unsymmetric)
	# X = matrix_unsymmetric
	Tcsr = minimum_spanning_tree(matrix_unsymmetric)
	arr= Tcsr.toarray().astype(float)
	return arr

# def unsymmetricMatrixToMST(matrix_unsymmetric):
# 	# X = csr_matrix(matrix_unsymmetric)
# 	# X = matrix_unsymmetric
# 	Tcsr = minimum_spanning_tree(matrix_unsymmetric)
# 	arr= Tcsr.toarray().astype(float)
# 	return arr


def main(args):
	

	time_start = time.time()

	testFile = args[0]

	start = -1
	end = -1
	if len(args) > 1:
		start,end = int(args[1]),int(args[2])


	filename = testFile.split("/")[-1].split(".")[0]

	# df = df.fillna(method='ffill')
	df = dataFrameFromFile(testFile)

	# tt(len(df))
	# sys.exit()

	if start !=  -1 and end != -1:
		df  = subFrame(df,start,end)

	# z norm
	df  = zNorm(df)

	# corr
	corr = df.corr(method='pearson',min_periods=30) # will be personalized with input

	matrix_unsymmetric = corrToDistance(corr)

	# tree
	arr = unsymmetricMatrixToMST(matrix_unsymmetric)
	
	g = addToColorGraph(arr,size=500)

	# g.show("tmp/"+filename+".html")
	g.show("demo_data/"+filename+".html")

	# time_lapse = time.time()-time_start
	# print("---\ntime span for total correlation",round(time_lapse,6))
	
if __name__ == '__main__':
  main(sys.argv[1:])
# print('HI from python')