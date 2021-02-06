"""
mst from inout timeseries data
example: 
python3 mstweb.py --csv demo_data/ita2001_2016.csv --ifCompleteGraph False

python3 mstweb.py --csv demo_data/ita2010_2020.csv
python3 mstweb.py --csv demo_data/ita2010_2020.csv --ifClassGraph True


python3 mstweb.py --csv demo_data/ita2001_2016_classified.csv --overlapThreshold -123


python3 mstweb.py --csv demo_data/longPaolo.csv


python3 mstweb.py --csv demo_data/ita2001_2016.csv --observedTickers  G.MI,ISP.MI,UCG.MI,FCA.MI,F.MI,ENI.MI,ENEL.MI,TIT.MI,BMED.MI,MB.MI
"""

import glob
import sys
import os 
import pandas as pd 
import numpy as np  
# from mst import addToColorGraph,corrToDistance,zNorm,subFrame,dataFrameFromFile,unsymmetricMatrixToMST,toJson,loadCsv,formatDF
from mst import *
import time
from pyvis.network import Network
import networkx as nx
from networkx.readwrite import json_graph
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import random
import os
import json
import collections
import networkx.algorithms.community as nx_comm

from sys import platform # system compatibility switch

def convertDataFrameToClasses(df,class_dic):
	symList = set(df.columns.values.tolist())

	# class_dic = classes.to_dict()
	class_symbols = {}
	# class_names = list(class_dic.keys())

	for tmp_sym in class_dic.keys():
		if tmp_sym not in symList:
			continue

		tmp_class = class_dic[tmp_sym]
		if tmp_class not in class_symbols:
			class_symbols[tmp_class] = [tmp_sym]
		else:
			class_symbols[tmp_class].append(tmp_sym)

	symbs = list(class_symbols.keys())

	for sym in symbs:
		df[sym] = df.loc[:,class_symbols[sym]].mean(axis=1)
	# df['total']=df.loc[:,list_name].sum(axis=1)


	df = df[symbs]

	return df

def analysis(edges,node_dic,class_dic={}):
	G=nx.Graph()
	

	for edge in edges:
		x = edge[0]
		y = edge[1]
		G.add_edge(x,y)
		# print(x,y)

	graph_degrees = G.degree()

	# print(platform)
	# print(graph_degrees)
	# print([y for x, y in graph_degrees.items()])
	# if platform == "linux" or platform == "linux2":
	# 	degree_sequence = sorted([y for x,y in G.degree().items()]) # ubuntu
	# elif platform == "darwin":
	# 	degree_sequence = sorted([d for n, d in G.degree()], reverse=False)  # mac
	# elif platform == "win32":
	# 	degree_sequence = sorted([d for n, d in G.degree()], reverse=False)  # mac default ?
	# else:
	degree_sequence = sorted([d for n, d in G.degree()], reverse=False)  # mac default ?
	# print(degree_sequence)
	# degree_sequence.sort()

	# print(degree_sequence)
	
	degreeCount = collections.Counter(degree_sequence)
	degreeDist = list(degreeCount.items())

	# print(degreeDist)
	# 
	# # remove disconnected nodes to calculate Radius and Graph centers.
	if not nx.is_empty(G):
		largest_cc = max(nx.connected_components(G))
		nodelist = G.nodes()
		to_remove = set(nodelist) - set(largest_cc)
		for nd in to_remove:
			G.remove_node(nd)


	radius,eccentricity,center = "","",""
	if not nx.is_empty(G):
		if nx.is_connected(G):
			radius = nx.radius(G)
			eccentricity = nx.eccentricity(G)
			center = nx.center(G) # The center is the set of nodes with eccentricity equal to radius.


	mod_score = -2
	if bool(class_dic.keys()):
		inv_dic = {}
		nodelist = G.nodes()
		for n in nodelist:
			converted_num = node_dic[n]
			class_tmp = class_dic[converted_num]
			if class_tmp not in inv_dic:
				inv_dic[class_tmp] = set([n])
			else:
				inv_dic[class_tmp].add(n)
		# print('hi')
		# G = nx.barbell_graph(3, 0)
		# print(G.degree())
		# mod_score = nx_comm.modularity(G, [{0, 1, 2}, {3, 4, 5}])
		# print('ms',mod_score)
		# print('okk',inv_dic)
		groups = [x for k,x in inv_dic.items()]
		# print(groups)
		mod_score = nx_comm.modularity(G, groups)
		# print(mod_score)

	# return 1,1,1,1,1
	return radius,eccentricity,center,degreeDist,mod_score


def elaborateData(df,inputCsv,start,end,observedTickers\
	,ifZNorm,presenceTreshold,ifClassGraph):

	debug_msg = ""
	debug_msg += "inputCsv: "+inputCsv+"; "

	hasClass = False
	class_dic = {}
	if "classes" in df.index.values:
		class_dic = df.loc["classes"].to_dict()
		df = df.drop("classes")
		hasClass = True
		debug_msg += "foundClasses: YES; "
	else:
		debug_msg += "foundClasses: NO; "

	df = formatDF(df)



	if start !=  -1 and end != -1:
		# df  = subFrame(df,start,end)
		df = df.iloc[start:end,:]
		debug_msg += "range: " + str(start)+"-"+str(end)+"; "
	elif start !=  -1: 		
		df = df.iloc[start:,:]
	elif end !=  -1: 
		df = df.iloc[:-end,:]
	else:
		debug_msg += "range: all; "

	if presenceTreshold:
		# further pruning since some symbols will be disqualify by the number of presence.
		size_pre = len(df.columns.values.tolist())
		# print(df)
		# print(presenceTreshold)
		
		df = delInsufficientTickers(df,presenceTreshold)
		# print(df)
		# sys.exit(111)
		size_after = len(df.columns.values.tolist())
		debug_msg += "pruned by percentage of NaNs: from"+str(size_pre)+" to "+str(size_after)+"; "
	else:
		debug_msg += "no pruning of NaNs"
		# if hasClass: # since some sym are eliminated, I have to adjust class_dic again.
		# 	new_class_dic = {}
		# 	for sym in df.columns.values:
		# 		new_class_dic[sym] = class_dic[sym]
		# 	class_dic = new_class_dic

	# z norm

	if ifZNorm != "no":
		# print('yes')
		df  = zNorm(df)

	
	# print(df)

	if hasClass and ifClassGraph:
		df = convertDataFrameToClasses(df,class_dic)
		debug_msg += "classGraph: Yes; "
	else:
		debug_msg += "classGraph: No; "

	# print(df)
	return df,hasClass,class_dic,debug_msg

def calculate_class_cover_ratio(edges,node_dic,class_dic):
	num = len(edges) if len(edges) !=0 else 1
	# print(num)
	count = 0
	for e in edges:
		n1 = class_dic[node_dic[int(e[0])]]
		n2 = class_dic[node_dic[int(e[1])]]
		if n1 == n2:
			count += 1

	return round(count/num,4)


def parseInputs():
	parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
	parser.add_argument("--csv", help="Input csv of time series")
	parser.add_argument("--start", help="start time point")
	parser.add_argument("--end", help="end time point")
	parser.add_argument("--edgeThreshold", help="Edge threshold")
	parser.add_argument("--observedTickers", help="Tickers to observe")
	parser.add_argument("--ifZNorm", help="Tickers to observe")
	# 13.05.2020
	parser.add_argument("--ifCompleteGraph", help="1 for complete 0 otherwise")
	parser.add_argument("--presenceTreshold", help="Percentage of available time points above which a symbol is qualified")
	parser.add_argument("--overlapThreshold", help="percentage of overlap for Pearson correlation")
	parser.add_argument("--classfile", help="A csv which defines the class of each ticker")
	parser.add_argument("--maxEdges", help="Max number of edges to show")
	parser.add_argument("--ifClassGraph", help="A graph of classes or symbols")

	# 13.06.2020
	parser.add_argument("--ifIncludeNegativeEdges", help="PLZ indicate true or false")

	args_ = parser.parse_args()

	inputCsv,start,end,edgeThreshold,observedTickers,ifZNorm,ifCompleteGraph,\
	presenceTreshold,classfile,maxEdges,ifClassGraph,overlapThreshold,ifIncludeNegativeEdges =\
	args_.csv,args_.start,args_.end,args_.edgeThreshold,args_.observedTickers,\
	args_.ifZNorm,args_.ifCompleteGraph,args_.presenceTreshold,args_.classfile,\
	args_.maxEdges,args_.ifClassGraph,args_.overlapThreshold,args_.ifIncludeNegativeEdges

	if None in [inputCsv]:
		print (__doc__)
		sys.exit(0)

	if not start:
		start = -1
	if not end:
		end = -1
	if not edgeThreshold:
		edgeThreshold = 0.0
	if not observedTickers:
		observedTickers = ""
	else:
		observedTickers = observedTickers.replace("[space]"," ").replace("[quote]","'")

	if presenceTreshold:
		presenceTreshold = float(presenceTreshold)


	if overlapThreshold:
		overlapThreshold = float(overlapThreshold)
		if  overlapThreshold <= 0.0 :
			overlapThreshold = 0.0
		elif overlapThreshold >= 1.0:
			overlapThreshold = 1.0
	else:
		overlapThreshold = 0.0

	if maxEdges:
		maxEdges = int(maxEdges)
	else:
		maxEdges = 1500

	if None not in [ifCompleteGraph]:
		ifCompleteGraph = True if ifCompleteGraph.lower() == "true" else False
	if None not in [ifClassGraph]:
		ifClassGraph = True if ifClassGraph.lower() == "true" else False
	if None not in [ifIncludeNegativeEdges]:
		ifIncludeNegativeEdges = True if ifIncludeNegativeEdges.lower() == "true" else False

	return inputCsv,int(start),int(end),float(edgeThreshold),\
	observedTickers,ifZNorm,ifCompleteGraph,presenceTreshold,\
	classfile,maxEdges,ifClassGraph,overlapThreshold,ifIncludeNegativeEdges


def delInsufficientTickers(df,presenceTreshold):
	ticker_ids = df.columns.values
	date_ids = df.index.values
	max_nam_allowed = df[ticker_ids[0]].size * (1-presenceTreshold)

	to_del = []
	for ticker in ticker_ids:
		num_nans = df[ticker].isna().sum() 
		if num_nans > max_nam_allowed:
			to_del.append(ticker)

	df=df.drop(columns=to_del)
	return df


def main(args):
	random.seed(0)

	inputCsv,start,end,edgeThreshold,observedTickers\
	,ifZNorm,ifCompleteGraph,presenceTreshold,\
	classfile,maxEdges,ifClassGraph,overlapThreshold,ifIncludeNegativeEdges = parseInputs()


	time_start = time.time()
	filename = inputCsv.split("/")[-1].split(".")[0]


	# df = df.fillna(method='ffill')

	# df,missing_ticker = dataFrameFromFile(inputCsv,observedTickers=observedTickers)
	df,missing_ticker = loadCsv(inputCsv,observedTickers=observedTickers)


	# print(df)
	# # sys.exit()
	end = end + 1 # date adjustment
	df,hasClass,class_dic,debug_msg = elaborateData(df,inputCsv,start,end,observedTickers\
	,ifZNorm,presenceTreshold,ifClassGraph)

	# print(df,ifZNorm,presenceTreshold)
	# # sys.exit()
	# # corr and mst

	min_periods = round(len(df) * overlapThreshold) if round(len(df) * overlapThreshold) != 0 else 1

	# # print(min_periods)
	# # sys.exit()
	corr = df.corr(method='pearson',min_periods=min_periods) # df, will be personalized with input
	
	# # testing!!
	# # corr['ch1_return_tank_level']['ch1_pump_speed'] = -0.99


	matrix_unsymmetric = corrToDistance(corr,edgeThreshold=edgeThreshold,ifIncludeNegativeEdges=ifIncludeNegativeEdges) # array like
	# print(matrix_unsymmetric)
	if ifCompleteGraph: # full network
		arr = matrix_unsymmetric
		debug_msg += "completeGraph: Yes; "
	else: # MST
		arr = unsymmetricMatrixToMST(matrix_unsymmetric) # df
		debug_msg += "completeGraph: No; "

	# print(corr)
	names = corr.columns.values
	# print(names)
	# print('Pre json')
	relevant_one,edges = toJson(corr,arr,names,trim=maxEdges)	
	# relevant_one,edges = toJson2(corr,arr,names)	
	# toJson2()
	# print(edges)

	# print(relevant_one,class_dic)
	radius,eccentricity,center,degreeDist,mod_score = analysis(edges,relevant_one,class_dic=class_dic)
	
	# print(radius,mod_score)
	# # sys.exit()


	rlt = {
		"edges":edges,
		"node_dic":relevant_one,
		# "degree_seq":degree_sequence,
		"degree_distribution":degreeDist,
		"radius":radius,
		"center":center,
		"eccentricity":eccentricity,
		"mod_score":round(mod_score,4),
		# 'missing_tickers':missing_tickers
		"debug":debug_msg
	}


	if hasClass and not ifClassGraph:
		rlt['symbol_class'] = class_dic
		rlt['classes'] = list(set(class_dic.values()))
		# print(class_dic)
		# print(relevant_one)

		class_cover = calculate_class_cover_ratio(edges,relevant_one,class_dic)

		rlt['class_correctness_ratio'] = class_cover

	print(json.dumps(rlt))



if __name__ == '__main__':
  main(sys.argv[1:])




