The tests are divided into two main parts:

 - Part 1: Compression and decompression changing the k parameter (sample size) that determines the "density" of the graph
 - Part 2: Comparing the method with tar.bz (or tar.gz) on the original files or the "separated files" format.
	
	
*** Part 1 ***

For part one the only thing you have to do is to run the compression_test.cpp file with the desired settings.

We used 5 iterations for each step, and k took the values 10%, 20%, 30%, ..., 100%. This was made using binaru encoding, relaxed encoding, traditional costs, branch_uniform as main method, ipfp as refinement method, and a refinement size of 1.

The results will be stored in a csv file.


*** Part 2 ***

Part 2 requires a little more work. First you need to manually organize two folders. One should have the original datasets (xml file + gxl files) and the other one the "separated files" format. We will call them folder_gxl_orig and folder_separate.

Important: To be able to compare the sizes of the collections in different formats, each representation should have the same information. From any representation one should be able to recreate the same collection with the same information. 

In my current folder organization, these two folders are located in 

	- folder_gxl_orig -> gedlib/compression/data/orig_datasets_to_tar
	- folder_separate -> gedlib/compression/data/output/separate_files (WITH attribute dictionaries in separate file)
	- folder_separate -> gedlib/compression/data/output/separate_files_2 (WITHOUT attribute dictionaries in separate file)
	
In each folder there should be one folder for each dataset that is going to be used.

	folder_gxl_orig: To create a copy of the original collection, a python file is distributed. The idea is that each dataset should have a folder containing the gxl files and the xml file listing the graphs. Pretty simple.
	
	folder_separate: Again, one folder per dataset. This folder will contain the files describing the collection. 
	
		dataset.edges for all the edges
		dataset.graph_idx to relate nodes and graphs
		dataset.graph_labels with the graph labels or classes
		dataset.node_labels with the node labels or attributes
		dataset.edge_labels with the edge labels or attributes

		If the attributes are encoded and a dictionary is being used, then 
		dataset.node_attr with the node attribute dictionary
		dataset.edge_attr with the edge attribute dictionary		
		
		To create this files, there is a function in the compress_test.cpp file. You just need to specify the datasets and the paths where to save the information.


After creating the folders and getting all the files in each representation, the next step is to compress the collections for each dataset. This process should be timed, to get an idea of the time taken in compression and decompression. This also needs to be done in the output folder where the datasets where compressed using the ABC method (with the folders encoded_bin and encoded_text).

Once you have all the files AND all the compressed files (.tar.bz), there is also a function in compress_test.cpp that will go and look for all the file and folder sizes to make a table comparing all the results. This will provide a comparisson between the abc method, abc method + tar.bz and tar.bz only.

**** So remember: You should have three (actually four) representations of each dataset: Original gxl files, separate files, and ABC compressed (in binary and text). And for each representation, you should also have a .tar.gz file compressing the representation.


I distribute some files to help with the tasks of copying the graphs, compressing, extracting, etc:

	- copy_graphs.py: only in folder_gxl_orig. It is used to create a copy of the gxl files from the main data source to a new folder
	
	- clean_folders.sh: This file will help you reset the folders by erasing everything and creating the empty folders with the dataset names.
	
	- compress_tar.sh: Depending on the folder the behaviour may change a bit, but the idea is that this file will help you compress the information of each dataset and store the .tar.bz (or gz) file somewhere to then extract it somewhere else. Each operation will be timed and results shown on the screen
	

