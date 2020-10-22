# To make the codes work:
1. Install gedlib by runing the install.py file from the gedlib folder
2. Go to gedlib/compression/build and delete everyting but the run.sh file
3. Open a terminal and execute the following commands:
	cmake ..
	make
4. Now you can execute the compression code using the run.sh file. Change the variable values and run it in the same terminal 
	bash run.sh

NOTE: The run.sh file has a structure that depends on the code for now. As I have been doing different tests, you will find variables that are not used.

The general structure of the compression folder is the following. Use the create_directories.sh file to create the data folder that has many folders inside

_compression
|
|-- CMakeLists.txt		------------------------ for building GEDLIB with the compression part
|-- create_directories.sh	-------------------- helper file to create_directories bin, build and data
|-- README.txt			------------------------ this file
|-- _bin
|	|-- abc				------------------------ the executable file to run the compression methods
|
|-- _build
|	|-- ...				------------------------ files used for building
|-- _data
|   |--_output			------------------------ where to write results of tests. One folder for each dataset to be tested
|	|	|-- acyclic 	------------------------ folder for the output of the acyclic collection
|	|	|-- AIDS
|	|	|-- other datasets ...
|	|	|-- separate_files   ------------------ folder to store collections in separate files format, and their .tar.bz files
|	|	|-- separate_files_2 ------------------ folder to store collections in separate files (2) format, and their .tar.bz files
|	|-- orig_datasets_to_tar ------------------ folder to store copies of the collections and compress them using .tar.bz
|-- _ext				----------------------- contains external implementations used in the project
|	|-- atofigh
|	|-- frangio68	
|-- _references			----------------------- some of the academic material used during the project
|	|-- ...
|-- _src
|	|-- arborescence_based_compression.cpp ---- main file with the compression methods
|	|-- arborescence_based_compression.hpp ---- header file
|-- _tests
|	|-- compression_test.cpp	--------------- file containing the tests. Used for compressing and decompressing
|	|-- README.txt		----------------------- text file explaining how to replicate all the results, including comparing with tar.gz
|-- _util_files
	|-- check_iso.py	----------------------- python file to test isomorphism between collections
	|-- copy_graphs.py	----------------------- python file to copy graph gxl collections
	|-- clean_folders.sh	------------------- bash file used to clean the dataset folders in different situations
	|-- compress_tar.sh	----------------------- bash file used to compress the collections present in the directory




