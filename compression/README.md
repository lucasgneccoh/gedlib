## To make it work:
1. Install gedlib by runing the install.py file from the gedlib folder
2. Go to the compression folder and run the create_directories.sh file. If you dont want or cant run it, you will have to create all the folders yourself to avoid problems with the output of the compression method.
3. Go to gedlib/compression/build
4. Open a terminal and execute the following commands:
```
cmake ..
make
```
	
5. Get the run.sh file from the util_files folder. This file is used to simplify execution. Define parameters, and execute the bash file using "bash run.sh".
Yout results will be available in the *compression/data/output* folder. I you use `run_2.sh` then you will have two files.


The general structure of the compression folder is the following. Use the create_directories.sh file to create the data folder that has many folders inside
```
compression
|
|-- CMakeLists.txt		------------------------ for building GEDLIB with the compression part
|-- create_directories.sh	-------------------- helper file to create_directories bin, build and data
|-- README.txt			------------------------ this file
|-- bin
	|-- abc				------------------------ the executable file to run the compression methods
|
|-- build
	|-- ...				------------------------ files used for building
|-- data
|   |--output			------------------------ where to write results of tests. One folder for each dataset to be tested
|	|-- acyclic 	------------------------ folder for the output of the acyclic collection
|	|-- AIDS
|	|-- other datasets ...
|	|-- separate_files   ------------------ folder to store collections in separate files format, and their .tar.bz files
|	|-- separate_files_2 ------------------ folder to store collections in separate files (2) format, and their .tar.bz files
|	|-- orig_datasets_to_tar ---------------------- folder to store copies of the collections and compress them using .tar.bz
|-- ext				----------------------- contains external implementations used in the project
|	|-- atofigh
|	|-- frangio68	
|-- references			----------------------- some of the academic material used during the project
|	|-- ...
|-- src
|	|-- arborescence_based_compression.cpp ---- main file with the compression methods
|	|-- arborescence_based_compression.hpp ---- header file
|-- tests
|	|-- compression_test.cpp	--------------- file containing the tests. Used for compressing and decompressing adn getting results
|	|-- README.txt		----------------------- text file explaining how to replicate all the results, including comparing with tar.gz
|-- util_files
	|-- check_iso.py	----------------------- python file to test isomorphism between collections
	|-- copy_graphs.py	----------------------- python file to copy graph gxl collections
	|-- clean_folders.sh	----------------------- bash file used to clean the dataset folders in different situations
	|-- compress_tar.sh	----------------------- bash file used to compress the
	|-- run.sh	------------------------------- bash file used to execute
	|-- run_2.sh	--------------------------- bash file used to execute and avoid system call problems with the tar call
```


