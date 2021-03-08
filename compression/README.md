# To make the compression part work:

1. Install `GEDLIB` by runing the `install.py` file from the gedlib folder. See the original gedlib **README.md** for details. You may need to use the `--lib gxl` flag.
2. Go to the compression folder and run the `setup.py` file from the terminal. This will create folders, copy some datasets and compress them using `tar`. This will be usefull to get the original collection sizes and the tar only size to compare the results.
3. Go to _gedlib/compression/build_.
4. Open a terminal and execute the following commands to compile the code:
```
cmake ..
make
```
5. Get one of the bash files from the *util_files* folder. I put them in the *build* folder, but you can do it from wherever you want, just update the paths in the file. These files are used to simplify execution. Define parameters, and execute the bash file using `bash file_name.sh`. There are different ways to run the code, and the tests were made using different execution files. I split the tests in different files because running them all together takes a lot of time.
Your results will be available in the _compression/data/output_ folder if you leave the parameters as they are.

The general structure of the compression folder is the following.

	compression
	|
	|-- setup.py			------------------- file to set up folders
	|-- CMakeLists.txt		------------------- for building GEDLIB with the compression part
	|-- create_directories.sh	------------------- helper file to create_directories bin, build and data
	|-- README.md			------------------- this file
	|-- bin
		|-- abc			------------------- the executable file to run the compression methods
	|
	|-- build
		|-- ...			------------------- files used for building. I run the run.sh file from here
	|-- data
	|   |--output			------------------- where to write results of tests. One folder for each dataset to be tested
	|	|-- acyclic 		------------------- folder for the output of the acyclic collection
	|	|-- AIDS
	|	|-- other datasets ...
	|	|-- orig_datasets_to_tar ------------------ folder to store copies of the collections and compress them using .tar.bz
	|-- ext				------------------- contains external implementations used in the project
	|	|-- atofigh
	|	|-- frangio68		------------------- used for the minimal spanning arborescence
	|	|-- stock_msts  	------------------- contains the code and the data used to create the msts collections
	|-- references			------------------- some of the academic material used during the project
	|	|-- ...
	|-- src
	|	|-- arborescence_based_compression.cpp ---- main file with the compression methods
	|	|-- arborescence_based_compression.hpp ---- header file
	|-- tests
	|	|-- compression_test.cpp ------------------ file containing the tests. Used for compressing, decompressing and getting results
	|	|-- README.txt		 ------------------ text file explaining how to replicate all the results, including comparing with tar.gz
	|-- util_files
		|-- check_iso.py	------------------- python file to test isomorphism between collections
		|-- copy_graphs.py	------------------- python file to copy graph gxl collections
		|-- clean_folders.sh	------------------- bash file used to clean the dataset folders in different situations
		|-- compress_tar.sh	------------------- bash file used to compress the original collections and get reference values
		|-- line_small.sh		------------------- example of a bash file used to excetude the script

