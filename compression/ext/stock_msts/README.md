# Installation

Please install required python libraries listed in ``requirements.txt``

# Create the folders with msts to use with GEDLIB

### 1. Use the stock informacion available at *./us_stock/archive/Stocks* to create a *csv* file with the time series in the right format

To do so, use the script `us_data_to_csv.py` available in the folder *./us_stock*. Here are two commands used to generate datasets.

	python3 us_data_to_csv.py --data_folder ./us_stock/archive/Stocks --start_date 1990-01-01 --max_stocks 3000 --min_values 7000 --out_file us_stocks_27y_7000.csv
	
	python3 us_data_to_csv.py --data_folder ./us_stock/archive/Stocks --start_date 1990-01-01 --max_stocks 4000 --min_values 4000 --out_file us_stocks_27y_4000.csv
	
This will generate a csv file in the current folder.
	
### 2. Use `mstEvolve.py` to generate the json file with the msts

Call the python script `mstEvolve.py` giving as input the csv file obtained in the previous step. Here are two examples.

	python3 mstEvolve.py --csv ./us_stocks_27y_7000.csv --step 3 --window 60 --edgeThreshold 0.1 --ifZNorm no --presenceTreshold 0.8 --overlapThreshold 0.8 --ifIncludeNegativeEdges false --out_file out_27y_7000.json
	
	python3 mstEvolve.py --csv ./us_stocks_27y_4000.csv --step 3 --window 60 --edgeThreshold 0.1 --ifZNorm no --presenceTreshold 0.8 --overlapThreshold 0.8 --ifIncludeNegativeEdges false --out_file out_27y_4000.json
	
This will generate a *json* file in the current folder called `out.json`. Change the name if you want to run the script again without overwriting it.

### 3. Generate the *xml* graphs from the *json* file.

Run the script `interpret_json.py` giving as input the *json* file returned in the previous step. Here are two examples assuming the *json* file names were changed.

	python3 interpret_json.py --json_file out_27y_7000.json --out_folder ./msts/
	python3 interpret_json.py --json_file out_27y_7000.json --out_folder ./msts/
	
The collections of *gxl* files representing the graphs should now be in the folder *./msts* in three folders, one for each type of edge weight. In *./msts* you should also find the corresponding *xml* files listing the collections
	

