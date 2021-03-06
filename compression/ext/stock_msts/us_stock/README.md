To create the csv from the US Stocks data. Data taken from 
https://www.kaggle.com/borismarjanovic/price-volume-data-for-all-us-stocks-etfs

# run

To create a *csv* file from these raw text files, use the script `us_data_to_csv.py`. The parameters are:

data_folder: Folder with the separate txt files
max_stocks: Max number of time series to consider. Will limit the number of columns of the csv file and will affect the number of edges in the msts
out_file: Path to the output file
start_date: Initial date to consider. All observations prior to this date are discarded. Format yyyy-mm-dd
end_date: Final date. All observations after this date are discarded. Format yyyy-mm-dd
min_values: Minimal number of observations in the given window of time to consider a time series. If the series has less values, it is discarded. 

#### Example to test

	python3 us_data_to_csv.py --data_folder archive/Stocks --max_stocks 50 --out_file us_stocks_test.csv --start_date 2015-01-01 --end_date 2015-12-31 --min_values 150


Here are the commands used to generate the two possible collections:

	python3 us_data_to_csv.py --data_folder archive/Stocks --start_date 1990-01-01 --max_stocks 3000 --min_values 7000 --out_file us_stocks_27y_7000.csv

	python3 us_data_to_csv.py --data_folder archive/Stocks --start_date 1990-01-01 --max_stocks 4000 --min_values 4000 --out_file us_stocks_27y_4000.csv


