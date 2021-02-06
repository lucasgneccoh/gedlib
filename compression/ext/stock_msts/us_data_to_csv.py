import pandas as pd
import os
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from progress.bar import Bar
import timeit

def parseInputs():
	parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
	parser.add_argument("--data_folder", help="Folder with the individual txt files with the series")
	parser.add_argument("--max_stocks", help="Maximum number of series to consider")
	parser.add_argument("--out_file", help="path to the output file")
	parser.add_argument("--start_date", help="intial date to consider records. Format is yyyy-mm-dd")
	parser.add_argument("--end_date", help="final date to consider records. Format is yyyy-mm-dd")	
	parser.add_argument("--min_values", help="minimal number of non empty observations needed in the desired window to consider a series")	
	args = parser.parse_args()
	return args


# START

args = parseInputs()
num = 0
files = os.listdir(args.data_folder)
m = int(args.max_stocks) if not args.max_stocks is None else len(files)
bar = Bar("Joining", max = min(len(files), m))
start = timeit.default_timer()
for f in files:
	try:
		if len(f)>=4 and f[-4:]==".txt":
			# Read as pandas
			df = pd.read_csv(args.data_folder + "/" + f, header=0,
				index_col = "Date",
				usecols=["Date", "Close"])
			df.index = pd.to_datetime(df.index)
			if not args.start_date is None:
				df = df.loc[df.index >= pd.to_datetime(args.start_date), :]
			if not args.end_date is None:
				df = df.loc[df.index <= pd.to_datetime(args.end_date), :]

			if not args.min_values is None and df["Close"].notna().sum() < int(args.min_values):
				bar.next()
				continue
			
			df.rename(columns={"Close":f[0:-4]}, inplace=True)
			
			if num == 0:
				total = df.copy()
			else:
				total = total.join(df, how='outer', sort=True)
			

			if not args.max_stocks is None and num+1==int(args.max_stocks): 
				break
			num += 1
			bar.next()
	except:
		print("Error", f)
		bar.next()

bar.finish()
print()
print(timeit.default_timer() - start)

# Write data to a csv

import matplotlib.pyplot as plt
plt.plot(total.index, total.notna().sum(axis=1))
plt.show()

total.to_csv(args.out_file, na_rep = '\\N')

