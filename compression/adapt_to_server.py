import shutil

# run.sh
with open("run.sh", "r") as source:
	with open("build/run.sh", "w") as dest:
		line = source.readline()
		while(line):
			line2 = line.replace("/home/lucas/Documents", "/home/lamsade/lgnecco")
			dest.write(line2)
			line = source.readline()

print("Done with run.sh")
# files with test cases
from os import listdir
from os.path import isfile, join
mypath = "data/test_collections"
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
for filename in onlyfiles:
	path = join(mypath, filename)
	temp = path + ".tmp"
	with open(path, "r") as source:
		with open(temp, "w") as dest:
			line = source.readline()
			while(line):
				line2 = line.replace("/home/lucas/Documents", "/home/lamsade/lgnecco")
				dest.write(line2)
				line = source.readline()
				
	shutil.move(temp, path)
	
print("Done with collection sets")
