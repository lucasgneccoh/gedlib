import shutil
from os import listdir
from os.path import isfile, join

# run files
mypath = "run_files"
build_path = "build"
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
	shutil.copy(path, build_path)

print("Done with run files")


# files with test cases

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
