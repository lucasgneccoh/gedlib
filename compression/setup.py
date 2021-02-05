import os
import subprocess
import shutil


# Create folders

os.makedirs("bin", exist_ok=True)
os.makedirs("build", exist_ok=True)
os.makedirs("data", exist_ok=True)

dirs="acyclic:AIDS:mao:pah:Mutagenicity:msts_no_w:msts_int_w:msts_float_w"
tar_folders="compressed:extract"
in_output="decoded_text:decoded_bin:encoded_bin:encoded_text:extract"

cwd = "data/orig_datasets_to_tar"

for d in dirs.split(":"):
	os.makedirs(cwd+'/'+d, exist_ok=True)

for d in tar_folders.split(":"):
	os.makedirs(cwd+'/'+d, exist_ok=True)

cwd = "data/orig_datasets_to_tar/separate_files_attr"

for d in dirs.split(":"):
	os.makedirs(cwd+'/'+d, exist_ok=True)

for d in tar_folders.split(":"):
	os.makedirs(cwd+'/'+d, exist_ok=True)

cwd = "data/orig_datasets_to_tar/separate_files_no_attr"

for d in dirs.split(":"):
	os.makedirs(cwd+'/'+d, exist_ok=True)

for d in tar_folders.split(":"):
	os.makedirs(cwd+'/'+d, exist_ok=True)

cwd = "data/output"
for d in dirs.split(":"):
	for folder in in_output.split(":"):
		os.makedirs(cwd+'/'+d+'/'+folder, exist_ok=True)


# Copy original datasets
py = None
ans = subprocess.call("python3 --version", shell=True)
if ans==0:
	py = "python3"
else:
	ans = subprocess.call("python --version", shell=True)
	if ans==0:
		py = "python"
	else:
		ans = subprocess.call("python2 --version", shell=True)
		if ans==0:
			py = "python2"		
	
if py is None:
	print("could not find the command to call python. Tried python, python2 and python3")
else:
	cmd = py + ' util_files/copy_graphs.py --target data/orig_datasets_to_tar --base_xml ../data/collections/ --base_datasets ../data/datasets/ --datasets ' + ' '.join(dirs.split(":"))
	ans = subprocess.call(cmd, shell=True)
	print("Original files copied")
	print("Launching compression")
	print("This may take some time (2 min max)")
	shutil.copy("util_files/compress_tar.sh", "data/orig_datasets_to_tar")
	os.chdir("data/orig_datasets_to_tar")
	cmd = "bash compress_tar.sh"
	ans = subprocess.call(cmd, shell=True)
	if ans==0:
		print("Original collections copied and compressed. See data/orig_datasets_to_tar/tar_times.txt for the times taken for compression and decompression using tar.bz")
	else:
		print("Something went wrong when compressing the original collections") 
	




