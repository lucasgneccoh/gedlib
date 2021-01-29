mkdir "bin"
mkdir "build"
mkdir "data"
cd "data"
dirs=("acyclic" "AIDS" "mao" "pah" "Mutagenicity" "Protein" "Letter" "msts_no_w" "msts_int_w" "msts_float_w" "compressed" "extract")
in_output=("decoded_text" "decoded_bin" "encoded_bin" "encoded_text" "extract")

mkdir "orig_datasets_to_tar"
cd "orig_datasets_to_tar"
for dir in ${dirs[@]}; do
	mkdir $dir
done

cd ".."
mkdir "output"
cd "output"


for dir in ${dirs[@]}; do
	mkdir $dir
	cd $dir
	for folder in ${in_output[@]}; do
		mkdir $folder
	done
	cd ".."
done

mkdir "separate_files"
cd "separate_files"
for dir in ${dirs[@]}; do
	mkdir $dir
done

cd ".."
mkdir "separate_files_2"
cd "separate_files_2"
for dir in ${dirs[@]}; do
	mkdir $dir
done

