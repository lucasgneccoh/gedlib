dirs=( "acyclic" "AIDS" "mao" "pah" "Mutagenicity" "msts_float_w" "msts_int_w" "msts_no_w")

eval "touch tar_times.txt"


for (( j = 1 ; j <= 1; j++ )) # In case we want to test many times
do
	echo "Compressing..."
	eval "rm -r compressed";
	eval "mkdir compressed";

	for dir in ${dirs[@]}; do
		echo $dir
		eval "echo "$dir" >> tar_times.txt";
		eval "{ time tar -cjf "$dir".tar.bz "$dir" ; } 2>> tar_times.txt";
		mv $dir".tar.bz" "compressed/"
	done

	eval "cd compressed"

	echo "Extracting..."
	eval "rm -r extract";
	eval "mkdir extract";

	for dir in ${dirs[@]}; do
		echo $dir
		eval "echo "$dir" >> tar_times.txt";
		eval "{ time tar -xjf "$dir".tar.bz -C extract ; } 2>> tar_times.txt";
	done
	
	eval "cd .."

done
eval "echo Done compressing and timing";
