dirs=( "acyclic" "AIDS" "mao" "pah" "Mutagenicity" "Protein" "Letter")

for dir in ${dirs[@]}; do
	echo $dir;
	eval "echo compress:";
	eval "time tar -cjf "$dir".tar.bz "$dir;
	mv $dir".tar.bz" "compressed/"
done

eval "cd compressed"

eval "rm -r extract";
eval "mkdir extract";

for dir in ${dirs[@]}; do
	echo $dir;
	eval "time tar -xjf "$dir".tar.bz -C extract";
	
	eval "echo --------------------------";
done
