for dir in */;
	do echo $dir;
	eval "rm -r ./"$dir"*"
done
