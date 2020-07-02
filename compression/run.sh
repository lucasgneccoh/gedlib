#From build dir
exec_path="src/"
exec_name="compression"
collection="/home/lucas/Documents/stage/gedlib/compression/data/test_collections/collections_list_muta.txt"
graph_dir="/home/lucas/Documents/stage/gedlib/compression/data/test_collections/graphs_dir_list_muta.txt"
output_file="/home/lucas/Documents/stage/gedlib/compression/data/output/results_compression_muta.csv"
stdout="true"
echo -ne '\007'
echo "execute"
eval "./"$exec_path$exec_name $collection $graph_dir $output_file $stdout
echo -ne '\007'

