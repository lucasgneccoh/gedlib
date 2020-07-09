#From build dir
exec_path="src/"
exec_name="compression"
collection="/home/lucas/Documents/stage/gedlib/compression/data/test_collections/collections_list_letters.txt"
graph_dir="/home/lucas/Documents/stage/gedlib/compression/data/test_collections/graphs_dir_list_letters.txt"
output_file="/home/lucas/Documents/stage/gedlib/compression/data/output/results_compression_test_letters.csv"
stdout="true"
ged_method="branch_fast"
ged_method_options="10"
ged_method_refinement="branch_uniform"
ged_method_refinement_options="10"
refinement_size="999"
train_set="/home/lucas/Documents/stage/gedlib/compression/data/collections/Letter-50g.xml"
train_path="/home/lucas/Documents/stage/gedlib/compression/data/training"
output_description_file="/home/lucas/Documents/stage/gedlib/compression/data/output"
output_file_name="alphabets.data"
echo -ne '\007'
echo "execute"
eval "./"$exec_path$exec_name $collection $graph_dir $output_file $stdout $ged_method $ged_method_options $ged_method_refinement $ged_method_refinement_options $refinement_size $train_set $train_path $output_description_file $output_file_name
echo -ne '\007'

