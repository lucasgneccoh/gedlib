#From build dir
exec_path="src/"
exec_name="compression"

stdout="0"

collection_file="/home/lucas/Documents/stage/gedlib/compression/data/test_collections/collections_list.txt"
graph_dir_file="/home/lucas/Documents/stage/gedlib/compression/data/test_collections/graphs_dir_list.txt"
output_root_file="/home/lucas/Documents/stage/gedlib/compression/data/test_collections/output_root.txt"
dataset_file="/home/lucas/Documents/stage/gedlib/compression/data/test_collections/dataset.txt"

output_results_file="/home/lucas/Documents/stage/gedlib/compression/data/output/compressed/results_compression.csv"

ged_method="branch_uniform"
ged_method_options="24"
ged_method_refinement="ipfp"
ged_method_refinement_options="24"
refinement_size="0"

encode="true"
decode="true"
write_decoded="true"

# for the ring method
train_set="Letter-500g"
train_path="/home/lucas/Documents/stage/gedlib/compression/data/training"
ring_method="LSAPE_OPTIMAL"

echo -ne '\007'
echo "execute"
eval "./"$exec_path$exec_name $stdout $collection_file $graph_dir_file $output_root_file $dataset_file $output_results_file $ged_method $ged_method_options $ged_method_refinement $ged_method_refinement_options $refinement_size $encode $decode $write_decoded $train_set $train_path $ring_method
echo -ne '\007'

