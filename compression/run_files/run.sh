#From build dir
exec_path="src/"
exec_name="compression"

stdout="1"

suffix="lucas"

collection_file="/home/lucas/Documents/stage/gedlib/compression/data/test_collections/collections_list_"$suffix".txt"
graph_dir_file="/home/lucas/Documents/stage/gedlib/compression/data/test_collections/graphs_dir_list_"$suffix".txt"
dataset_file="/home/lucas/Documents/stage/gedlib/compression/data/test_collections/dataset_"$suffix".txt"


k_sample_file="/home/lucas/Documents/stage/gedlib/compression/data/test_collections/k_sample_file.txt"
output_root_file="/home/lucas/Documents/stage/gedlib/compression/data/test_collections/output_root.txt"


output_results_file="/home/lucas/Documents/stage/gedlib/compression/data/output/compressed/results_compression.csv"

ged_method="branch_uniform"
ged_method_options="24"
graph_sample_size="20"
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

edit_cost_type="second"
relaxed_coding="false"

echo -ne '\007'
echo "execute"
eval "./"$exec_path$exec_name $stdout $collection_file $graph_dir_file $output_root_file $dataset_file $output_results_file $ged_method $ged_method_options $graph_sample_size $ged_method_refinement $ged_method_refinement_options $refinement_size $encode $decode $write_decoded $train_set $train_path $ring_method $edit_cost_type $relaxed_coding $k_sample_file
echo -ne '\007'

