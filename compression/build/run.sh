#From build dir
exec_path="src/"
exec_name="test_abc"

suffix="lucas"

# Parameters

stdout="1"

collection_file="/home/lucas/Documents/stage/gedlib/compression/data/test_collections/collections_list_"$suffix".txt"
graph_dir="/home/lucas/Documents/stage/gedlib/compression/data/test_collections/graphs_dir_list_"$suffix".txt"
output_root="/home/lucas/Documents/stage/gedlib/compression/data/test_collections/output_root.txt"
file_preffix="/home/lucas/Documents/stage/gedlib/compression/data/test_collections/dataset_"$suffix".txt"

k_sample_file="/home/lucas/Documents/stage/gedlib/compression/data/test_collections/k_sample_file.txt"



output_results_file="/home/lucas/Documents/stage/gedlib/compression/data/output/compressed/results_compression.csv"

ged_method="branch_uniform"
ged_method_options="24"
graph_sample_size="20"
ged_method_refinement="branch_fast"
ged_method_refinement_options="24"
refinement_size="999"

write_ged_matrix="true"
write_arb="true"
write_results="true"

# for the ring method
train_set="Letter-500g"
train_path="/home/lucas/Documents/stage/gedlib/compression/data/training"
ring_method="LSAPE_OPTIMAL"

edit_cost_type="mod"
relaxed_coding="false"

echo -ne '\007'
echo "execute"
eval "./"$exec_path$exec_name $stdout $collection_file $graph_dir $output_root $file_preffix $output_results_file $ged_method $ged_method_options $graph_sample_size $ged_method_refinement $ged_method_refinement_options $refinement_size $write_ged_matrix $write_arb $write_results $train_set $train_path $ring_method $edit_cost_type $relaxed_coding $k_sample_file
echo -ne '\007'

