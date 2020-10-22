#From build dir
# python3 install.py --tests abc_compression --clean --lib gxl
exec_path="compression/bin/"
exec_name="abc"
root_dir="/home/lucas/Documents"
#root_dir="/home/lamsade/lgnecco"

suffix="pah"
num_trials="1"
test_mode="false"

graph_sample_type="%"

# Parameters

stdout="1"

collection_file=$root_dir"/stage/gedlib/compression/data/test_collections/collections_list_"$suffix".txt"
graph_dir=$root_dir"/stage/gedlib/compression/data/test_collections/graphs_dir_list_"$suffix".txt"
output_root=$root_dir"/stage/gedlib/compression/data/test_collections/output_root.txt"
file_preffix=$root_dir"/stage/gedlib/compression/data/test_collections/dataset_"$suffix".txt"

k_sample_file=$root_dir"/stage/gedlib/compression/data/test_collections/k_sample_file.txt"

output_results_file=$root_dir"/stage/gedlib/compression/data/output/compressed/test_results_new_format.csv"

ged_method="branch_uniform"
ged_method_options="24"
ged_method_refinement="ipfp"
ged_method_refinement_options="24"
refinement_size="999"

write_ged_matrix="f"
write_arb="f"
write_results="true"


edit_cost_type="trad"
relaxed_coding="true"

echo -ne '\007'
echo "execute"
eval "./"$exec_path$exec_name $stdout $collection_file $graph_dir $output_root $file_preffix $output_results_file $ged_method $ged_method_options $ged_method_refinement $ged_method_refinement_options $refinement_size $write_ged_matrix $write_arb $write_results $edit_cost_type $relaxed_coding $k_sample_file $num_trials $test_mode $graph_sample_type
echo -ne '\007'

