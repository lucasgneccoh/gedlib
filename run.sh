#From build dir
# python3 install.py --tests abc_compression --clean --lib gxl
exec_path="compression/bin/"
exec_name="abc"
root_dir="/home/lucas/Documents"
#root_dir="/home/lamsade/lgnecco"

suffix="pah"
num_trials="1"
test_mode="false"

# Parameters

stdout="1"

collection_file=$root_dir"/stage/gedlib/compression/data/test_collections/collections_list_"$suffix".txt"
graph_dir=$root_dir"/stage/gedlib/compression/data/test_collections/graphs_dir_list_"$suffix".txt"
output_root=$root_dir"/stage/gedlib/compression/data/test_collections/output_root.txt"
file_preffix=$root_dir"/stage/gedlib/compression/data/test_collections/dataset_"$suffix".txt"

k_sample_file=$root_dir"/stage/gedlib/compression/data/test_collections/k_sample_file.txt"

output_results_file=$root_dir"/stage/gedlib/compression/data/output/compressed/results_refinement_2.csv"

ged_method="branch_fast"
ged_method_options="24"
ged_method_refinement="ipfp"
ged_method_refinement_options="24"
refinement_size="999999"

write_ged_matrix="true"
write_arb="true"
write_results="true"


edit_cost_type="trad"
relaxed_coding="true"

echo -ne '\007'
echo "execute"
eval "./"$exec_path$exec_name $stdout $collection_file $graph_dir $output_root $file_preffix $output_results_file $ged_method $ged_method_options $ged_method_refinement $ged_method_refinement_options $refinement_size $write_ged_matrix $write_arb $write_results $edit_cost_type $relaxed_coding $k_sample_file $num_trials $test_mode
echo -ne '\007'

