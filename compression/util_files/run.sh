#From build dir
exec_path="../bin/"
exec_name="abc"
# ------------------
# PARAMETERS. For arrays, use ":" as separator
test_mode="complete"
stdout="1"
datasets_names="msts_no_w:msts_int_w:msts_float_w:msts_no_w_un:msts_int_w_un:msts_float_w_un"
num_trials="1"
path_structure="true"
output_root="../data/output"
output_results_file="testing_msts.csv"
ged_method="branch_uniform"
ged_method_options="24"
ged_method_refinement="ipfp"
ged_method_refinement_options="24"
refinement_size="0"
write_ged_matrix="true"
write_arb="true"
write_results="true"
edit_cost_type="mod"
relaxed_compression="true"
graph_sample_sizes="100"
graph_sample_type="%"
binary_encoding="true"
decomp_only="false"
match_node_map="true"
match_node_map_by="stock"

# --- EXECUTION
echo -ne '\007'
eval $exec_path$exec_name $test_mode $stdout $datasets_names $num_trials $path_structure $output_root $output_results_file $ged_method $ged_method_options $ged_method_refinement $ged_method_refinement_options $refinement_size $write_ged_matrix $write_arb $write_results $edit_cost_type $relaxed_compression $graph_sample_sizes $graph_sample_type $binary_encoding $decomp_only $match_node_map $match_node_map_by
echo -ne '\007'

