#From build dir
exec_path="../bin/"
exec_name="abc"
# ------------------
# PARAMETERS. For arrays, use ":" as separator
test_mode="complete"
stdout="0"
datasets_names="acyclic:mao:pah"
num_trials="30"
path_structure="false"
output_root="../data/output"
output_results_file="line_small.csv"
ged_method="branch_uniform"
ged_method_options="24"
ged_method_refinement="ipfp"
ged_method_refinement_options="24"
refinement_size="1"
write_ged_matrix="false"
write_arb="false"
write_results="true"
edit_cost_type="trad"
relaxed_compression="true"
graph_sample_sizes="0:10:20:30:40:50:60:70:80:90:100"
graph_sample_type="%"
binary_encoding="true"
decomp_only="false"
match_node_map="false"
match_node_map_by="stock"
write_headers="1"

# --- EXECUTION
echo -ne '\007'
eval $exec_path$exec_name $test_mode $stdout $datasets_names $num_trials $path_structure $output_root $output_results_file $ged_method $ged_method_options $ged_method_refinement $ged_method_refinement_options $refinement_size $write_ged_matrix $write_arb $write_results $edit_cost_type $relaxed_compression $graph_sample_sizes $graph_sample_type $binary_encoding $decomp_only $match_node_map $match_node_map_by $write_headers
echo -ne '\007'


