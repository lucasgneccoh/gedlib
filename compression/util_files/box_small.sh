#From build dir
exec_path="../bin/"
exec_name="abc"
# ------------------
# PARAMETERS. For arrays, use ":" as separator

# Test mode: "complete" for compression and decompression. See code for other modes
test_mode="complete"
decomp_only="false"

stdout="0"
datasets_names="acyclic:mao:pah"
num_trials="15"

# Paths for output
output_root="../data/output"
output_results_file="box_small.csv"

# GED computation
ged_method="branch_uniform"
ged_method_options="24"
ged_method_refinement="ipfp"
ged_method_refinement_options="24"
refinement_size="1"

# Handle outputs
write_ged_matrix="false"
write_arb="false"
write_results="true"

# Type of edit costs/encoding
edit_cost_type="trad"
relaxed_compression="true"
graph_sample_sizes="40"
graph_sample_type="%"
binary_encoding="true"

# Add ordered induced edges
path_structure="false"

# Compute node maps using a key attribute
match_node_map="false"
match_node_map_by="none"

write_headers="1"

# --- EXECUTION
echo -ne '\007'
eval $exec_path$exec_name $test_mode $stdout $datasets_names $num_trials $path_structure $output_root $output_results_file $ged_method $ged_method_options $ged_method_refinement $ged_method_refinement_options $refinement_size $write_ged_matrix $write_arb $write_results $edit_cost_type $relaxed_compression $graph_sample_sizes $graph_sample_type $binary_encoding $decomp_only $match_node_map $match_node_map_by $write_headers
echo -ne '\007'


