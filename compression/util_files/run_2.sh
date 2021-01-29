#From build dir
exec_path="../bin/"
exec_name="abc"
# ------------------
# PARAMETERS. For arrays, use ":" as separator

stdout="0"

path_structure="true"
output_root="../data/output"
output_results_file="test_new_run.csv"
output_results_file_tar_size="test_new_run_tar_size.csv"
ged_method="branch_uniform"
ged_method_options="24"
ged_method_refinement="ipfp"
ged_method_refinement_options="24"
refinement_size="0"
write_ged_matrix="true"
write_arb="true"
write_results="true"
edit_cost_type="traditional"
relaxed_compression="true"

graph_sample_type="%"
binary_encoding="true"
decomp_only="false"
match_node_map="false"
match_node_map_by="stock"

# Define the graph sample here, let num trials in 1!
ks=("0" "10" "20" "30" "40" "50" "60" "70" "80" "90" "100")
ks=("10" "20")
datasets=("msts_no_w" "msts_int_w" "msts_float_w")
datasets=("acyclic" "mao")
for d in ${datasets[@]}; do
	for k in ${ks[@]}; do

		for (( j = 1 ; j <= 3; j++ ))
		do
			# in build
			eval $exec_path$exec_name "complete" $stdout $d "1" $path_structure $output_root $output_results_file $ged_method $ged_method_options $ged_method_refinement $ged_method_refinement_options $refinement_size $write_ged_matrix $write_arb $write_results $edit_cost_type $relaxed_compression $k $graph_sample_type $binary_encoding $decomp_only $match_node_map $match_node_map_by
		
			# move to dataset output folder
			eval "cd ../data/output/"$d
			eval "tar -cjf encoded_bin.tar.bz encoded_bin";
		
			# go back to build
			cd "../../../build"
			
			eval $exec_path$exec_name "tar_size" $stdout $d "1" $path_structure $output_root $output_results_file_tar_size $ged_method $ged_method_options $ged_method_refinement $ged_method_refinement_options $refinement_size $write_ged_matrix $write_arb $write_results $edit_cost_type $relaxed_compression $k $graph_sample_type $binary_encoding $decomp_only $match_node_map $match_node_map_by
		done
	done
done


