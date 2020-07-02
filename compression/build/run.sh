#From build dir
exec_path="src/"
exec_name="compression"
size="100"
girth="1"
# add_node, add_edge, rel_node, rel_edge, delete_node, delete_edge
prob="0 1 0 0 0 0"
keep_c="0"
stdout="0"
cost_version="0"
print_normal="0"
print_empty="0"
print_itself="0"
#collection="/home/lamsade/lgnecco/stage/gedlib/data/collections/Mutagenicity.xml"
collection="../../data/collections/Letter_A_Z.xml"
graph_dir="../../../data/datasets/Letter/MED"
output_file="../../data/output/results_compression.csv"
#graph_dir="/home/lamsade/lgnecco/stage/gedlib/compression/data/datasets/Mutagenicity/data"
#echo "Cmake"
#cmake ..
#echo "make"
#make
echo -ne '\007'
echo "execute"
#eval "./bin/Debug/stage_cpp" $size $girth $prob $keep_c $stdout
#eval "./"$exec_path$exec_name $prob $girth $collection $graph_dir $cost_version $print_normal $print_empty $print_itself
eval "./"$exec_path$exec_name $collection $graph_dir $stdout
echo -ne '\007'

