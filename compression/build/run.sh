#From build dir
exec_path="src/"
exec_name="compression"
size="100"
girth="1"
# add_node, add_edge, rel_node, rel_edge, delete_node, delete_edge
prob="0 0 0 0 0 1"
keep_c="0"
stdout="2"
collection="/home/lucas/Documents/stage/gedlib/data/collections/Letter.xml"
graph_dir="/home/lucas/Documents/stage/gedlib/data/datasets/Letter/HIGH"
#echo "Cmake"
#cmake ..
#echo "make"
#make
echo -ne '\007'
echo "execute"
#eval "./bin/Debug/stage_cpp" $size $girth $prob $keep_c $stdout
eval "./"$exec_path$exec_name $prob $girth $collection $graph_dir
echo -ne '\007'

