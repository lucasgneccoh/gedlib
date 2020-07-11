#From build dir
exec_path="src/"
exec_name="compression"
collection="/home/lucas/Documents/stage/gedlib/compression/data/test_collections/collections_list_train_ring.txt"
graph_dir="/home/lucas/Documents/stage/gedlib/compression/data/test_collections/graphs_dir_list_train_ring.txt"
suffix="/home/lucas/Documents/stage/gedlib/compression/data/test_collections/suffix_train_ring.txt"
train_path="/home/lucas/Documents/stage/gedlib/compression/data/test_collections/train_path_train_ring.txt"
echo -ne '\007'
echo "execute"
eval "./"$exec_path$exec_name $collection $graph_dir $suffix $train_path
echo -ne '\007'

