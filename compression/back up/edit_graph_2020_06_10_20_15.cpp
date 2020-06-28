#include <iostream>
#include "/home/lucas/Documents/stage_gedlibpy/gedlib/gedlib/src/env/ged_env.hpp"
#include "/home/lucas/Documents/stage_gedlibpy/gedlib/gedlib/median/src/median_graph_estimator.hpp"
#include <random>
#include <string>

#include "edit_graph.h"

void describe_graph(ged::ExchangeGraph<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> exchange_graph){
	std::cout<<"NODES:\n";
	for (std::size_t i{0}; i < exchange_graph.num_nodes; i++) {
		// TODO: for sobre los labels
		std::cout << "Node " << exchange_graph.original_node_ids.at(i) <<":\n";
		for(auto const &l : exchange_graph.node_labels.at(i)){
			std::cout<<"\t"<<l.first<<": " <<l.second << "\n";
		}
		
	}
	std::cout<<"EDGES:\n";
	for (std::size_t i{0}; i < exchange_graph.num_nodes; i++) {
		for (std::size_t j{i + 1}; j < exchange_graph.num_nodes; j++) {
			if (exchange_graph.adj_matrix[i][j] == 1) {
				std::cout << "Edge:  " << exchange_graph.original_node_ids.at(i);
				std::cout << " -- " << exchange_graph.original_node_ids.at(j) << "\n";
				for(auto const &e : exchange_graph.edge_labels.at(std::make_pair(i, j))){
					std::cout<<"\t"<<e.first<<": " <<e.second << "\n";
				} 
			}
		}
	}
}

//TODO: Add templates here
std::map<std::string, std::map<std::string, std::vector<std::string>>> 
	get_graphs_structure(ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> &env){

	std::map<std::string, std::vector<std::string>> node_attr;
	std::map<std::string, std::vector<std::string>> edge_attr;

	ged::ExchangeGraph<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> g;
	std::pair<ged::GEDGraph::GraphID, ged::GEDGraph::GraphID> graph_ids = env.graph_ids();
	for (ged::GEDGraph::GraphID g_id{graph_ids.first}; g_id<graph_ids.second; g_id++){
		g = env.get_graph(g_id);
		for (std::size_t i{0}; i < g.num_nodes; i++) {
			for(auto l : g.node_labels.at(i)){
				if(node_attr.count(l.first)==0){
					node_attr.emplace(std::make_pair(l.first, std::vector<std::string>()));
				}		
				node_attr.at(l.first).emplace_back(l.second);
			}
		}
		for (std::size_t i{0}; i < g.num_nodes; i++) {
			for (std::size_t j{i + 1}; j < g.num_nodes; j++) {
				if (g.adj_matrix[i][j] == 1) {				
					for(auto e : g.edge_labels.at(std::make_pair(i, j))){
						if(edge_attr.count(e.first)==0){
							edge_attr.emplace(std::make_pair(e.first, std::vector<std::string>()));
						}
						edge_attr.at(e.first).emplace_back(e.second);
					} 
				}
			}
		}
	}
	std::map<std::string, std::map<std::string, std::vector<std::string>>> result;
	result.emplace(std::make_pair("node_attr", node_attr));
	result.emplace(std::make_pair("edge_attr", edge_attr));	
	return result;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
edit_add_node(ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &exchange_graph,
		UserNodeID node_id, std::map<std::string, std::map<std::string, std::vector<std::string>>> structure){
	
	UserNodeLabel label;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<int> distr(0, 1);

	std::cout<<".";
	// Create the node label according to the graph structure
	for(auto const attr : structure.at("node_attr")){
		distr = std::uniform_int_distribution<int>(0, attr.second.size()-1);		
		label.emplace(std::make_pair(attr.first, attr.second.at(distr(gen))));
	}

	std::cout<<".";
	exchange_graph.original_node_ids.emplace_back(node_id);
	
	std::cout<<".";
	exchange_graph.node_labels.emplace_back(label);

	std::cout<<".";
	std::vector<std::size_t> last_line;
	for (std::size_t i{0}; i < exchange_graph.num_nodes; i++){
		exchange_graph.adj_matrix.at(i).emplace_back(0);
		last_line.emplace_back(0);
	}

	std::cout<<".";
	last_line.emplace_back(0);

	std::cout<<".";
	exchange_graph.adj_matrix.emplace_back(last_line);

	//env.add_node(g_id, node_id, label);
	std::cout<<".";
	exchange_graph.num_nodes++;
	return;

}


template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
edit_add_edge(ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &exchange_graph,
		std::map<std::string, std::map<std::string, std::vector<std::string>>> structure, bool ignore_duplicates){
	// We create the UserNodeLabel (GXLLabel in our case) and use the existing env.add_node method
	if(exchange_graph.num_nodes<2){
		return;
	}
	std::cout<<".";
	UserEdgeLabel label;
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<std::size_t> distr(0, 1);

	std::cout<<".";
	// Create the edge label according to the graph structure
	for(auto const attr : structure.at("edge_attr")){
		distr = std::uniform_int_distribution<std::size_t>(0, attr.second.size()-1);
		label.emplace(std::make_pair(attr.first, attr.second.at(distr(gen))));
	}

	std::cout<<".";
	std::pair<std::size_t, std::size_t> e_add;
	if(!ignore_duplicates){
		std::vector<std::pair<std::size_t, std::size_t>> non_edges;
		std::cout<<".";
		// Get non edges
		for (std::size_t i{0}; i < exchange_graph.num_nodes; i++) {
			for (std::size_t j{i + 1}; j < exchange_graph.num_nodes; j++) {
				std::cout<<".";
				if (exchange_graph.adj_matrix[i][j] == 0) {	
					non_edges.emplace_back(std::make_pair(i,j));
				}
			}
		}
		std::cout<<".";
		// Chose the non edge to add
		distr = std::uniform_int_distribution<std::size_t>(0, non_edges.size()-1);
		e_add = non_edges.at(distr(gen));
	}
	else{
		std::cout<<".";
		// Does not matter if there are multi-edges. We pick two random nodes
		distr = std::uniform_int_distribution<std::size_t>(0, exchange_graph.num_nodes-1);
		std::size_t i = distr(gen);
		std::size_t j = distr(gen);
		// For now, avoid loops
		while(i != j){
			std::cout<<".";
			std::size_t j = distr(gen);
		}
		std::cout<<".";
		e_add = std::make_pair(i,j);
	}

	//UserNodeID from_id = exchange_graph.original_node_ids.at(e_add.first);
	//UserNodeID to_id = exchange_graph.original_node_ids.at(e_add.second);

	//env.add_edge(g_id, from_id, to_id, label, ignore_duplicates);

	// Modify adj_matrix and edge_list
	std::cout<<".";
	exchange_graph.adj_matrix.at(e_add.first).at(e_add.second) = 1;
	std::cout<<".";
	exchange_graph.edge_list.emplace_back(std::make_pair(e_add, label));
	std::cout<<".";
	exchange_graph.edge_labels.emplace(std::make_pair(e_add, label));
	std::cout<<".";
	exchange_graph.num_edges++;
	return;

}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
edit_transform_node(ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &exchange_graph,
		std::map<std::string, std::map<std::string, std::vector<std::string>>> structure){

	std::cout<<".";
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<std::size_t> distr(0, exchange_graph.num_nodes-1);
	std::size_t sel_node = distr(gen);

	std::cout<<".";
	UserNodeLabel label = exchange_graph.node_labels.at(sel_node);
	std::cout<<".";
	// Transform label according to the structure of the graphs. Only modify existing attributes
	for(auto  attr : exchange_graph.node_labels.at(sel_node)){
		std::cout<<".";
		distr = std::uniform_int_distribution<std::size_t>(0, structure.at("node_attr").at(attr.first).size()-1);		
		exchange_graph.node_labels.at(sel_node).at(attr.first) = structure.at("node_attr").at(attr.first).at(distr(gen));
	}

	// Update graph in env. Keep same name and class
	//env.load_exchange_graph(exchange_graph, g_id, ged::Options::ExchangeGraphType::ADJ_MATRIX ,env.get_graph_name(g_id) , env.get_graph_class(g_id));
	return;
}


template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
edit_transform_edge(ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &exchange_graph,
		std::map<std::string, std::map<std::string, std::vector<std::string>>> structure){

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<std::size_t> distr;
	distr = std::uniform_int_distribution<std::size_t>(0, exchange_graph.num_edges-1);
	std::size_t sel_edge = distr(gen);
	std::size_t cont{0};
	std::cout<<".";
	typename std::list<std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel>>::iterator iter;
	std::cout<<".";
	iter = exchange_graph.edge_list.begin();
	std::cout<<".";
	advance(iter, sel_edge);
	std::cout<<".";
	std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel> edge = *iter;
	std::cout<<".";
	// Transform label according to the structure of the graphs. Only modify existing attributes
	for(auto attr : edge.second){
		std::cout<<".";
		distr = std::uniform_int_distribution<std::size_t>(0, structure.at("edge_attr").at(attr.first).size()-1);		
		edge.second.at(attr.first) = structure.at("edge_attr").at(attr.first).at(distr(gen));
		exchange_graph.edge_labels.at(edge.first) = edge.second;
	}

	//env.load_exchange_graph(exchange_graph, g_id, ged::Options::ExchangeGraphType::EDGE_LIST ,env.get_graph_name(g_id) , env.get_graph_class(g_id));
	
}


template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
edit_remove_node(ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &exchange_graph){
	
	if(exchange_graph.num_nodes==0){
		return;
	}
	// Remove even if it has edges??
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<std::size_t> distr(0, exchange_graph.num_nodes-1);
	std::size_t n = distr(gen);
	std::cout<<".";
	exchange_graph.original_node_ids.erase(exchange_graph.original_node_ids.begin() + n);
	std::cout<<".";
	exchange_graph.node_labels.erase(exchange_graph.node_labels.begin() + n);

	std::cout<<".";
	// TODO: May cause problems if vectors become empty??
	typename std::vector<std::size_t>::iterator iter_row;
	for (ged::GEDGraph::NodeID i{0}; i < exchange_graph.num_nodes; i++){
		iter_row = exchange_graph.adj_matrix.at(i).begin();
		advance(iter_row, n);
		exchange_graph.adj_matrix.at(i).erase(iter_row);	
	}
	std::cout<<".";
	typename std::vector<std::vector<std::size_t>>::iterator iter_col;
	std::cout<<".";
	iter_col = exchange_graph.adj_matrix.begin();
	std::cout<<".";
	advance(iter_col, n);
	std::cout<<".";
	exchange_graph.adj_matrix.erase(iter_col);
	std::cout<<".";
	typename std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel> entry;
	for(typename std::list<std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel>>::iterator iter = exchange_graph.edge_list.begin();
		iter != exchange_graph.edge_list.end(); iter++){
		std::cout<<".";
		entry = *iter;
		if(entry.first.first == n | entry.first.second == n){			
			exchange_graph.edge_labels.erase(entry.first);
			iter = exchange_graph.edge_list.erase(iter);
			iter--;
		}
	}
	std::cout<<".";
	// re-index nodes in edge list
	typename std::map<std::pair<std::size_t, std::size_t>, UserEdgeLabel> new_edge_labels;
	for(typename std::list<std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel>>::iterator iter = exchange_graph.edge_list.begin();
		iter != exchange_graph.edge_list.end(); iter++){
		std::cout<<".";
		entry = *iter;		
		if(entry.first.first > n){(*iter).first.first--;}
		if(entry.first.second > n){(*iter).first.second--;}
		new_edge_labels.emplace(entry);
	}

	exchange_graph.edge_labels = new_edge_labels;
	exchange_graph.num_nodes--;
	exchange_graph.num_edges = exchange_graph.edge_list.size();
	//env.load_exchange_graph(exchange_graph, g_id, ged::Options::ExchangeGraphType::ADJ_MATRIX ,env.get_graph_name(g_id) , env.get_graph_class(g_id));
	
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
edit_remove_edge(ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &exchange_graph){
	
	if(exchange_graph.num_edges==0){
		return;
	}
	std::cout<<".";
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<std::size_t> distr(0, exchange_graph.num_edges-1);
	std::size_t n = distr(gen);	
	std::cout<<".";
	typename std::list<std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel>>::iterator iter;
	typename std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel> edge;
	iter = exchange_graph.edge_list.begin();
	advance(iter,n);
	edge = *iter;
	std::cout<<".";
	exchange_graph.adj_matrix.at(edge.first.first).at(edge.first.second)=0;
	std::cout<<".";
	exchange_graph.edge_labels.erase(edge.first);
	std::cout<<".";
	exchange_graph.edge_list.erase(iter);
	std::cout<<".";
	exchange_graph.num_edges--;
	//env.load_exchange_graph(exchange_graph, g_id, ged::Options::ExchangeGraphType::ADJ_MATRIX ,env.get_graph_name(g_id) , env.get_graph_class(g_id));
	
}

int weight_sample(std::vector<double> prob){
	double sum = 0.0;
	for(auto x : prob){
		sum+=x;
	}
	// Uniform sample
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> distr(0.0, sum);
	double cut = distr(gen);
	double acu = 0;
	int pos=0;
	for(auto x : prob){
		acu += x;
		if(cut<=acu){return(pos);		}
		pos++;
	}
	return -1;	
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>
make_blobs(ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env, int size, int girth, std::vector<double> p,
		std::map<std::string, std::map<std::string, std::vector<std::string>>> structure, bool keep_centers, bool ignore_duplicates=false){
	
	ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> res_env;
	// Each graph in env is a center
	// p is a vector of size 6 containing the probability of each edit operation: add_node, add_edge, transform_node, transform_edge, remove_node, remove_edge

	std::pair<ged::GEDGraph::GraphID, ged::GEDGraph::GraphID>ids = env.graph_ids();
	int pos;
	std::string node_id;
	ged::ExchangeGraph<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> center;
	for(ged::GEDGraph::GraphID id{ids.first}; id<ids.second; id++){
		if(keep_centers){
			res_env.load_exchange_graph(env.get_graph(id, true, true, true), ged::undefined(), ged::Options::ExchangeGraphType::EDGE_LIST, env.get_graph_name(id), env.get_graph_class(id));
		}
		// Make blob for the center
		for(int i = 0; i<size; i++){
			center = env.get_graph(id, true, true, true);
			//Make modifications to the center graph
			int new_node_ids = 0;		
			for(int j=0; j<girth; j++){				
				//Choose operation to make
				pos = weight_sample(p);
				switch(pos){
					case 0:
						node_id = "new" + std::to_string(new_node_ids);
						edit_add_node<UserNodeID, UserNodeLabel, UserEdgeLabel>(center, node_id, structure);					
						new_node_ids++;
						break;
					case 1:
						edit_add_edge<UserNodeID, UserNodeLabel, UserEdgeLabel>(center, structure, ignore_duplicates);
						break;
					case 2:
						edit_transform_node<UserNodeID, UserNodeLabel, UserEdgeLabel>(center, structure);
						break;
					case 3:
						edit_transform_edge<UserNodeID, UserNodeLabel, UserEdgeLabel>(center, structure);
						break;
					case 4:
						edit_remove_node<UserNodeID, UserNodeLabel, UserEdgeLabel>(center);
						break;
					case 5:
						edit_remove_edge<UserNodeID, UserNodeLabel, UserEdgeLabel>(center);
						break;
					default:
						std::cout<<"Error in make_blobs: index "<<pos<<" out of bounds"<<"\n";
				}
			}
			res_env.load_exchange_graph(center, ged::undefined(), ged::Options::ExchangeGraphType::EDGE_LIST, env.get_graph_name(id), env.get_graph_class(id));
		}
	}
	return(res_env);
}




// TAKEN FROM EXAMPLES TO MAKE THE CODE WORK
std::unordered_set<std::string> irrelevant_node_attributes(const std::string & dataset){
	std::unordered_set<std::string> irrelevant_attributes;
	if (dataset == "AIDS") {
		irrelevant_attributes.insert({"x", "y", "symbol", "charge"});
	}
	return irrelevant_attributes;
}

bool constant_node_costs(const std::string & dataset) {
	if (dataset == "Letter") {
		return false;
	}
	else {
		return true;
	}
}

int main(int argc, char* argv[]){

	std::cout<<"MAIN\n";

	/* Test # 1. Make GEDLIB work, manipulate a graph and calculate a median for a known dataset
	 * 
	 * Works

	std::cout<<"--------------START-----------------\n";
	
	std::string gedlib_root("/home/lucas/Documents/stage_gedlibpy/gedlib/gedlib");
	std::string project_root("/home/lucas/Documents/stage_gedlibpy/stage/cpp");

	std::cout<<"--------------LOAD GRAPHS-----------------\n";
	
	
	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;

	std::string dataset{"Letter"};
	std::string class_test{"A"};
	std::string extra_dir{"/LOW"};

	std::string collection_file(project_root + "/data/collections/" + dataset + "_" + class_test +".xml");
	std::string graph_dir(gedlib_root + "/data/datasets/" + dataset + extra_dir);

	std::cout<<"Dataset: "<<dataset<<"\n";
	std::cout<<"Collection file: "<<collection_file<<"\n";
	std::cout<<"Graph directory: "<<graph_dir<<"\n";

	
	env.set_edit_costs(ged::Options::EditCosts::CONSTANT, {});
	
	
	std::vector<ged::GEDGraph::GraphID> graph_ids(env.load_gxl_graphs(graph_dir, collection_file,
			ged::Options::GXLNodeEdgeType::LABELED, ged::Options::GXLNodeEdgeType::LABELED, irrelevant_node_attributes(dataset)));
	ged::GEDGraph::GraphID median_id{env.add_graph("median")};
	env.init(ged::Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES);
	
	std::string ipfp_options("--threads 6 --initial-solutions 5 --initialization-method RANDOM");
	env.set_method(ged::Options::GEDMethod::IPFP, ipfp_options);



	std::cout<<"--------------DESCRIBE GRAPH-----------------\n";
	describe_graph(env.get_graph(graph_ids.at(0), true, false, false));

	
	std::cout<<"--------------FIND MEDIAN FOR LOADED GRAPHS-----------------\n";
	
	// Set up the estimator.
	ged::MedianGraphEstimator<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> mge(&env, constant_node_costs(dataset));
	mge.set_refine_method(ged::Options::GEDMethod::IPFP, "--threads 6 --initial-solutions 10 --ratio-runs-from-initial-solutions .5");
	

	// Varied estimator parameters.
	std::vector<std::string> init_types{"RANDOM", "MAX", "MIN", "MEAN", "MEDOID"};
	std::vector<std::string> nums_inits{"16", "1", "2", "4", "8", "32"};

	// Varied algorithm parameters.
	std::vector<ged::Options::GEDMethod> algos{ged::Options::GEDMethod::IPFP, ged::Options::GEDMethod::BRANCH_FAST, ged::Options::GEDMethod::REFINE};
	std::vector<std::string> algo_options_suffixes{" --initial-solutions 10 --ratio-runs-from-initial-solutions .5", "", " --initial-solutions 10 --ratio-runs-from-initial-solutions .5"};

	// Select the GED algorithm.
	std::size_t algo_id{0};

	std::string init_type = "RANDOM";

	std::string num_inits = "16";


	ged::Options::GEDMethod algo{algos.at(algo_id)};
	std::string algo_options("--threads 6" + algo_options_suffixes.at(algo_id));
	
	std::string mge_options("--time-limit 600 --stdout 0 --init-type " + init_type);
	if (init_type != "RANDOM" and num_inits != "1") {
	}
	else {
		std::random_device rng;
		mge_options += " --random-inits " + num_inits + " --randomness PSEUDO --seed " + std::to_string(rng());
	}


	mge.set_options(mge_options);
	mge.set_init_method(algo, algo_options);
	mge.set_descent_method(algo, algo_options);

	// Run the estimator.
	mge.run(graph_ids, median_id);

	// Write the results.
	std::cout << "\n Algo details" << init_type << ", " << num_inits << ", " << algo;
	std::cout << "\n Runtime" << mge.get_runtime() << ", " << mge.get_runtime(ged::Options::AlgorithmState::INITIALIZED) << ", " << mge.get_runtime(ged::Options::AlgorithmState::CONVERGED) <<"\n";
	std::cout << "\n SOD" << mge.get_sum_of_distances() << ", " << mge.get_sum_of_distances(ged::Options::AlgorithmState::INITIALIZED) << ", " << mge.get_sum_of_distances(ged::Options::AlgorithmState::CONVERGED)<<"\n";

	std::cout<<"--------------END-----------------\n";
	*/


	/* Test # 2 work with library edit_graph
	 *	Get graph structure and domain - works
	 *	Edit graph operations
	 *	Make_blobs
	 * 	Calculate median for blob
	 *	Calculate ged between blob and center, blob and median, median and center
	 *
	*/

	std::cout<<"--------------START-----------------\n";
	
	std::string gedlib_root("/home/lucas/Documents/stage_gedlibpy/gedlib/gedlib");
	std::string project_root("/home/lucas/Documents/stage_gedlibpy/stage/cpp");

	
	/*
	std::cout<<"--------------LOAD GRAPHS-----------------\n";
	
	
	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;

	std::string dataset{"Letter"};
	std::string class_test{"_A"};
	std::string extra_dir{"/LOW"};

	std::string collection_file(project_root + "/data/collections/" + dataset + class_test + ".xml");
	std::string graph_dir(gedlib_root + "/data/datasets/" + dataset + extra_dir);

	std::cout<<"Dataset: "<<dataset<<"\n";
	std::cout<<"Collection file: "<<collection_file<<"\n";
	std::cout<<"Graph directory: "<<graph_dir<<"\n";

	
	env.set_edit_costs(ged::Options::EditCosts::CONSTANT, {});
	
	
	std::vector<ged::GEDGraph::GraphID> graph_ids(env.load_gxl_graphs(graph_dir, collection_file,
			ged::Options::GXLNodeEdgeType::LABELED, ged::Options::GXLNodeEdgeType::LABELED, irrelevant_node_attributes(dataset)));
	ged::GEDGraph::GraphID median_id{env.add_graph("median")};
	env.init(ged::Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES);

	
	std::map<std::string, std::map<std::string, std::vector<std::string>>> 
	structure = get_graphs_structure(env);

	std::cout<<"--------------GRAPH STRUCTURE-----------------\n";
	for(auto t: structure){
		std::cout<<t.first<<"\n";
		for(auto att: structure.at(t.first)){
			std::cout<<"\t"<<att.first<<":\n\t\t";
			for(std::size_t i{0}; i<5; i++){
				std::cout<<att.second.at(i)<<" ";
			}
			std::cout<<"\n";
		}
	}
	*/
	/*
	ged::ExchangeGraph<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> g;

	std::cout<<"--------------EDIT_ Add Node-----------------\n";
	std::cout<<"Before edit\n";
	g = env.get_graph(0, true, true, true);
	describe_graph(g);
	edit_add_node<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>(g, "new_node_1", structure);
	std::cout<<"After edit\n";
	describe_graph(g);
	
	
	
	std::cout<<"--------------EDIT_ Add Edge-----------------\n";
	std::cout<<"Before edit\n";
	g = env.get_graph(0, true, true, true);
	describe_graph(g);
	edit_add_edge<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>(g, structure, false);
	std::cout<<"After edit\n";
	describe_graph(g);

	
	std::cout<<"--------------EDIT_ Transform Node-----------------\n";
	std::cout<<"Before edit\n";
	g = env.get_graph(0, true, true, true);
	describe_graph(g);
	edit_transform_node<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>(g, structure);
	std::cout<<"After edit\n";
	describe_graph(g);
	

	
	std::cout<<"--------------EDIT_ Transform Edge-----------------\n";
	std::cout<<"Before edit\n";
	g = env.get_graph(0, true, true, true);
	describe_graph(g);
	edit_transform_edge<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>(g, structure);
	std::cout<<"After edit\n";
	describe_graph(g);
	

	std::cout<<"--------------EDIT_ Remove_Node-----------------\n";
	std::cout<<"Before edit\n";
	g = env.get_graph(0, true, true, true);
	describe_graph(g);
	edit_remove_node<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>(g);
	std::cout<<"After edit\n";
	describe_graph(g);


	std::cout<<"--------------EDIT_ Remove_Edge-----------------\n";
	std::cout<<"Before edit\n";
	g = env.get_graph(0, true, true, true);
	describe_graph(g);
	edit_remove_edge<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>(g);
	std::cout<<"After edit\n";
	describe_graph(g);
	*/

	std::cout<<"--------------MAKE BLOBS-----------------\n";
	

	// For the centers
	class_test = "_A_one";
	collection_file = project_root + "/data/collections/" + dataset + class_test + ".xml";
	graph_dir = gedlib_root + "/data/datasets/" + dataset + extra_dir;

	std::cout<<"Dataset: "<<dataset<<"\n";
	std::cout<<"Collection file: "<<collection_file<<"\n";
	std::cout<<"Graph directory: "<<graph_dir<<"\n";

	env = ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>();
	env.set_edit_costs(ged::Options::EditCosts::CONSTANT, {});
	
	
	graph_ids = env.load_gxl_graphs(graph_dir, collection_file,
			ged::Options::GXLNodeEdgeType::LABELED, ged::Options::GXLNodeEdgeType::LABELED, irrelevant_node_attributes(dataset));
	env.init(ged::Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES);

	int size = 1;
	int girth = 1;
	std::vector<double> p {1.0,1.0,1.0,1.0,1.0,1.0};
	bool keep_centers = false;
	if(argc>1){
		size = std::stoi(argv[1]);
	}
	if(argc>2){
		girth = std::stoi(argv[2]);
	}
	if(argc>3){
		p = std::vector<double>();
		for(int k=0; k<6;k++){	
			p.emplace_back(std::stoi(argv[k+3]) + 0.0);
		}
	}
	if(argc>9){
		if(std::stoi(argv[9])==1){
			keep_centers = true;
		}
	}
	bool ignore_duplicates = false;

	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> blobs;

	ged::GEDGraph::GraphID id_center = graph_ids.at(0);


	// Create blobs without centers
	
	keep_centers = false;
	std::cout<<"\t\t\t\tNumber of graphs BEFORE: "<<env.num_graphs()<<"\n";
	blobs = make_blobs<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>(env, size, girth, p, structure, keep_centers, ignore_duplicates);
	env.init(ged::Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES);
	std::cout<<"\t\t\t\tNumber of graphs AFTER: "<<blobs.num_graphs()<<"\n";


	std::cout<<"--------------COMPUTE MEDIAN-----------------\n";
	
	blobs.init(ged::Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES);
	std::pair<ged::GEDGraph::GraphID, ged::GEDGraph::GraphID> ids_lim = blobs.graph_ids();
	graph_ids = std::vector<ged::GEDGraph::GraphID>();
	for(ged::GEDGraph::GraphID i{ids_lim.first}; i<ids_lim.second; i++){
		graph_ids.emplace_back(i);
	}
	median_id = blobs.add_graph("median");
	blobs.init(ged::Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES);


// Set up the estimator.
	ged::MedianGraphEstimator<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> mge(&blobs, constant_node_costs(dataset));
	mge.set_refine_method(ged::Options::GEDMethod::IPFP, "--threads 6 --initial-solutions 10 --ratio-runs-from-initial-solutions .5");
	

	// Varied estimator parameters.
	std::vector<std::string> init_types{"RANDOM", "MAX", "MIN", "MEAN", "MEDOID"};
	std::vector<std::string> nums_inits{"16", "1", "2", "4", "8", "32"};

	// Varied algorithm parameters.
	std::vector<ged::Options::GEDMethod> algos{ged::Options::GEDMethod::IPFP, ged::Options::GEDMethod::BRANCH_FAST, ged::Options::GEDMethod::REFINE};
	std::vector<std::string> algo_options_suffixes{" --initial-solutions 10 --ratio-runs-from-initial-solutions .5", "", " --initial-solutions 10 --ratio-runs-from-initial-solutions .5"};

	// Select the GED algorithm.
	std::size_t algo_id{0};

	std::string init_type = "RANDOM";

	std::string num_inits = "16";


	ged::Options::GEDMethod algo{algos.at(algo_id)};
	std::string algo_options("--threads 6" + algo_options_suffixes.at(algo_id));
	
	std::string mge_options("--time-limit 600 --stdout 0 --init-type " + init_type);
	if (init_type != "RANDOM" and num_inits != "1") {
	}
	else {
		std::random_device rng;
		mge_options += " --random-inits " + num_inits + " --randomness PSEUDO --seed " + std::to_string(rng());
	}


	mge.set_options(mge_options);
	mge.set_init_method(algo, algo_options);
	mge.set_descent_method(algo, algo_options);

	// Run the estimator.
	mge.run(graph_ids, median_id);

	// Write the results.
	std::cout << "\n Algo details" << init_type << ", " << num_inits << ", " << algo;
	std::cout << "\n Runtime" << mge.get_runtime() << ", " << mge.get_runtime(ged::Options::AlgorithmState::INITIALIZED) << ", " << mge.get_runtime(ged::Options::AlgorithmState::CONVERGED) <<"\n";
	std::cout << "\n SOD" << mge.get_sum_of_distances() << ", " << mge.get_sum_of_distances(ged::Options::AlgorithmState::INITIALIZED) << ", " << mge.get_sum_of_distances(ged::Options::AlgorithmState::CONVERGED)<<"\n";

}




