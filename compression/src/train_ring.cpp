#include <iostream>
#define GXL_GEDLIB_SHARED
#undef GXL_GEDLIB_SHARED
#define COMPRESS_EDIT_COST
//#include "/home/lucas/Documents/stage/gedlib/src/env/ged_env.hpp"
#include "src/env/ged_env.hpp"
#include "median/src/median_graph_estimator.hpp"
#undef COMPRESS_EDIT_COST
#include "compression/ext/frangio68/Minimal-Spanning-Arborescence-master/MSArbor.C"
#include <random>
#include <string>
#include <cmath>
#define BITS_IN_BYTE 8
#define bytes(num) ceil(log2(num+1)/BITS_IN_BYTE)
#define NBTHREADS_COMPR 2



template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void get_graphs_structure(ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env,
	std::map<std::string, std::map<std::string, std::vector<std::string>>> &distribution,
	std::map<std::string, std::map<std::string, std::set<std::string>>> &alphabets,
	double &b_ni, double &b_na, double &b_ei, double &b_ea){


	std::map<std::string, std::vector<std::string>> node_attr;
	std::map<std::string, std::vector<std::string>> edge_attr;
	std::map<std::string, std::set<std::string>> node_attr_set;
	std::map<std::string, std::set<std::string>> edge_attr_set;

	std::size_t max_nodes = 0;
	std::size_t max_edges = 0;

	ged::ExchangeGraph<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> g;
	std::pair<ged::GEDGraph::GraphID, ged::GEDGraph::GraphID> graph_ids = env.graph_ids();
	for (ged::GEDGraph::GraphID g_id{graph_ids.first}; g_id<graph_ids.second; g_id++){
		g = env.get_graph(g_id);
		max_nodes = max(max_nodes, g.num_nodes);
		max_edges = max(max_edges, g.num_edges);
		for (std::size_t i{0}; i < g.num_nodes; i++) {
			for(auto l : g.node_labels.at(i)){
				if(node_attr.count(l.first)==0){
					node_attr.emplace(std::make_pair(l.first, std::vector<std::string>()));
					node_attr_set.emplace(std::make_pair(l.first, std::set<std::string>()));
				}		
				node_attr.at(l.first).emplace_back(l.second);
				node_attr_set.at(l.first).insert(l.second);
			}
		}
		for (std::size_t i{0}; i < g.num_nodes; i++) {
			for (std::size_t j{i + 1}; j < g.num_nodes; j++) {
				if (g.adj_matrix[i][j] == 1) {				
					for(auto e : g.edge_labels.at(std::make_pair(i, j))){
						if(edge_attr.count(e.first)==0){
							edge_attr.emplace(std::make_pair(e.first, std::vector<std::string>()));
							edge_attr_set.emplace(std::make_pair(e.first, std::set<std::string>()));
						}
						edge_attr.at(e.first).emplace_back(e.second);
						edge_attr_set.at(e.first).insert(e.second);
					} 
				}
			}
		}
	}

	distribution.clear();
	distribution.emplace(std::make_pair("node_attr", node_attr));
	distribution.emplace(std::make_pair("edge_attr", edge_attr));

	alphabets.clear();
	alphabets.emplace(std::make_pair("node_attr", node_attr_set));
	alphabets.emplace(std::make_pair("edge_attr", edge_attr_set));

	b_ni = bytes(max_nodes);
	b_ei = bytes(max_edges);
	b_na = 0;
	b_ea = 0;

	for(auto const& a: alphabets.at("node_attr")){
		b_na += bytes(a.second.size());
	}
	for(auto const& a: alphabets.at("edge_attr")){
		b_ea += bytes(a.second.size());
	}
	
}

std::string init_options(const std::string & path, const std::string & data_suffix = "", bool save_train = false, bool load_train = false, std::size_t threads = 10) {
	std::string options("--threads ");
	options += std::to_string(threads);
	if (save_train) {
		if (load_train) {
			throw ged::Error("Training data cannot be both saved and loaded.");
		}
		options += " --save " + path + "/ring_" + data_suffix + ".data";
	}
	if (load_train) {
		options += " --load " + path + "/ring_" + data_suffix + ".data";
	}
	return options;
}

void train_ring(std::string graph_dir, std::string collection, std::string train_path, std::string suffix){

	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env_train;

	std::cout<<"Loading graphs"<<std::endl;
	std::cout<<"Graph directory: "<<graph_dir<<std::endl;
	std::cout<<"Collection: "<<collection<<std::endl;
	std::vector<ged::GEDGraph::GraphID> graph_ids_train(env_train.load_gxl_graphs(graph_dir, 
	 collection,
	  ged::Options::GXLNodeEdgeType::LABELED, ged::Options::GXLNodeEdgeType::UNLABELED));


	std::map<std::string, std::map<std::string, std::vector<std::string>>> distribution;
	std::map<std::string, std::map<std::string, std::set<std::string>>> alphabets;
	
	double b_ni;
	double b_na;
	double b_ei; 
	double b_ea;

	get_graphs_structure<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>(env_train, distribution, alphabets, b_ni, b_na, b_ei, b_ea);

	std::vector<double> comp_costs;
	double c_nd, c_ni, c_ns, c_ed, c_ei, c_es, c_es_id;
	
	// Compression costs: first version (eq 22 - 25)
	c_ni = b_na;
	c_nd = b_ni;
	c_ns = b_ni + b_na;
	c_ei = 2*b_ni + b_ea;
	c_ed = b_ei;
	c_es = b_ei + b_ea;
	c_es_id = 0;
	comp_costs.emplace_back(c_ni);
	comp_costs.emplace_back(c_nd);
	comp_costs.emplace_back(c_ns);
	comp_costs.emplace_back(c_ei);
	comp_costs.emplace_back(c_ed);
	comp_costs.emplace_back(c_es);
	comp_costs.emplace_back(c_es_id);

	env_train.set_edit_costs(ged::Options::EditCosts::COMPRESSION, comp_costs);
	env_train.init(ged::Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES);


	// Learn the parameters.
	std::vector<std::string> led_methods{"GAMMA", "LSAPE_GREEDY", "LSAPE_OPTIMAL"};
	for (auto led_method : led_methods) {
		std::cout << "\n=== " << led_method << " ===\n";
		env_train.set_method(ged::Options::GEDMethod::RING, init_options(train_path, led_method + "_" + suffix, true, false, 24) +  " --led-method " + led_method);
		env_train.init_method();
	
	}					

}

		
int main(int argc, char* argv[]){

	std::map<std::string, std::string> args;
	
	std::string graph_dir;
	std::string collection;
	std::string suffix;
	std::string train_path;

	
	if(argc>1) collection = argv[1];
	if(argc>2) graph_dir = argv[2];
	if(argc>3) suffix = argv[3];
	if(argc>4) train_path = argv[4];

	
	std::ifstream input_collection(collection.c_str());
	
	std::ifstream input_graph_dir(graph_dir.c_str());
	
	std::ifstream input_suffix(suffix.c_str());
	
	std::ifstream input_train_path(train_path.c_str());


    if (!input_graph_dir || !input_collection || !input_suffix || !input_train_path) {
        std::cout << "Unable to open files"<<std::endl;
        exit(1); // terminate with error
    }

    while (input_graph_dir >> graph_dir) {
        input_collection >> collection;
        input_suffix >> suffix;
        input_train_path >> train_path;
        train_ring(graph_dir, collection, train_path, suffix);
    }
   
	return 0;
	
}


