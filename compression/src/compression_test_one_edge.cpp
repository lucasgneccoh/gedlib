#include <iostream>
#define GXL_GEDLIB_SHARED
#undef GXL_GEDLIB_SHARED
#define COMPRESS_EDIT_COST
//#include "/home/lucas/Documents/stage/gedlib/src/env/ged_env.hpp"
#include "src/env/ged_env.hpp"
#include "median/src/median_graph_estimator.hpp"
#undef COMPRESS_EDIT_COST
#include <random>
#include <string>
#include <cmath>
#define BITS_IN_BYTE 8
#define bytes(num) ceil(log2(num+1)/BITS_IN_BYTE)


void describe_graph(ged::ExchangeGraph<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> exchange_graph){
	std::cout<<"NODES:"<<std::endl;
	for (std::size_t i{0}; i < exchange_graph.num_nodes; i++) {
		std::cout << "Node " << i << " -- (" <<exchange_graph.original_node_ids.at(i) <<") :"<<std::endl;
		for(auto const &l : exchange_graph.node_labels.at(i)){
			std::cout<<"\t"<<l.first<<": " <<l.second << std::endl;
		}
		
	}
	std::cout<<"EDGES:"<<std::endl;
	for (std::size_t i{0}; i < exchange_graph.num_nodes; i++) {
		for (std::size_t j{i + 1}; j < exchange_graph.num_nodes; j++) {
			if (exchange_graph.adj_matrix[i][j] == 1) {
				std::cout << "Edge:  " <<i << " -- "<<j<< " -> " << exchange_graph.original_node_ids.at(i);
				std::cout << " -- " << exchange_graph.original_node_ids.at(j) << std::endl;
				for(auto const &e : exchange_graph.edge_labels.at(std::make_pair(i, j))){
					std::cout<<"\t"<<e.first<<": " <<e.second << std::endl;
				} 
			}
		}
	}
}

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

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
edit_add_node(ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &exchange_graph,UserNodeID node_id, std::map<std::string, std::map<std::string, std::vector<std::string>>> structure){
	
	UserNodeLabel label;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<int> distr(0, 1);

		// Create the node label according to the graph structure
	for(auto const attr : structure.at("node_attr")){
		distr = std::uniform_int_distribution<int>(0, attr.second.size()-1);		
		label.emplace(std::make_pair(attr.first, attr.second.at(distr(gen))));
	}

		exchange_graph.original_node_ids.emplace_back(node_id);
		exchange_graph.node_labels.emplace_back(label);

		std::vector<std::size_t> last_line;
	for (std::size_t i{0}; i < exchange_graph.num_nodes; i++){
		exchange_graph.adj_matrix.at(i).emplace_back(0);
		last_line.emplace_back(0);
	}

		last_line.emplace_back(0);

		exchange_graph.adj_matrix.emplace_back(last_line);

	//env.add_node(g_id, node_id, label);
	exchange_graph.num_nodes++;
	return;

}


template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
edit_add_edge(ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &exchange_graph,std::map<std::string, std::map<std::string, std::vector<std::string>>> structure, bool ignore_duplicates){
	// We create the UserNodeLabel (GXLLabel in our case) and use the existing env.add_node method
	if(exchange_graph.num_nodes<2){
		return;
	}
	UserEdgeLabel label;
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<std::size_t> distr(0, 1);

		// Create the edge label according to the graph structure
	for(auto const attr : structure.at("edge_attr")){
		distr = std::uniform_int_distribution<std::size_t>(0, attr.second.size()-1);
		label.emplace(std::make_pair(attr.first, attr.second.at(distr(gen))));
	}

	std::pair<std::size_t, std::size_t> e_add;
	if(!ignore_duplicates){
		std::vector<std::pair<std::size_t, std::size_t>> non_edges;
				// Get non edges
		for (std::size_t i{0}; i < exchange_graph.num_nodes; i++) {
			for (std::size_t j{i + 1}; j < exchange_graph.num_nodes; j++) {
				if (exchange_graph.adj_matrix[i][j] == 0) {	
					non_edges.emplace_back(std::make_pair(i,j));
				}
			}
		}
				// Chose the non edge to add
		distr = std::uniform_int_distribution<std::size_t>(0, non_edges.size()-1);
		e_add = non_edges.at(distr(gen));
	}
	else{
				// Does not matter if there are multi-edges. We pick two random nodes
		distr = std::uniform_int_distribution<std::size_t>(0, exchange_graph.num_nodes-1);
		std::size_t i = distr(gen);
		std::size_t j = distr(gen);
		// For now, avoid loops
		while(i != j){
			j = distr(gen);
		}
		e_add = std::make_pair(i,j);
	}



	//env.add_edge(g_id, from_id, to_id, label, ignore_duplicates);

	// Modify adj_matrix and edge_list
	exchange_graph.adj_matrix.at(e_add.first).at(e_add.second) = 1;
	exchange_graph.edge_list.emplace_back(std::make_pair(e_add, label));
	exchange_graph.edge_labels.emplace(std::make_pair(e_add, label));
	exchange_graph.num_edges++;
	return;

}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
edit_transform_node(ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &exchange_graph, std::map<std::string, std::map<std::string, std::vector<std::string>>> structure){

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<std::size_t> distr(0, exchange_graph.num_nodes-1);
	std::size_t sel_node = distr(gen);

	UserNodeLabel label = exchange_graph.node_labels.at(sel_node);
	// Transform label according to the structure of the graphs. Only modify existing attributes
	for(auto  attr : exchange_graph.node_labels.at(sel_node)){
		distr = std::uniform_int_distribution<std::size_t>(0, structure.at("node_attr").at(attr.first).size()-1);		
		exchange_graph.node_labels.at(sel_node).at(attr.first) = structure.at("node_attr").at(attr.first).at(distr(gen));
	}

	// Update graph in env. Keep same name and class
	//env.load_exchange_graph(exchange_graph, g_id, ged::Options::ExchangeGraphType::ADJ_MATRIX ,env.get_graph_name(g_id) , env.get_graph_class(g_id));
	return;
}


template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
edit_transform_edge(ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &exchange_graph, std::map<std::string, std::map<std::string, std::vector<std::string>>> structure){

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<std::size_t> distr;
	distr = std::uniform_int_distribution<std::size_t>(0, exchange_graph.num_edges-1);
	std::size_t sel_edge = distr(gen);

	typename std::list<std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel>>::iterator iter;
	iter = exchange_graph.edge_list.begin();
	advance(iter, sel_edge);
	std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel> edge = *iter;
	// Transform label according to the structure of the graphs. Only modify existing attributes
	for(auto attr : edge.second){
		distr = std::uniform_int_distribution<std::size_t>(0, structure.at("edge_attr").at(attr.first).size()-1);		
		(*iter).second.at(attr.first) = structure.at("edge_attr").at(attr.first).at(distr(gen));
		exchange_graph.edge_labels.at(edge.first) = edge.second;
	}

}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
edit_remove_node(ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &exchange_graph){
	
	if(exchange_graph.num_nodes==0){
		return;
	}

	if(exchange_graph.num_nodes==1){
		exchange_graph.num_nodes=0;
		exchange_graph.num_edges=0;
		exchange_graph.original_node_ids.clear();
		exchange_graph.node_labels.clear();
		exchange_graph.adj_matrix.clear();
		exchange_graph.edge_labels.clear();
		exchange_graph.adj_lists.clear();
		exchange_graph.edge_list.clear();
		return;
	}

	// Remove even if it has edges??
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<std::size_t> distr(0, exchange_graph.num_nodes-1);
	std::size_t n = distr(gen);
		exchange_graph.original_node_ids.erase(exchange_graph.original_node_ids.begin() + n);
		exchange_graph.node_labels.erase(exchange_graph.node_labels.begin() + n);

		// TODO: May cause problems if vectors become empty??
	typename std::vector<std::size_t>::iterator iter_row;
	for (ged::GEDGraph::NodeID i{0}; i < exchange_graph.num_nodes; i++){
		iter_row = exchange_graph.adj_matrix.at(i).begin();
		advance(iter_row, n);
		exchange_graph.adj_matrix.at(i).erase(iter_row);	
	}
		typename std::vector<std::vector<std::size_t>>::iterator iter_col;
		iter_col = exchange_graph.adj_matrix.begin();
		advance(iter_col, n);
		exchange_graph.adj_matrix.erase(iter_col);
		typename std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel> entry;
	for(typename std::list<std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel>>::iterator iter = exchange_graph.edge_list.begin();
		iter != exchange_graph.edge_list.end(); iter++){
				entry = *iter;
		if((entry.first.first == n) | (entry.first.second == n)){
			exchange_graph.edge_labels.erase(entry.first);
			iter = exchange_graph.edge_list.erase(iter);
			iter--;
		}
	}
		// re-index nodes in edge list
	typename std::map<std::pair<std::size_t, std::size_t>, UserEdgeLabel> new_edge_labels;
	for(typename std::list<std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel>>::iterator iter = exchange_graph.edge_list.begin();
		iter != exchange_graph.edge_list.end(); iter++){
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
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<std::size_t> distr(0, exchange_graph.num_edges-1);
	std::size_t n = distr(gen);	
	typename std::list<std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel>>::iterator iter;
	typename std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel> edge;
	iter = exchange_graph.edge_list.begin();
	advance(iter,n);
	edge = *iter;
	exchange_graph.adj_matrix.at(edge.first.first).at(edge.first.second)=0;
	exchange_graph.edge_labels.erase(edge.first);
	exchange_graph.edge_list.erase(iter);
	exchange_graph.num_edges--;
	//env.load_exchange_graph(exchange_graph, g_id, ged::Options::ExchangeGraphType::ADJ_MATRIX ,env.get_graph_name(g_id) , env.get_graph_class(g_id));
	
}

void weight_sample(std::vector<double> &prob, int &pos, double &cut, double &acu){	
	acu = 0;
	pos=0;
	for(std::size_t k=0; k<prob.size(); k++){
		acu = acu + prob.at(k);
		if(cut<=acu){
			
			return;
		}
		pos++;
	}
	pos = -1;
	return;	
}

// TODO: make blobs inside env, not a copy of it
// there was a problem with clear_graph
template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>
make_blobs(ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env, int size, int girth, std::vector<double> p,std::map<std::string, std::map<std::string, std::vector<std::string>>> structure, bool keep_centers, bool ignore_duplicates=false){
	
	//std::size_t tasks_todo = size * env.num_graphs();
	ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> res_env;

	std::random_device rd;
	std::mt19937 gen(rd());
	double sum = 0.0;
	for(std::size_t i=0; i<p.size(); i++){
		sum = sum + p.at(i);
	}
	std::uniform_real_distribution<double> distr(0.0, sum);
	double cut=0.0;
	double acu=0.0;
	int pos=0;
	
	// Each graph in env is a center
	// p is a vector of size 6 containing the probability of each edit operation: add_node, add_edge, transform_node, transform_edge, remove_node, remove_edge

	std::pair<ged::GEDGraph::GraphID, ged::GEDGraph::GraphID>ids = env.graph_ids();
	
	std::string node_id;
	ged::ExchangeGraph<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> center;
	//ged::ProgressBar progress(size * env.num_graphs() * girth);
	//std::cout << "Making blobs: " <<tasks_todo <<" graphs to do..." << std::endl;
	for(ged::GEDGraph::GraphID id{ids.first}; id<ids.second; id++){
		if(keep_centers){
			res_env.load_exchange_graph(env.get_graph(id, true, true, true), ged::undefined(), ged::Options::ExchangeGraphType::EDGE_LIST, env.get_graph_name(id), env.get_graph_class(id));
		}
		// Make blob for the center
		for(int graph_made = 0; graph_made<size; graph_made++){
			center = env.get_graph(id, true, true, true);
			//Make modifications to the center graph
			int new_node_ids = 0;		
			for(int modification=0; modification<girth; modification++){				
				//Choose operation to make
				//std::cout<<"Sample... ";
				cut = distr(gen);
				//std::cout<<cut;
				pos = 0;
				weight_sample(p, pos, cut, acu);
				//std::cout<<". Chosen: "<<pos<<" ... ";
				switch(pos){
					case 0:						
						node_id = "new_" + std::to_string(new_node_ids);
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
						std::cout<<"Error in make_blobs: index "<<pos<<" out of bounds"<<std::endl;
						break;
				}				
			}			
			res_env.load_exchange_graph(center, ged::undefined(), ged::Options::ExchangeGraphType::EDGE_LIST, env.get_graph_name(id), env.get_graph_class(id));
		}
	}
	//std::cout << "all done"<<std::endl;
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

bool compare_label(std::map<std::string, std::string> label_1, std::map<std::string, std::string> label_2){
	
	if (label_1.size()!=label_2.size()){
		return(false);
	}
	for(auto const& entry : label_1){
		if(label_2.count(entry.first)<1){
			return(false);
		}
		if(label_2.at(entry.first) != entry.second){
			return(false);
		}
	}
	return(true);
}

void get_deleted_nodes(ged::NodeMap node_map, vector<ged::GEDGraph::NodeID> &v_d){
	v_d.clear();
	std::vector<ged::NodeMap::Assignment> relation;
	node_map.as_relation(relation);
	for(auto assignment : relation){
		if(assignment.second == ged::GEDGraph::dummy_node()){
			v_d.emplace_back(assignment.first);
		}
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void get_added_nodes(ged::NodeMap node_map, vector<ged::GEDGraph::NodeID> &v_i, vector<UserNodeLabel> &varphi_i,
	ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &g2){
	v_i.clear();
	varphi_i.clear();
	std::vector<ged::NodeMap::Assignment> relation;
	node_map.as_relation(relation);
	for(auto assignment : relation){
		if(assignment.first == ged::GEDGraph::dummy_node()){
			v_i.emplace_back(assignment.second);
			varphi_i.emplace_back(g2.node_labels.at(assignment.second));
		}
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void get_substituted_nodes(ged::NodeMap node_map, vector<ged::GEDGraph::NodeID> &v_s, vector<ged::GEDGraph::NodeID> &v_is, vector<UserNodeLabel> &varphi_s,
	ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &g1, ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &g2){
	v_s.clear();
	v_is.clear();
	varphi_s.clear();
	std::vector<ged::NodeMap::Assignment> relation;
	node_map.as_relation(relation);
	for(auto assignment : relation){
		if(assignment.first != ged::GEDGraph::dummy_node() &&  assignment.second != ged::GEDGraph::dummy_node()){
			if(compare_label(g1.node_labels.at(assignment.first),g2.node_labels.at(assignment.second))){
				v_is.emplace_back(assignment.first);
			}			
			else{
				v_s.emplace_back(assignment.first);
				varphi_s.emplace_back(g2.node_labels.at(assignment.second));
			}			
		}
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void get_deleted_edges(ged::NodeMap node_map, vector<std::pair<ged::GEDGraph::NodeID, ged::GEDGraph::NodeID>> &e_nd,
	vector<std::pair<ged::GEDGraph::NodeID, ged::GEDGraph::NodeID>> &e_ed,
	ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &g1, ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &g2){
	e_nd.clear();
	e_ed.clear();
	//std::cout<<"UNDEFINED NODE: "<<ged::GEDGraph::undefined_node()<<std::endl;
	//std::cout<<"DUMMY NODE: "<<ged::GEDGraph::dummy_node()<<std::endl;
	for(auto edge : g1.edge_list){
		//std::cout<<"\tEdge "<<edge.first.first<<" - "<<edge.first.second<<std::endl;
		//std::cout<<"\t"<<node_map.image(edge.first.first)<<" - " <<node_map.image(edge.first.second)<<std::endl;
		if(node_map.image(edge.first.first)==ged::GEDGraph::dummy_node() || node_map.image(edge.first.second)==ged::GEDGraph::dummy_node()){
			e_nd.emplace_back(edge.first);
		}
		else{
			if(g2.adj_matrix[node_map.image(edge.first.first)][node_map.image(edge.first.second)]==0){
				e_ed.emplace_back(edge.first);
			}			
		}
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void get_added_edges(ged::NodeMap node_map, vector<std::pair<ged::GEDGraph::NodeID, ged::GEDGraph::NodeID>> &e_ni,
	vector<std::pair<ged::GEDGraph::NodeID, ged::GEDGraph::NodeID>> &e_ei, vector<UserEdgeLabel> &phi_ni, vector<UserEdgeLabel> &phi_ei,
	ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &g1, ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &g2){
	e_ni.clear();
	e_ei.clear();
	phi_ni.clear();
	phi_ei.clear();
	for(auto edge : g2.edge_list){
		//std::cout<<"\tEdge "<<edge.first.first<<" - "<<edge.first.second<<std::endl;
		//std::cout<<"\t"<<node_map.image(edge.first.first)<<" - " <<node_map.image(edge.first.second)<<std::endl;
		if(node_map.pre_image(edge.first.first)==ged::GEDGraph::dummy_node() || node_map.pre_image(edge.first.second)==ged::GEDGraph::dummy_node()){
			e_ni.emplace_back(edge.first);
			phi_ni.emplace_back(edge.second);			
		}
		else{
			if(g1.adj_matrix[node_map.pre_image(edge.first.first)][node_map.pre_image(edge.first.second)]==0){
				e_ei.emplace_back(edge.first);
				phi_ei.emplace_back(edge.second);
			}
		}
	}	
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void get_substituted_edges(ged::NodeMap node_map, vector<std::pair<ged::GEDGraph::NodeID, ged::GEDGraph::NodeID>> &e_s,
	vector<std::pair<ged::GEDGraph::NodeID, ged::GEDGraph::NodeID>> &e_is, vector<UserEdgeLabel> &phi_s,
	ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &g1, ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &g2){	
	e_s.clear();
	e_is.clear();
	phi_s.clear();
	for(auto edge : g1.edge_list){
		//std::cout<<"\tEdge "<<edge.first.first<<" - "<<edge.first.second<<std::endl;
		//std::cout<<"\t"<<node_map.image(edge.first.first)<<" - " <<node_map.image(edge.first.second)<<std::endl;
		if(node_map.image(edge.first.first)!=ged::GEDGraph::dummy_node() && node_map.image(edge.first.second)!=ged::GEDGraph::dummy_node()){
			if(g2.adj_matrix[node_map.image(edge.first.first)][node_map.image(edge.first.second)]==1){
				if(compare_label(edge.second,g2.edge_labels.at(edge.first))){
					e_is.emplace_back(edge.first);
				}
				else{
					e_s.emplace_back(edge.first);
					phi_s.emplace_back(edge.second);
				}
			}			
		
		}
	}
}


template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void get_all_compression_sets(ged::NodeMap node_map,
	vector<ged::GEDGraph::NodeID> &v_d,
	vector<ged::GEDGraph::NodeID> &v_i,
	vector<UserNodeLabel> &varphi_i,
	vector<ged::GEDGraph::NodeID> &v_s,
	vector<ged::GEDGraph::NodeID> &v_is,
	vector<UserNodeLabel> &varphi_s,
	vector<std::pair<ged::GEDGraph::NodeID, ged::GEDGraph::NodeID>> &e_nd,
	vector<std::pair<ged::GEDGraph::NodeID, ged::GEDGraph::NodeID>> &e_ed,
	vector<std::pair<ged::GEDGraph::NodeID, ged::GEDGraph::NodeID>> &e_ni,
	vector<std::pair<ged::GEDGraph::NodeID, ged::GEDGraph::NodeID>> &e_ei,
	vector<UserEdgeLabel> &phi_ni,
	vector<UserEdgeLabel> &phi_ei,
	vector<std::pair<ged::GEDGraph::NodeID, ged::GEDGraph::NodeID>> &e_s,
	vector<std::pair<ged::GEDGraph::NodeID, ged::GEDGraph::NodeID>> &e_is,
	vector<UserEdgeLabel> &phi_s,
	ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &g1,
	ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &g2){
	
	//Nodes
	std::vector<ged::NodeMap::Assignment> relation;
	node_map.as_relation(relation);
	v_d.clear();
	v_i.clear();
	varphi_i.clear();
	v_s.clear();
	v_is.clear();
	varphi_s.clear();
	for(auto assignment : relation){
		if(assignment.second == ged::GEDGraph::dummy_node()){
			v_d.emplace_back(assignment.first);
		}
		if(assignment.first == ged::GEDGraph::dummy_node()){
			v_i.emplace_back(assignment.second);
			varphi_i.emplace_back(g2.node_labels.at(assignment.second));
		}
		if(assignment.first != ged::GEDGraph::dummy_node() &&  assignment.second != ged::GEDGraph::dummy_node()){
			if(compare_label(g1.node_labels.at(assignment.first),g2.node_labels.at(assignment.second))){
				v_is.emplace_back(assignment.first);
			}			
			else{
				v_s.emplace_back(assignment.first);
				varphi_s.emplace_back(g2.node_labels.at(assignment.second));
			}			
		}
	}
	
	// Edges
	e_nd.clear();
	e_ed.clear();
	e_ni.clear();
	e_ei.clear();
	phi_ni.clear();
	phi_ei.clear();
	e_s.clear();
	e_is.clear();
	phi_s.clear();
	for(auto const &edge : g1.edge_list){
		if(node_map.image(edge.first.first)==ged::GEDGraph::dummy_node() || node_map.image(edge.first.second)==ged::GEDGraph::dummy_node()){
			e_nd.emplace_back(edge.first);
		}
		else{
			if(g2.adj_matrix[node_map.image(edge.first.first)][node_map.image(edge.first.second)]==0){
				e_ed.emplace_back(edge.first);
			}
			else{
				if(compare_label(edge.second,g2.edge_labels.at(std::make_pair(node_map.image(edge.first.first),node_map.image(edge.first.second))))){
					e_is.emplace_back(edge.first);
				}
				else{
					e_s.emplace_back(edge.first);
					phi_s.emplace_back(edge.second);
				}
			}		
		}
	}
	
	for(auto const &edge2 : g2.edge_list){
		//std::cout<<"\tEdge "<<edge.first.first<<" - "<<edge.first.second<<std::endl;
		//std::cout<<"\t"<<node_map.image(edge.first.first)<<" - " <<node_map.image(edge.first.second)<<std::endl;
		if(node_map.pre_image(edge2.first.first)==ged::GEDGraph::dummy_node() || node_map.pre_image(edge2.first.second)==ged::GEDGraph::dummy_node()){
			e_ni.emplace_back(edge2.first);
			phi_ni.emplace_back(edge2.second);			
		}
		else{
			if(g1.adj_matrix[node_map.pre_image(edge2.first.first)][node_map.pre_image(edge2.first.second)]==0){
				e_ei.emplace_back(edge2.first);
				phi_ei.emplace_back(edge2.second);
			}
		}
	}
}
	
template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
std::pair<std::size_t, std::size_t> compression_size(ged::NodeMap node_map,
	ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &g1,
	ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &g2,
	double &b_ni, double &b_na,
	double &b_ei, double &b_ea){

	vector<ged::GEDGraph::NodeID> v_d;
	vector<ged::GEDGraph::NodeID> v_i;
	vector<UserNodeLabel> varphi_i;
	vector<ged::GEDGraph::NodeID> v_s;
	vector<ged::GEDGraph::NodeID> v_is;
	vector<UserNodeLabel> varphi_s;
	vector<std::pair<ged::GEDGraph::NodeID, ged::GEDGraph::NodeID>> e_nd;
	vector<std::pair<ged::GEDGraph::NodeID, ged::GEDGraph::NodeID>> e_ed;
	vector<std::pair<ged::GEDGraph::NodeID, ged::GEDGraph::NodeID>> e_ni;
	vector<std::pair<ged::GEDGraph::NodeID, ged::GEDGraph::NodeID>> e_ei;
	vector<UserEdgeLabel> phi_ni;
	vector<UserEdgeLabel> phi_ei;
	vector<std::pair<ged::GEDGraph::NodeID, ged::GEDGraph::NodeID>> e_s;
	vector<std::pair<ged::GEDGraph::NodeID, ged::GEDGraph::NodeID>> e_is;
	vector<UserEdgeLabel> phi_s;

	get_all_compression_sets(node_map,v_d,
		v_i,varphi_i,
		v_s,v_is,varphi_s,
		e_nd,e_ed,
		e_ni,e_ei,phi_ni,phi_ei,
		e_s,e_is,phi_s,
		g1,g2);
	
	std::size_t v_size=0;
	std::size_t e_size=0;
	std::size_t aux=0;
	// Vertex //Eq 44
	if(g1.num_nodes > g2.num_nodes){
		aux = b_ni * min(v_d.size(), v_is.size());
	}
	v_size = b_ni + b_na*g2.num_nodes - (b_na*v_is.size() - b_ni*(1+v_s.size() - aux));

	// Edges 
	//Eq 53
	e_size = 3*b_ei + b_ei*min(e_ed.size(), e_is.size()) + (b_ei + b_ea)*e_s.size() + (2*b_ni + b_ea)*(g2.num_edges-e_is.size()-e_s.size());
	//Eq 55
	/*
	if(e_ed.size() > e_is.size()){
		aux = (2*b_ni - b_ei + b_ea) * e_is.size() - 2*b_ei;
	}
	else{
		aux = (2*b_ni + b_ea) * e_is.size() - b_ei*e_ed.size() - 2*b_ei;	
	}
	e_size = b_ei + (2*b_ni + b_ea)*g2.num_edges - ( (2*b_ni - b_ei)*e_s.size() + aux );
	*/

	return(std::make_pair(v_size, e_size));
}


void print_compression_sets(
	vector<ged::GEDGraph::NodeID> &v_d,
	vector<ged::GEDGraph::NodeID> &v_i,
	vector<ged::GEDGraph::NodeID> &v_is,
	vector<ged::GEDGraph::NodeID> &v_s,
	vector<std::pair<ged::GEDGraph::NodeID, ged::GEDGraph::NodeID>> &e_nd,
	vector<std::pair<ged::GEDGraph::NodeID, ged::GEDGraph::NodeID>> &e_ed,
	vector<std::pair<ged::GEDGraph::NodeID, ged::GEDGraph::NodeID>> &e_ni,
	vector<std::pair<ged::GEDGraph::NodeID, ged::GEDGraph::NodeID>> &e_ei,
	vector<std::pair<ged::GEDGraph::NodeID, ged::GEDGraph::NodeID>> &e_is,
	vector<std::pair<ged::GEDGraph::NodeID, ged::GEDGraph::NodeID>> &e_s
	){
	std::cout<<"--------------NODES-----------------"<<std::endl;

	std::cout<<"Deleted nodes:"<<std::endl;
	for(auto n : v_d){
		std::cout<<n<<", ";
	}
	std::cout<<std::endl;


	std::cout<<"Added nodes:"<<std::endl;
	for(auto n : v_i){
		std::cout<<n<<", ";
	}
	std::cout<<std::endl;


	std::cout<<"Identically substituted nodes:"<<std::endl;
	for(auto n : v_is){
		std::cout<<n<<", ";
	}
	std::cout<<std::endl;


	std::cout<<"NON Identically substituted nodes:"<<std::endl;
	for(auto n : v_s){
		std::cout<<n<<", ";
	}
	std::cout<<std::endl;

	std::cout<<"--------------EDGES-----------------"<<std::endl;

	std::cout<<"Deleted edges, node deletion (e_nd):"<<std::endl;
	for(auto e : e_nd){
		std::cout<<e.first<<" -- "<<e.second<<", ";
	}
	std::cout<<std::endl;

	std::cout<<"Deleted edges, other (e_ed):"<<std::endl;
	for(auto e : e_ed){
		std::cout<<e.first<<" -- "<<e.second<<", ";
	}
	std::cout<<std::endl;

	std::cout<<"Added edges, node insertion (e_ni):"<<std::endl;
	for(auto e : e_ni){
		std::cout<<e.first<<" -- "<<e.second<<", ";
	}
	std::cout<<std::endl;

	std::cout<<"Added edges, other (e_ei):"<<std::endl;
	for(auto e : e_ei){
		std::cout<<e.first<<" -- "<<e.second<<", ";
	}
	std::cout<<std::endl;

	std::cout<<"Substituted edges, Identically (e_is):"<<std::endl;
	for(auto e : e_is){
		std::cout<<e.first<<" -- "<<e.second<<", ";
	}
	std::cout<<std::endl;

	std::cout<<"Substituted edges, NOT Identically (e_s):"<<std::endl;
	for(auto e : e_s){
		std::cout<<e.first<<" -- "<<e.second<<", ";
	}
	std::cout<<std::endl;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
ged::GEDGraph::GraphID add_shuffled_copy(
	ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env,
	ged::GEDGraph::GraphID graph_id
	){
	ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> g_ex = env.get_graph(graph_id, true, true, true);
	std::vector<std::pair<UserNodeID, UserNodeLabel>> nodes;
	for (std::size_t i = 0; i<g_ex.num_nodes; i++) {
		nodes.emplace_back(std::make_pair(g_ex.original_node_ids.at(i), g_ex.node_labels.at(i)));
	}
	ged::GEDGraph::GraphID copied_graph_id = env.add_graph(env.get_graph_name(graph_id) + "_shuffled_copy", env.get_graph_class(graph_id));
	std::random_device rng_g;
	std::mt19937 urng_g(rng_g());
	std::shuffle(nodes.begin(), nodes.end(), urng_g);
	for (const auto & node : nodes) {
		env.add_node(copied_graph_id, node.first, node.second);
	}

	for(typename std::list<std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel>>::iterator iter = g_ex.edge_list.begin();
		iter != g_ex.edge_list.end(); iter++){
		env.add_edge(copied_graph_id, 
			g_ex.original_node_ids.at((*iter).first.first),
			g_ex.original_node_ids.at((*iter).first.second), (*iter).second);
	}

	return copied_graph_id;
}


int main(int argc, char* argv[]){

	std::cout<<"MAIN\n";

	std::cout<<"--------------START-----------------\n";
	
	std::string gedlib_root("/home/lucas/Documents/stage_gedlibpy/gedlib/gedlib");
	std::string project_root("/home/lucas/Documents/stage_gedlibpy/stage/cpp");
	std::string dataset{"Letter"};
	std::string class_test{""};
	std::string extra_dir{"/HIGH"};

	std::string collection_file(gedlib_root + "/data/collections/" + dataset  + class_test +".xml");
	std::string graph_dir(gedlib_root + "/data/datasets/" + dataset + extra_dir);

	if(argc>8){
		collection_file = argv[8];
	}
	if(argc>9){
		graph_dir = argv[9];
	}

	std::cout<<"Collection file: "<<collection_file<<std::endl;
	std::cout<<"Graph directory: "<<graph_dir<<std::endl;


	std::cout<<"--------------LOAD GRAPHS, GET GRAPH STRUCTURE-----------------"<<std::endl;

	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;

		
	std::vector<ged::GEDGraph::GraphID> graph_ids(env.load_gxl_graphs(graph_dir, collection_file,
			ged::Options::GXLNodeEdgeType::LABELED, ged::Options::GXLNodeEdgeType::LABELED, irrelevant_node_attributes(dataset)));

	std::map<std::string, std::map<std::string, std::vector<std::string>>> distribution;
	std::map<std::string, std::map<std::string, std::set<std::string>>> alphabets;
	double b_ni;
	double b_na;
	double b_ei; 
	double b_ea;

	get_graphs_structure<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>(env, distribution, alphabets, b_ni, b_na, b_ei, b_ea);

	std::cout<<"Structure:"<<std::endl;
	std::cout<<"Nodes: "<<b_ni<<", "<<b_na<<std::endl;
	for(auto const& a: alphabets.at("node_attr")){
		std::cout<<"\t"<<a.first<<": "<<a.second.size()<<" values"<<std::endl;	
	}
	std::cout<<"Edges: "<<b_ei<<", "<<b_ea<<std::endl;
	for(auto const& a: alphabets.at("edge_attr")){
		std::cout<<"\t"<<a.first<<": "<<a.second.size()<<" values"<<std::endl;	
	}
	
	std::cout<<"b_ni: "<<b_ni<<", b_na: "<<b_na<<std::endl;
	std::cout<<"b_ei: "<<b_ei<<", b_ea: "<<b_ea<<std::endl;

	ged::ExchangeGraph<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> g_ex_1_orig;
	ged::ExchangeGraph<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> g_ex_1_copy;
	ged::ExchangeGraph<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> g_ex;
	ged::ExchangeGraph<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> g_ex_2;
	ged::ExchangeGraph<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> g_ex_3;
	ged::ExchangeGraph<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> g_empty;



	std::cout<<"--------------NODE MAP SECTION-----------------"<<std::endl;
	ged::GEDGraph::GraphID g1_id_orig = 0;
	ged::GEDGraph::GraphID g1_id = 0;
	ged::GEDGraph::GraphID g2_id = 0;
	ged::GEDGraph::GraphID g3_id = 0;
	ged::GEDGraph::GraphID empty_graph_id;

	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> blobs;
	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> aux_env;
	std::pair<ged::GEDGraph::GraphID, ged::GEDGraph::GraphID> ids;
	// Shuffle de graph
	g_ex_1_orig = env.get_graph(g1_id, true, true, true);
	
	// add_node, add_edge, transform_node, transform_edge, delete_node, delete_edge 
	std::vector<double> p{1,1,1,1,1,1};

	bool keep_centers = true;
	bool ignore_duplicates=false;
	int size = 1;
	int girth = 3;
	
	if(argc>1){
		p = std::vector<double>();
		for(int k=0; k<6;k++){	
			p.emplace_back(std::stoi(argv[k+1]) + 0.0);
		}
	}
	if(argc>7){
		girth = std::stoi(argv[7]);
	}
	int cost_version;
	int print_normal=0;
	int print_empty=0;
	int print_itself=0;

	
	if(argc>10){
		cost_version = std::stoi(argv[10]);
	}
	if(argc>11){
		print_normal = std::stoi(argv[11]);
	}
	if(argc>12){
		print_empty = std::stoi(argv[12]);
	}
	
	if(argc>13){
		print_itself = std::stoi(argv[13]);
	}
	int times = 1;
	
	if(argc>14){
		times = std::stoi(argv[14]);
	}
	int good = 0;
	int good_itself=0;
	int good_empty=0;
	std::vector<double> from_g1;
	std::vector<double> from_empty;
	std::vector<double> from_itself;

	std::vector<double> comp_costs;
	double c_nd, c_ni, c_ns, c_ed, c_ei, c_es, c_es_id;
	if(cost_version==0){
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
	
	}
	else{
		// Compression costs: second version (eq 34 - 38)
		double omega = 9999999999;

		if(g_ex.num_nodes > g_ex_2.num_nodes){
			c_nd = 0;
		}
		else{
			c_nd = omega;
		}

		if(g_ex.num_nodes < g_ex_2.num_nodes){
			c_ni = 0;
			c_ei = 0;
			c_es = b_ei - 2*b_ni;
			c_es_id = -(2*b_ni + b_ea);

		}
		else{
			c_ni = omega;
			c_ei = 2*b_ni + b_ea;
			c_es = b_ei + b_ea;	
			c_es_id = 0;
			
		}

		c_ns = b_ni + b_na;
		c_ed = b_ei;

		comp_costs.emplace_back(c_ni);
		comp_costs.emplace_back(c_nd);
		comp_costs.emplace_back(c_ns);
		comp_costs.emplace_back(c_ei);
		comp_costs.emplace_back(c_ed);
		comp_costs.emplace_back(c_es);
		comp_costs.emplace_back(c_es_id);
	}

	std::pair<std::size_t, std::size_t> comp_sizes;	


	vector<ged::GEDGraph::NodeID> v_d;
	vector<ged::GEDGraph::NodeID> v_i;
	vector<ged::GEDGraph::NodeID> v_s;
	vector<ged::GEDGraph::NodeID> v_is;
	vector<ged::GXLLabel> varphi_i;
	vector<ged::GXLLabel> varphi_s;
	vector<std::pair<ged::GEDGraph::NodeID, ged::GEDGraph::NodeID>> e_nd;
	vector<std::pair<ged::GEDGraph::NodeID, ged::GEDGraph::NodeID>> e_ed;
	vector<std::pair<ged::GEDGraph::NodeID, ged::GEDGraph::NodeID>> e_ni;
	vector<std::pair<ged::GEDGraph::NodeID, ged::GEDGraph::NodeID>> e_ei;
	vector<std::pair<ged::GEDGraph::NodeID, ged::GEDGraph::NodeID>> e_s;
	vector<std::pair<ged::GEDGraph::NodeID, ged::GEDGraph::NodeID>> e_is;
	vector<ged::GXLLabel> phi_ni;
	vector<ged::GXLLabel> phi_ei;
	vector<ged::GXLLabel> phi_s;
	vector<ged::GXLLabel> phi_is;


	
	for(int time=0; time<times; time++){

		std::cout<<time<<"."<<std::endl;
		ged::GEDGraph::GraphID shuffled_copy_id = 
		add_shuffled_copy<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>(env, g1_id_orig);
		
		g_ex_1_copy  = env.get_graph(shuffled_copy_id, true, true, true); 

		aux_env = ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>();

		env.clear_graph(shuffled_copy_id);

		g2_id = aux_env.load_exchange_graph(g_ex_1_copy, ged::undefined(), ged::Options::ExchangeGraphType::EDGE_LIST, "base", env.get_graph_class(shuffled_copy_id));

		
		blobs = ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>();
		
		blobs = make_blobs<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>(aux_env, size, girth, p, distribution, keep_centers, ignore_duplicates);
		
		// blobs has the two graphs
		ids = blobs.graph_ids();
		
		g2_id = ids.first;
		
		g3_id = ids.second-1;
		
		
		g1_id = blobs.load_exchange_graph(g_ex_1_orig, ged::undefined(), ged::Options::ExchangeGraphType::EDGE_LIST, env.get_graph_name(g1_id_orig), env.get_graph_class(g1_id_orig));

		g_ex  = blobs.get_graph(g1_id, true, true, true);
		g_ex_2  = blobs.get_graph(g2_id, true, true, true);
		g_ex_3  = blobs.get_graph(g3_id, true, true, true);
		if(times==1){
			blobs.save_as_gxl_graph(g1_id, "/home/lucas/Documents/stage/gedlib/compression/data/output/graph1.gxl");
			blobs.save_as_gxl_graph(g2_id, "/home/lucas/Documents/stage/gedlib/compression/data/output/graph2.gxl");
			blobs.save_as_gxl_graph(g3_id, "/home/lucas/Documents/stage/gedlib/compression/data/output/graph3.gxl");
		}
		// Describe the graphs:
		if(times==1){
			std::cout<<"G1: "<<std::endl;
			describe_graph(g_ex);
		}

		if(times==1){
			std::cout<<"G2: "<<std::endl;
			describe_graph(g_ex_2);	
		} 
		
		// Describe the two graphs:
		if(times==1){
			std::cout<<"G3: "<<std::endl;
			describe_graph(g_ex_3);
		}

		if(times==1) std::cout<<"--------------SET COMPRESSION EDIT COST-----------------"<<std::endl;

		
		if(times==1) std::cout<<"--------------INIT-----------------"<<std::endl;
		blobs.set_edit_costs(ged::Options::EditCosts::COMPRESSION, comp_costs);
		blobs.init(ged::Options::InitType::LAZY_WITHOUT_SHUFFLED_COPIES);


		// Set method
		//std::string bu_options("--lsape-model FLWC --threads 6 --greedy-method REFINED --optimal TRUE --centrality-method PAGERANK --max-num-solutions 4");
		std::string bu_options = "";
		blobs.set_method(ged::Options::GEDMethod::BRANCH_UNIFORM, bu_options);



		if(times==1) std::cout<<"--------------RUN METHOD-----------------"<<std::endl;
		blobs.run_method(g1_id, g2_id);
		blobs.run_method(g1_id, g3_id);

		

		if(times==1) std::cout<<"--------------GET NODE MAP-----------------"<<std::endl;
		ged::NodeMap node_map_itself = blobs.get_node_map(g1_id, g2_id);
		ged::NodeMap node_map = blobs.get_node_map(g1_id, g3_id);
		if(times == 1){
			std::cout<<"G_1 to G_2: "<<node_map_itself<<std::endl<<std::endl;
			std::cout<<"G_1 to G_3: "<<node_map<<std::endl<<std::endl;
		}

		if(times==1) std::cout<<"--------------GET INFO-----------------"<<std::endl;
		

		if(print_normal){
			std::cout<<"--------------COMPRESSION FROM G_1 to G_2-----------------"<<std::endl;
			get_all_compression_sets<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>(node_map_itself, v_d, v_i, varphi_i, 
				v_s, v_is, varphi_s, e_nd,e_ed, e_ni,e_ei, phi_ni, phi_ei, e_s,e_is, phi_s, g_ex, g_ex_2);

			print_compression_sets(v_d,v_i,v_is, v_s, e_nd, e_ed, e_ni, e_ei, e_is, e_s);
		}

		/*
		if(print_empty){
			std::cout<<"--------------COMPRESSION FROM EMPTY to G_2-----------------"<<std::endl;
			get_all_compression_sets<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>(node_map_empty, v_d, v_i, varphi_i, 
				v_s, v_is, varphi_s, e_nd,e_ed, e_ni,e_ei, phi_ni, phi_ei, e_s,e_is, phi_s, g_empty, g_ex_2);

			print_compression_sets(v_d,v_i,v_is, v_s, e_nd, e_ed, e_ni, e_ei, e_is, e_s);
		}
		*/

		if(print_normal){
			std::cout<<"--------------COMPRESSION FROM G_1 to G_3------------------"<<std::endl;
			get_all_compression_sets<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>(node_map, v_d, v_i, varphi_i, 
				v_s, v_is, varphi_s, e_nd,e_ed, e_ni,e_ei, phi_ni, phi_ei, e_s,e_is, phi_s, g_ex, g_ex_3);

			print_compression_sets(v_d,v_i,v_is, v_s, e_nd, e_ed, e_ni, e_ei, e_is, e_s);
		}

		


	 	
	 	if(print_normal){
			std::cout<<"--------------FROM G_1 to G_3-----------------"<<std::endl;
		 	comp_sizes = compression_size<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>(node_map, g_ex, g_ex_3, b_ni,b_na, b_ei, b_ea);
		 	std::cout<<"Vertex compression: "<<comp_sizes.first<<std::endl;
		 	std::cout<<"Edge compression: "<<comp_sizes.second<<std::endl;
			std::cout<<"Total compression: "<<comp_sizes.first + comp_sizes.second<<std::endl;

			std::cout<<"GED G_1 to G_3: "<<blobs.get_upper_bound(g1_id, g3_id)<<std::endl;
		}
		from_g1.emplace_back(blobs.get_upper_bound(g1_id, g3_id));
		if(blobs.get_upper_bound(g1_id, g3_id) == 3){
			good++;
		}
		
		/*
		if(print_empty){
			std::cout<<"--------------FROM  to G_2-----------------"<<std::endl;
			
			comp_sizes = compression_size<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>(node_map_empty, g_empty, g_ex_2, b_ni,b_na, b_ei, b_ea);
		 	std::cout<<"Vertex compression: "<<comp_sizes.first<<std::endl;
		 	std::cout<<"Edge compression: "<<comp_sizes.second<<std::endl;
			std::cout<<"Total compression: "<<comp_sizes.first + comp_sizes.second<<std::endl;

			std::cout<<"GED EMPTY to G_2: "<<blobs.get_upper_bound(empty_graph_id, g2_id)<<std::endl;
		}

		from_empty.emplace_back(blobs.get_upper_bound(empty_graph_id, g2_id));
		if(blobs.get_upper_bound(empty_graph_id, g2_id) == 123){
			good_empty++;
		}
		*/

		if(print_itself){
			std::cout<<"--------------FROM G_1 to G_2-----------------"<<std::endl;
			
			comp_sizes = compression_size<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>(node_map_itself, g_ex, g_ex_2, b_ni,b_na, b_ei, b_ea);
		 	std::cout<<"Vertex compression: "<<comp_sizes.first<<std::endl;
		 	std::cout<<"Edge compression: "<<comp_sizes.second<<std::endl;
			std::cout<<"Total compression: "<<comp_sizes.first + comp_sizes.second<<std::endl;

			std::cout<<"GED G_1 to G_2: "<<blobs.get_upper_bound(g1_id, g2_id)<<std::endl;
		}

		from_itself.emplace_back(blobs.get_upper_bound(g1_id, g2_id));
		if(blobs.get_upper_bound(g1_id, g2_id) == 0){
			good_itself++;
		}

		if(times==1) std::cout<<"--------------END FOR NOW-----------------"<<std::endl;
	}


	std::cout<<"Good g1 to g2: "<<good<<std::endl;	
	//std::cout<<"Good empty to g2: "<<good_empty<<std::endl;	
	std::cout<<"Good g2 to itself: "<<good_itself<<std::endl;	

	std::cout<<"--------------WRITE RESULTS-----------------"<<std::endl;
	// Generate the result file.
	std::string result_filename("/home/lucas/Documents/stage/gedlib/compression/data/output/");
	result_filename += "comp_cost_one_edge.csv";
	std::ofstream result_file(result_filename.c_str());
	result_file << "g1_to_g3,g1_to_g2\n";
	for(int time=0; time<times; time++){
		//result_file << from_g1.at(time) << "," << from_empty.at(time) << "," << from_itself.at(time) << "\n";
		result_file << from_g1.at(time) << ","<< from_itself.at(time) << "\n";
	}
	result_file.close();


	
	return 0;
	
}


