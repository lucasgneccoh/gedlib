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
#include <algorithm>
#define BITS_IN_BYTE 8
#define bytes(num) ceil(log2(num+1)/BITS_IN_BYTE)

// For the spanning arborescence
template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void describe_graph(ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> g){
	std::cout<<"NODES: "<<g.node_labels.size()<<std::endl;
	for (std::size_t i{0}; i < g.num_nodes; i++) {
		std::cout << "Node " << i << " -------- (" <<g.original_node_ids.at(i) <<") :"<<std::endl;
		for(auto const &l : g.node_labels.at(i)){
			std::cout<<"\t"<<l.first<<": " <<l.second << std::endl;
		}
		
	}
	std::cout<<"EDGES: "<<g.edge_list.size()<<std::endl;
	std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel> edge;

	typename std::list<std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel>>::iterator iter;
	for(iter = g.edge_list.begin(); iter != g.edge_list.end(); iter++){
		edge = (*iter);
		std::cout << "Edge:  " <<edge.first.first << " -- "<<edge.first.second << std::endl;
		for(auto const &e : edge.second){
			std::cout<<"\t"<<e.first<<": " <<e.second << std::endl;
		} 
	}

	std::cout<<"END"<<std::endl;
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

template<class UserEdgeLabel>
bool compare_edges(std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel> &e1, std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel> &e2){
	if(e1.first.first < e2.first.first)	return true;
	if(e1.first.first == e2.first.first && e1.first.second < e2.first.second) return true;
	return false;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void get_all_compression_sets(ged::NodeMap node_map,
	vector<ged::GEDGraph::NodeID> &v_d,
	vector<ged::GEDGraph::NodeID> &v_i,
	vector<UserNodeLabel> &varphi_i,
	vector<ged::GEDGraph::NodeID> &v_s,
	vector<ged::GEDGraph::NodeID> &v_is,
	vector<UserNodeLabel> &varphi_s,
	vector<std::size_t> &e_nd,
	vector<std::size_t> &e_ed,
	vector<std::size_t> &e_ni,
	vector<std::size_t> &e_ei,
	vector<UserEdgeLabel> &phi_ni,
	vector<UserEdgeLabel> &phi_ei,
	vector<std::size_t> &e_s,
	vector<std::size_t> &e_is,
	vector<UserEdgeLabel> &phi_s,
	ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &g1,
	ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &g2){
	
	std::size_t index;
	//Nodes
	//std::cout<<" compression sets: Nodes"<<std::endl;
	std::vector<ged::NodeMap::Assignment> relation;
	node_map.as_relation(relation);
	v_d.clear();
	v_i.clear();
	varphi_i.clear();
	v_s.clear();
	v_is.clear();
	varphi_s.clear();
	for(auto const assignment : relation){
		if(assignment.second == ged::GEDGraph::dummy_node()){
			v_d.emplace_back(assignment.first); // node index in G1
		}
		if(assignment.first == ged::GEDGraph::dummy_node()){
			v_i.emplace_back(assignment.second); // node index in G2
			varphi_i.emplace_back(g2.node_labels.at(assignment.second));
		}
		if(assignment.first != ged::GEDGraph::dummy_node() &&  assignment.second != ged::GEDGraph::dummy_node()){
			if(compare_label(g1.node_labels.at(assignment.first),g2.node_labels.at(assignment.second))){
				v_is.emplace_back(assignment.first); // node index in G1
			}			
			else{
				v_s.emplace_back(assignment.first); // node index in G1
				varphi_s.emplace_back(g2.node_labels.at(assignment.second));
			}			
		}
	}
	
	// Edges
	//std::cout<<" compression sets: Edges"<<std::endl;
	e_nd.clear();
	e_ed.clear();
	e_ni.clear();
	e_ei.clear();
	phi_ni.clear();
	phi_ei.clear();
	e_s.clear();
	e_is.clear();
	phi_s.clear();
	g1.edge_list.sort(compare_edges<UserEdgeLabel>);
	g2.edge_list.sort(compare_edges<UserEdgeLabel>);

	std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel> edge;
	typename std::list<std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel>>::iterator iter;
	index = 0;
	for(iter = g1.edge_list.begin(); iter != g1.edge_list.end(); iter ++ ){
		
		edge = (*iter);
		if(node_map.image(edge.first.first)==ged::GEDGraph::dummy_node() || node_map.image(edge.first.second)==ged::GEDGraph::dummy_node()){
			e_nd.emplace_back(index); // in G1.edge_list
		}
		else{
			if(g2.adj_matrix[node_map.image(edge.first.first)][node_map.image(edge.first.second)]==0){
				e_ed.emplace_back(index); // in G1.edge_list
			}
			else{
				if(compare_label(edge.second,g2.edge_labels.at(std::make_pair(node_map.image(edge.first.first),node_map.image(edge.first.second))))){
					e_is.emplace_back(index); // in G1.edge_list
				}
				else{
					e_s.emplace_back(index); // in G1.edge_list
					phi_s.emplace_back(g2.edge_labels.at(std::make_pair(node_map.image(edge.first.first),node_map.image(edge.first.second))));
				}
			}		
		}
		index++;
	}
	index = 0;
	std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel> edge2;

	for(iter = g2.edge_list.begin(); iter != g2.edge_list.end(); iter ++ ){
		edge2 = (*iter);
		if(node_map.pre_image(edge2.first.first)==ged::GEDGraph::dummy_node() || node_map.pre_image(edge2.first.second)==ged::GEDGraph::dummy_node()){
			e_ni.emplace_back(index); // in G2.edge_list
			phi_ni.emplace_back(edge2.second);			
		}
		else{
			if(g1.adj_matrix[node_map.pre_image(edge2.first.first)][node_map.pre_image(edge2.first.second)]==0){
				e_ei.emplace_back(index); // in G2.edge_list
				phi_ei.emplace_back(edge2.second);
			}
		}
		index++;
	}
	//std::cout<<" compression sets: END"<<std::endl;
}
	
std::pair<double, double> simplified_compression_size_from_empty_graph(
	std::size_t num_nodes, std::size_t num_edges, 
	double &b_ni, double &b_na,
	double &b_ei, double &b_ea){

	double v_size=0;
	double e_size=0;
	// Vertex //Eq 44

	v_size = 2 * b_ni + b_na*num_nodes;

	// Edges 
	//Eq 53
	e_size = 3 * b_ei  + (2*b_ni + b_ea) * num_edges;
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

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
double calculate_total_compression_cost(ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env, ged::GEDGraph::GraphID &empty_id, double &b_ni, double &b_na, double &b_ei, double &b_ea){
	double total_cost_compression = 0;
	std::pair<double,double> aux_pair;
	ged::ExchangeGraph<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> g_ex;
	std::pair<ged::GEDGraph::GraphID, ged::GEDGraph::GraphID> limits = env.graph_ids();
	for(ged::GEDGraph::GraphID i{limits.first}; i<limits.second; i++){
		if(i!=empty_id){
			g_ex = env.get_graph(i);
			aux_pair = simplified_compression_size_from_empty_graph(g_ex.num_nodes, g_ex.num_edges, b_ni, b_na, b_ei, b_ea);
			total_cost_compression += (aux_pair.first + aux_pair.second);
		}
	}

	return total_cost_compression;
}


template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
std::pair<double, double> compression_size(ged::NodeMap node_map,
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
	vector<std::size_t> e_nd;
	vector<std::size_t> e_ed;
	vector<std::size_t> e_ni;
	vector<std::size_t> e_ei;
	vector<UserEdgeLabel> phi_ni;
	vector<UserEdgeLabel> phi_ei;
	vector<std::size_t> e_s;
	vector<std::size_t> e_is;
	vector<UserEdgeLabel> phi_s;

	get_all_compression_sets<UserNodeID, UserNodeLabel, UserEdgeLabel>(node_map,v_d,
		v_i,varphi_i,
		v_s,v_is,varphi_s,
		e_nd,e_ed,
		e_ni,e_ei,phi_ni,phi_ei,
		e_s,e_is,phi_s,
		g1,g2);
	
	double v_size=0;
	double e_size=0;
	double aux=0;
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

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void print_compression_sets(
	vector<ged::GEDGraph::NodeID> &v_d,
	vector<ged::GEDGraph::NodeID> &v_i,
	vector<ged::GEDGraph::NodeID> &v_s,
	vector<ged::GEDGraph::NodeID> &v_is,
	vector<std::size_t> &e_nd,
	vector<std::size_t> &e_ed,
	vector<std::size_t> &e_ni,
	vector<std::size_t> &e_ei,	
	vector<std::size_t> &e_s,
	vector<std::size_t> &e_is,
	ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &g1,
	ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &g2
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

	g1.edge_list.sort(compare_edges<UserEdgeLabel>);
	g2.edge_list.sort(compare_edges<UserEdgeLabel>);

	std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel> edge;
	typename std::list<std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel>>::iterator iter;
	
	std::cout<<"Deleted edges, node deletion (e_nd):"<<std::endl;
	for(auto e : e_nd){
		iter = g1.edge_list.begin();
		std::advance (iter, e);
		edge = (*iter);
		std::cout<<edge.first.first<<" -- "<<edge.first.second<<", ";
	}
	std::cout<<std::endl;

	std::cout<<"Deleted edges, other (e_ed):"<<std::endl;
	for(auto e : e_ed){
		iter = g1.edge_list.begin();
		std::advance (iter, e);
		edge = (*iter);
		std::cout<<edge.first.first<<" -- "<<edge.first.second<<", ";
	}
	std::cout<<std::endl;

	std::cout<<"Substituted edges, Identically (e_is):"<<std::endl;
	for(auto e : e_is){
		iter = g1.edge_list.begin();
		std::advance (iter, e);
		edge = (*iter);
		std::cout<<edge.first.first<<" -- "<<edge.first.second<<", ";
	}
	std::cout<<std::endl;

	std::cout<<"Substituted edges, NOT Identically (e_s):"<<std::endl;
	for(auto e : e_s){
		iter = g1.edge_list.begin();
		std::advance (iter, e);
		edge = (*iter);
		std::cout<<edge.first.first<<" -- "<<edge.first.second<<", ";
	}
	std::cout<<std::endl;

	std::cout<<"Added edges, node insertion (e_ni):"<<std::endl;
	for(auto e : e_ni){
		iter = g2.edge_list.begin();
		std::advance (iter, e);
		edge = (*iter);
		std::cout<<edge.first.first<<" -- "<<edge.first.second<<", ";
	}
	std::cout<<std::endl;

	std::cout<<"Added edges, other (e_ei):"<<std::endl;
	for(auto e : e_ei){
		iter = g2.edge_list.begin();
		std::advance (iter, e);
		edge = (*iter);
		std::cout<<edge.first.first<<" -- "<<edge.first.second<<", ";
	}
	std::cout<<std::endl;

}

void spanning_arborescence_of_minimum_weight(
	std::vector<std::size_t> &tree,
	double &cost,
	std::vector<std::pair<std::size_t, std::size_t>> &edges,
	std::vector<std::vector<double>> &w,
	std::size_t &root,
	bool stdout = false
	){

	// FROM FISCHETTI 1993
	//https://github.com/frangio68/Minimal-Spanning-Arborescence
	if(stdout) std::cout<<"START ARBOR"<<std::endl;
	MSA_di_unipi_it::MSArbor::Index n = w.size();
	if(stdout) std::cout<<"Create csts"<<std::endl;
	MSA_di_unipi_it::MSArbor::CRow csts = new MSA_di_unipi_it::MSArbor::CNumber[ n * ( n - 1 ) ];

	if(stdout) std::cout<<"Fill csts"<<std::endl;
	for( MSA_di_unipi_it::MSArbor::Index i = 0 ; i < n ; i++ ){
		for( MSA_di_unipi_it::MSArbor::Index j = 0 ; j < n-1 ; j++ ){
			if(i == j){
				csts[ n * j + i ] = MSA_di_unipi_it::MSArbor::C_INF;
			}
			else{
				csts[ n * j + i ] = w.at(i).at(j);
			}
		}
	}
	if(stdout) std::cout<<"csts filled"<<std::endl;

	MSA_di_unipi_it::MSArbor MSA( n );

	cost = MSA.Solve( csts );

	if(stdout) std::cout<<"Solving"<<std::endl;
	if(stdout) std::cout << "n: " << MSA.GetN() << "\t Cost: " << cost << endl;


	tree.clear();
	if(stdout) std::cout<<"Creating results"<<std::endl;
		for(unsigned int i = 0 ; i < n - 1 ; i++ ){
		 	tree.emplace_back(MSA.ReadPred()[ i ]);
		}

	delete[] csts;

	return;
} 

void get_arborescence_info(
	std::vector<std::string> &headers,
	std::vector<std::string> &values,
	std::map<std::size_t, std::vector<std::size_t>> &depth_degrees,
	std::vector<std::size_t> &arborescence,
	std::size_t root){
	
	depth_degrees.clear();

	std::map<std::size_t, std::vector<std::size_t>> children;
	std::list<std::pair<std::size_t,std::size_t>> todo_depth;

	std::size_t num_nodes = arborescence.size()+1; 
	for(std::size_t i=0; i<num_nodes; i++){
		children.emplace(std::make_pair(i, std::vector<std::size_t>()));
	}	

	for(std::size_t i=0; i<arborescence.size(); i++){
		children.at(arborescence.at(i)).emplace_back(i);			
	}

	double avg_degree=0;
	std::size_t min_degree= std::numeric_limits<std::size_t>::max();
	std::size_t max_degree=0;
	double avg_depth=0;
	std::size_t max_depth=0;
	std::size_t leafs=0;

	for(auto const &entry : children){
		if(entry.second.size()==0){
			leafs++;
		}
		else{
			avg_degree += entry.second.size();
			if(entry.second.size()>max_degree) max_degree = entry.second.size();
			if(entry.second.size()<min_degree) min_degree = entry.second.size();	
		}
		
	}
	avg_degree = avg_degree/(num_nodes-leafs);


	std::size_t parent;
	std::size_t depth;

	for(std::size_t i=0; i<num_nodes; i++){
		parent = i;
		depth = 0;
		while(parent != root){
			parent = arborescence.at(parent);
			depth++;
		}
		avg_depth+=depth;
		if(depth>max_depth) max_depth = depth;
	}

	for(std::size_t i=0; i<max_depth; i++){
		depth_degrees.emplace(std::make_pair(i, std::vector<std::size_t>()));
	}

	todo_depth.push_front(std::make_pair(root, 0));
	std::pair<std::size_t,std::size_t> aux;
	while(! todo_depth.empty()){
		aux = todo_depth.front();
		todo_depth.pop_front();
		if(children.at(aux.first).size()==0){
			//leaf
			continue;
		}
		else{
			for(auto const &c : children.at(aux.first)){
				todo_depth.push_front(std::make_pair(c, aux.second+1));
			}
			depth_degrees.at(aux.second).emplace_back(children.at(aux.first).size());	
		}
	}

	avg_depth = avg_depth/num_nodes;

	headers.emplace_back("num_nodes");
	values.emplace_back(std::to_string(num_nodes));

	headers.emplace_back("num_leafs");
	values.emplace_back(std::to_string(leafs));

	headers.emplace_back("avg_degree");
	values.emplace_back(std::to_string(avg_degree));

	headers.emplace_back("min_degree");
	values.emplace_back(std::to_string(min_degree));

	headers.emplace_back("max_degree");
	values.emplace_back(std::to_string(max_degree));
	
	headers.emplace_back("root_degree");
	values.emplace_back(std::to_string(children.at(root).size()));

	headers.emplace_back("avg_depth");
	values.emplace_back(std::to_string(avg_depth));

	headers.emplace_back("max_depth");
	values.emplace_back(std::to_string(max_depth));

	/*
	std::cout<<"Avg degree by depth: "<<avg_depth<<std::endl;
	double sum=0;
	for(auto const &entry : depth_degrees){
		sum=0;

		for(auto const &d : entry.second){
			sum += d;
		}
		std::cout<<"Depth "<< entry.first <<" : "<< sum/entry.second.size() << std::endl;

	}
	*/

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

void write_to_file(std::ofstream &file, std::vector<std::string> &values){
	
	for(std::size_t j=0; j<values.size(); j++){
		file << values.at(j);
		if(j < values.size()-1){
			file << ","; 	
		} 
		else{
			file << "\n"; 	
		}
	}
}  


template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void encode_environment(ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env_coded, 
	std::map<std::string, std::map<std::string, std::map<std::string,std::string>>> &encoded_attributes,	
	std::string path, 
	std::string dataset, 
	ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env, 
	std::map<std::string, std::map<std::string, std::set<std::string>>> &alphabets,
	std::vector<std::size_t> &arborescence,
	std::size_t root,
	char escape_char = '#',
	char separator = '\n'
	){



	// Write alphabets in order and transform graphs
	std::ofstream ofile;
	ofile.open(path + "/" + dataset + ".info_file", ios::out | ios::binary);
	std::size_t cont=0;

	encoded_attributes.emplace("graphs", std::map<std::string, std::map<std::string,std::string>> ());
	encoded_attributes.emplace("node_attr", std::map<std::string, std::map<std::string,std::string>> ());
	encoded_attributes.emplace("edge_attr", std::map<std::string, std::map<std::string,std::string>> ());
	
	std::pair<ged::GEDGraph::GraphID, ged::GEDGraph::GraphID> limits;
	limits = env.graph_ids();
	
	if(ofile.is_open()){

		// Dataset name
		ofile<<dataset<<separator;


		// # of graphs // dont take into account the empty graph
		ofile<<env.num_graphs()-1<<separator;

		// type of edges. 0-undirected, 1-directed, 2-mixed
		ofile<<0<<separator;

		// graph names (and edge type if edgetype is 2-mixed)
		//ofile<<"#graph_names"<<"\n";
		encoded_attributes.at("graphs").emplace("name", std::map<std::string,std::string> ());
		cont=0;
		for(ged::GEDGraph::GraphID i{limits.first}; i<limits.second; i++){
			if(i == root) continue;
			encoded_attributes.at("graphs").at("name").emplace(env.get_graph_name(i), std::to_string(cont));
			ofile<<env.get_graph_name(i)<<separator;
			cont++;
		}

		// arborescence
		for(std::size_t k=0; k<arborescence.size();k++){
			ofile<<arborescence.at(k)<<separator;
		}

		// # of attributes
		// for nodes
		ofile<<alphabets.at("node_attr").size()<<separator;
		// for edges		
		ofile<<alphabets.at("edge_attr").size()<<separator;

		// Type of attributes
		for(auto const &entry : alphabets.at("node_attr")){
			ofile<< entry.first + "_type"<<separator;	
		}
		
		for(auto const &entry : alphabets.at("edge_attr")){
			ofile<< entry.first + "_type"<<separator;	
		}

		// Alphabets
		//ofile<<"#node_attr"<<"\n";		
		for(auto const attr : alphabets.at("node_attr")){
			// attr name
			ofile<<escape_char<<attr.first<<separator;
			encoded_attributes.at("node_attr").emplace(attr.first, std::map<std::string,std::string> ());
			cont=0;
			for(auto const val : attr.second){
				ofile<<val<<separator;
				encoded_attributes.at("node_attr").at(attr.first).emplace(val, std::to_string(cont));
				cont++;
			}
		}

		//ofile<<"#edge_attr"<<"\n";		
		for(auto const attr : alphabets.at("edge_attr")){
			ofile<< escape_char <<attr.first<<separator;
			encoded_attributes.at("edge_attr").emplace(attr.first, std::map<std::string,std::string> ());
			cont=0;
			for(auto const val : attr.second){
				ofile<<val<<separator;
				encoded_attributes.at("edge_attr").at(attr.first).emplace(val, std::to_string(cont));
				cont++;
			}
			
		}
		ofile.close();	
	}
	else{
		std::cout<<"Error when opening output file"<<std::endl;
		exit(1);
	}


	// code the environment into a new one
	ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> g;
	std::pair<std::pair<std::size_t, std::size_t>, ged::GXLLabel> edge;
	for(ged::GEDGraph::GraphID i{limits.first}; i<limits.second; i++){
		g = env.get_graph(i, true, true, true);
		for(std::size_t n=0; n< g.node_labels.size(); n++){
			for(auto label: g.node_labels.at(n)){
				g.node_labels.at(n).at(label.first) = encoded_attributes.at("node_attr").at(label.first).at(label.second);
			}
		}

		typename std::list<std::pair<std::pair<std::size_t, std::size_t>, ged::GXLLabel>>::iterator iter;
		
		for(iter = g.edge_list.begin(); iter != g.edge_list.end(); iter++){
			edge = (*iter);
			for(auto const l: edge.second){
				(*iter).second.at(l.first) = encoded_attributes.at("edge_attr").at(l.first).at(l.second);
			}
		}

		env_coded.load_exchange_graph(g, ged::undefined(), ged::Options::ExchangeGraphType::EDGE_LIST, env.get_graph_name(i)+"_coded", env.get_graph_class(i));

	}

}

/*
	Get the canonical node map used to decompress
*/
ged::NodeMap get_aux_node_map(std::size_t v1, std::size_t v2, vector<ged::GEDGraph::NodeID> &v_d, vector<ged::GEDGraph::NodeID> &v_i){
	if(v_d.size()>0){
		std::sort(v_d.begin(), v_d.end());
	}

	if(v_i.size()>0){
		std::sort(v_i.begin(), v_i.end());	
	}

	ged::NodeMap np = ged::NodeMap(v1, v2);
	std::size_t decr = 0;
	//std::cout<<"Start v_d for"<<std::endl;
	for(std::size_t i=0; i<v1; i++){
		if(std::find(v_d.begin(), v_d.end(), i) != v_d.end()){
			np.add_assignment(i, ged::GEDGraph::dummy_node());
			decr++;
		}
		else{
			np.add_assignment(i, i-decr);	
		}
	}
	//std::cout<<"Start v_i for"<<std::endl;
	for(std::size_t i=0; i<v_i.size(); i++){
		np.add_assignment(ged::GEDGraph::dummy_node(), i + v1);	
	}

	return np;

}


ged::NodeMap get_id_node_map(std::size_t  num_nodes, ged::NodeMap node_map, ged::NodeMap node_map_aux, ged::NodeMap node_map_id){
	std::cout<<"get_id_node_map"<<std::endl;
	std::cout<<"node_map"<<std::endl;
	std::cout<<node_map<<std::endl;
	std::cout<<"node_map_aux"<<std::endl;
	std::cout<<node_map_aux<<std::endl;
	std::cout<<"node_map_id"<<std::endl;
	std::cout<<node_map_id<<std::endl;
	ged::NodeMap res = ged::NodeMap(num_nodes,num_nodes);
	std::vector<ged::GEDGraph::GraphID> v_i;
	std::size_t cont=0;
	for(std::size_t  n = 0; n<num_nodes; n++){
		if(node_map.pre_image(n) != ged::GEDGraph::dummy_node()){
			res.add_assignment(n, node_map_aux.image(node_map_id.image(node_map.pre_image(n))));
			cont++;	
		}
		else{
			v_i.emplace_back(n);
		}
	}
	for(std::size_t  n = 0; n<v_i.size(); n++){
		res.add_assignment(v_i.at(n), cont);
		cont++;

	}
	std::cout<<"res"<<std::endl;
	std::cout<<res<<std::endl;
	return res;
}

/*

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void change_graph_version(ged::NodeMap permutation, ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &g){
	
	std::vector<UserNodeID> old_ids(g.original_node_ids);
	std::vector<UserNodeLabel> old_labels(g.node_labels);

	for(std::size_t i = 0; i < old_labels.size();i ++){
		g.node_labels.at(permutation.image(i)) = old_labels.at(i);
		g.original_node_ids.at(permutation.image(i)) = old_ids.at(i);
	}


	std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel> edge;
	std::list<std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel>> new_list;
	typename std::list<std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel>>::iterator iter;

	if(g.num_edges>0){
		for(iter = g.edge_list.begin(); iter != g.edge_list.end(); iter++){
			edge = *iter;
			*iter.first.first = permutation.image(*iter.first.first);
			*iter.first.second = permutation.image(*iter.first.second);
			//new_list.push_front(std::make_pair(std::make_pair(permutation.image(edge.first.first)
			//	permutation.image(edge.first.second)), edge.second)); 
		}
		//g.edge_list = new_list;
	}


	std::cout<<" **************** "<<std::endl;

}

*/
void permute_nodes(ged::NodeMap permutation, vector<ged::GEDGraph::NodeID> &v){
	for(std::size_t i = 0; i < v.size();i ++){
		v.at(i) = permutation.image(v.at(i));
	}
}


template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void encode_arborescence(
	std::string path, 
	std::string file_name, 
	ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env_coded,
	std::map<std::string, std::map<std::string, std::map<std::string,std::string>>> &encoded_attributes,
	std::vector<std::size_t> &arborescence,
	std::size_t root,
	bool stdout=false,
	char escape_char = '#',
	char separator = '\n'
	){

	std::ofstream ofile;
	std::ofstream ocollection;
	
	if(stdout) std::cout<<"Create children"<<std::endl;
	std::map<std::size_t, std::vector<std::size_t>> children;
	std::size_t num_nodes = arborescence.size()+1; 
	for(std::size_t i=0; i<num_nodes; i++){
		children.emplace(std::make_pair(i, std::vector<std::size_t> ()));
	}	

	for(std::size_t i=0; i<arborescence.size(); i++){
		children.at(arborescence.at(i)).emplace_back(i);			
	}


	ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> g1;
	ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> g2;

	std::vector<ged::GEDGraph::NodeID> v_d;
	std::vector<ged::GEDGraph::NodeID> v_i;
	std::vector<UserNodeLabel> varphi_i;
	std::vector<ged::GEDGraph::NodeID> v_s;
	std::vector<ged::GEDGraph::NodeID> v_is;
	std::vector<UserNodeLabel> varphi_s;
	std::vector<std::size_t> e_nd;
	std::vector<std::size_t> e_ed;
	std::vector<std::size_t> e_ni;
	std::vector<std::size_t> e_ei;
	std::vector<UserEdgeLabel> phi_ni;
	std::vector<UserEdgeLabel> phi_ei;
	std::vector<std::size_t> e_s;
	std::vector<std::size_t> e_is;
	std::vector<UserEdgeLabel> phi_s;

	std::vector<ged::GEDGraph::NodeID> v_i_aux;
	std::vector<UserNodeLabel> varphi_i_aux;
	std::map<std::size_t,std::size_t> v_i_map;

	std::string graph_file_name;
	std::vector<std::string> graph_file_names;

	ged::NodeMap node_map_aux = ged::NodeMap(1,1);
	ged::NodeMap node_map_id = ged::NodeMap(1,1);
	ged::NodeMap node_map_id_before = ged::NodeMap(1,1);
	std::map<ged::GEDGraph::GraphID, ged::NodeMap> graph_permutations;
	node_map_id_before.add_assignment(0,0);
	graph_permutations.emplace(std::make_pair(root, node_map_id_before));


	std::size_t from=0;
	std::size_t to=0;

	// Write the collection file
	
	std::pair<ged::GEDGraph::GraphID, ged::GEDGraph::GraphID> limits;
	limits = env_coded.graph_ids();
	
	graph_file_name = path + "/" + file_name + ".collection";	
	ocollection.open(graph_file_name.c_str(), ios::out);
	if(!ocollection.is_open()){
		std::cout<<"Unable to open collection file to decode"<<std::endl;
		exit(1);

	}
	if(stdout) std::cout<<"Write collection file:"<<std::endl;
	for(ged::GEDGraph::GraphID i{limits.first}; i<limits.second; i++){
		if(i==root) continue;
		graph_file_name = path + "/encoded/" + env_coded.get_graph_name(i) + "_" + std::to_string(i) + ".graph";
		graph_file_names.emplace_back(graph_file_name);
		ocollection<<graph_file_name<<separator;
		if(stdout) std::cout<<"File: "<<graph_file_name<<std::endl;
	}
	ocollection.close();


	if(stdout) std::cout<<"Start encode arborescence"<<std::endl;
	std::list<std::size_t> to_do;
	std::size_t parent_num;
	to_do.emplace_front(root);
	while(!to_do.empty()){
		
		parent_num = to_do.front();
		to_do.pop_front();
		
		for(auto const child : children.at(parent_num)){
			//std::cout<<"START FOR: "<<child<<std::endl;
		to_do.emplace_front(child);
			//if(stdout) std::cout<<"Get NodeMap"<<std::endl;
			ged::NodeMap node_map = env_coded.get_node_map(parent_num, child);
			//if(stdout) std::cout<<"Get g1"<<std::endl;
			g1 = env_coded.get_graph(parent_num, true, true, true);
			//if(stdout) std::cout<<"Get g2"<<std::endl;
			g2 = env_coded.get_graph(child, true, true, true);
			
			if(stdout) std::cout<<"Get compression sets Firts time"<<std::endl;

			get_all_compression_sets<UserNodeID, UserNodeLabel, UserEdgeLabel>(node_map,v_d,v_i,varphi_i,v_s,v_is,varphi_s,e_nd,e_ed,e_ni,e_ei,
				phi_ni,phi_ei,e_s,e_is,phi_s,g1,g2);
			if (stdout) print_compression_sets(v_d,v_i,v_s,v_is,e_nd,e_ed,e_ni,e_ei,e_s,e_is, g1, g2);


			std::cout<<std::endl<<"START-------------------------"<<std::endl;
			std::cout<<parent_num<<" -> "<<child<<std::endl;


			node_map_id_before = graph_permutations.at(parent_num);
			node_map_aux.clear();

			if(stdout) std::cout<<"Get aux node map"<<std::endl;
			permute_nodes(node_map_id_before, v_d); // Now in the compressed version
			permute_nodes(node_map_id_before, v_s); // Now in the compressed version
			permute_nodes(node_map_id_before, v_is); // Now in the compressed version

			if(stdout) std::cout<<"With modified vertex"<<std::endl;
			if (stdout) print_compression_sets(v_d,v_i,v_s,v_is,e_nd,e_ed,e_ni,e_ei,e_s,e_is, g1, g2);

			node_map_aux = get_aux_node_map(g1.num_nodes, g2.num_nodes, v_d, v_i); // Between the compressed versions	

			if(stdout) std::cout<<"Get id node map"<<std::endl;
			node_map_id = get_id_node_map(g2.num_nodes, node_map, node_map_aux, node_map_id_before);
			graph_permutations.emplace(std::make_pair(child, node_map_id));

			if(stdout) std::cout<<"Id BEFORE: "<<node_map_id_before<<std::endl;
			if(stdout) std::cout<<"Id CURRENT: "<<node_map_id<<std::endl;
			
			// Check that deletions and insertions are mutually exclusive
			if(v_d.size()>0 && v_i.size()>0) std::cout<<"ERROR: Vertex deletions and insertions at the same time"<<std::endl;
			if((e_nd.size()+e_ed.size())>0 && (e_ni.size()+e_ei.size())>0) std::cout<<"ERROR: Edge deletions and insertions at the same time"<<std::endl;


			// Now the compression sets are in terms of the compression graphs
			graph_file_name = graph_file_names.at(child);
			if(stdout) std::cout<< "File: " << graph_file_name << std::endl;

			ofile.open(graph_file_name, ios::out | ios::binary);
			if(! ofile.is_open()){
				std::cout<< "Unable to open file " << graph_file_name << std::endl;
				exit(1);				
			}
			//if(stdout) std::cout<<"Start writing to file"<<std::endl;
			
			// Arborescence is now in the info_file
			/*
			if(parent.first==root){
				// Coding a root graph
					// Use the number of the graph acording to the encoding
				ofile<<"#from:"<<"empty"<<"\n";
				
			}		
			else{
				ofile<<"#from:"<<std::to_string(parent.first)<<"\n";
			}
			*/	

			//Eq 46 and 47
			//Rec V
			//if(stdout) std::cout<<"REC V"<<std::endl;
			ofile<<g2.num_nodes<<separator;
			ofile<<v_s.size()<<separator;

			if(v_d.size()>0){
				std::sort(v_d.begin(), v_d.end());	
			}
			if(v_i.size()>0){
				std::sort(v_i.begin(), v_i.end());	
			}

			// Order the inserted nodes
			v_i_aux.clear();
			varphi_i_aux.clear();
			v_i_map.clear();
			for(std::size_t k = 0; k<v_i.size(); k++){
				v_i_aux.emplace_back(k + g1.num_nodes);
				varphi_i_aux.emplace_back(varphi_i.at(k));
				v_i_map.emplace(v_i.at(k), k + g1.num_nodes);
			}

			// Type of record. Cases in eq 46
			if(g1.num_nodes < g2.num_nodes){
								
				// Insertions
				for(std::size_t i=0; i<v_i.size(); i++){
					ofile<< escape_char <<v_i_aux.at(i)<< separator;
					for(auto const label : varphi_i_aux.at(i)){
						ofile<<label.first<<":"<<label.second<< separator;
					}
					
				}	
			}
			else{
				if(v_d.size() < v_is.size()){
					// Deletions				
					for(std::size_t i=0; i<v_d.size(); i++){
						// AQUI
						//ofile<<node_map_aux.image(v_d.at(i))<<separator;
						ofile<<v_d.at(i)<<separator; // Already transformed into tilda (node_map_id.image)
						
					}
				}
				else{
					// Identical Substitutions
					// Identical substitutions?? The image??				
					for(std::size_t i=0; i<v_is.size(); i++){
						// AQUI
						//ofile<<node_map_aux.image(v_is.at(i))<<separator;
						ofile<<v_is.at(i)<<separator; // Already transformed into tilda (node_map_id.image
						
					}					
				}
			}

			// Substitutions
				
			for(std::size_t i=0; i<v_s.size(); i++){
				// AQUI
				ofile<<escape_char<<v_s.at(i)<<separator; // Already transformed into tilda (node_map_id.image
				//ofile<<escape_char<<node_map_aux.image(v_s.at(i))<<separator;
				
				for(auto const label : varphi_s.at(i)){
					ofile<<label.first<<":"<<label.second<<separator;
				}
			}

			//Rec E

			// Get the right indices for the edges... 
			// Remember that edges are ordered using node ids. If they change, so does this ordering
			
			//std::cout<<"Edges"<<std::endl;


			// Must use the correct indices
			
			std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel> edge;
			typename std::list<std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel>>::iterator iter;

			ofile<<g2.num_edges<<separator;
			ofile<<e_s.size()<<separator;

			g1.edge_list.sort(compare_edges<UserEdgeLabel>);
			g2.edge_list.sort(compare_edges<UserEdgeLabel>);

			//std::cout<<"E_d o E_is"<<std::endl;
			if(e_ed.size() <= e_is.size()){
				ofile<<0<<separator;
				ofile<<e_ed.size()<<separator;
				for(auto e : e_ed){
					iter = g1.edge_list.begin();
					std::advance (iter,e);
					edge = (*iter);
					from = node_map_aux.image(node_map_id_before.image(edge.first.first));
					to = node_map_aux.image(node_map_id_before.image(edge.first.second));
					ofile<<escape_char<<from<<","<<to<<separator;
				}
			}
			else{
				ofile<<1<<separator;
				ofile<<e_is.size()<<"\n";
				for(auto e : e_is){
					iter = g1.edge_list.begin();
					std::advance (iter,e);
					edge = (*iter);
					from = node_map_aux.image(node_map_id_before.image(edge.first.first));
					to = node_map_aux.image(node_map_id_before.image(edge.first.second));
					ofile<<escape_char<<from<<","<<to<<separator;
				}
			}

			//std::cout<<"E_s"<<std::endl;
			// Substitutions			
			for(std::size_t i=0; i<e_s.size(); i++){
				iter = g1.edge_list.begin();
				std::advance (iter, e_s.at(i));
				edge = (*iter);
				from = node_map_aux.image(node_map_id_before.image(edge.first.first));
				to = node_map_aux.image(node_map_id_before.image(edge.first.second));
				ofile<<escape_char<<from<<","<<to<<separator;
				for(auto const label : phi_s.at(i)){
					ofile<<label.first<<":"<<label.second<<separator;
				}
				
			}


			//std::cout<<"Insertions"<<std::endl;
			// Insertions
			for(std::size_t i=0; i<e_ni.size(); i++){
				iter = g2.edge_list.begin();
				std::advance (iter, e_ni.at(i));
				edge = (*iter);	
				from = node_map_id.image(edge.first.first);
				to = node_map_id.image(edge.first.second);
				
					//from = v_i_map.at(edge.first.first);
					//to = v_i_map.at(edge.first.second);
				

				ofile<<escape_char<<from<<","<<to<<separator;
				for(auto const label : phi_ni.at(i)){
					ofile<<label.first<<":"<<label.second<<separator;
				}
				
			}

			for(std::size_t i=0; i<e_ei.size(); i++){
				// in e_ei, the edges are supposed to be inserted in existing nodes
				iter = g2.edge_list.begin();
				std::advance (iter, e_ei.at(i));
				edge = (*iter);	
				from = node_map_id.image(edge.first.first);
				to = node_map_id.image(edge.first.second);
				// AQUI
				//ofile<<escape_char<<node_map_aux.image(node_map.pre_image(edge.first.first))<<","<<node_map_aux.image(node_map.pre_image(edge.first.second))<<separator;
				//ofile<<escape_char<<node_map_id_before.image(edge.first.first)<<","<<node_map_id_before.image(edge.first.second)<<separator;
				ofile<<escape_char<<from<<","<<to<<separator;
				for(auto const label : phi_ei.at(i)){
					ofile<<label.first<<":"<<label.second<<separator;
				}
			}
			ofile.close();
			if(stdout) std::cout<<"Wrote: "<<graph_file_name<<std::endl;
		}
	}
	if(stdout) std::cout<<"End encode arborescence"<<std::endl;
}


template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void decode_collection(
	ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env,
	std::map<std::string, std::string> &args,
	char escape_char = '#'
	){


	std::string path = args.at("output_root_file") + "/" + args.at("dataset_file"); 
	std::string file_name = args.at("dataset_file"); 

	bool stdout = false;
	if(args.count("stdout")>0 && args.at("stdout") == "true") stdout = true;

	if(stdout) std::cout<<"Start DECODING"<<std::endl;

	
	// Read info_file to get the decoding structure
	std::map<std::string, std::map<std::string, std::map<std::string,std::string>>> decoded_attributes;
	decoded_attributes.emplace("graphs", std::map<std::string, std::map<std::string,std::string>> ());
	decoded_attributes.emplace("node_attr", std::map<std::string, std::map<std::string,std::string>> ());
	decoded_attributes.emplace("edge_attr", std::map<std::string, std::map<std::string,std::string>> ());
	
	std::ifstream in;
	std::string line;
	std::string attr;
	std::map<std::string,std::string> value_map;
	std::string dataset;
	std::size_t num_graphs;
	int num_node_attr;
	int num_edge_attr;
	short type_edges;
	std::vector<std::size_t> arborescence;

	int cont=0;
	int cont2=0;
	std::string info_file = path + "/" + file_name + ".info_file";
	in.open(info_file.c_str());

	if(in.is_open()){
		
		// 1. dataset name
		in>>line;
		if(stdout) std::cout<<line<<std::endl;
		dataset = line;

		// 2. # of graphs 
		in>>line;
		if(stdout) std::cout<<line<<std::endl;
		num_graphs = std::stoi(line);

		// 3. Type of edges
		in>>line;
		if(stdout) std::cout<<line<<std::endl;
		type_edges = std::stoi(line);
		type_edges++;

		// graph names
		decoded_attributes.at("graphs").emplace("name", std::map<std::string,std::string> ());
		for(std::size_t i=0; i<num_graphs; i++){
			in>>line;
			if(stdout) std::cout<<line<<std::endl;
			decoded_attributes.at("graphs").at("name").emplace(std::to_string(i), line);
		}

		// arborescence
		arborescence.clear();
		for(std::size_t i=0; i<num_graphs; i++){ // no empty graph in env
			in>>line;
			if(stdout) std::cout<<line<<std::endl;
			arborescence.emplace_back(std::stoi(line));
		}

		// # of attributes
		//nodes
		in>>line;
		if(stdout) std::cout<<line<<std::endl;
		num_node_attr = std::stoi(line);
		//edges
		in>>line;
		if(stdout) std::cout<<line<<std::endl;
		num_edge_attr = std::stoi(line);
		

		// type of attr
		for(int i=0; i<num_node_attr; i++){ // no empty graph in env
			in>>line;
			if(stdout) std::cout<<line<<std::endl;
			// Type of node attribute i
		}

		for(int i=0; i<num_edge_attr; i++){ // no empty graph in env
			in>>line;
			if(stdout) std::cout<<line<<std::endl;
			// Type of edge attribute i
		}

		// Read alphabets
		cont = -1; // number of attributes to pass to edge attributes
		cont2=0;
		value_map.clear();
		in>>line;
		if(stdout) std::cout<<line<<std::endl;
		value_map.clear();
		while(cont<num_node_attr){

			if(line.at(0)==escape_char || in.eof()){
				//add previous node if information exists
				if (!value_map.empty()){
					decoded_attributes.at("node_attr").emplace(attr, value_map);
					value_map.clear();
				}
				cont2 = 0;
				attr = line.substr(1);
				cont++;
				if(cont==num_node_attr) break;

			}
			else{
				value_map.emplace(std::to_string(cont2), line);
				cont2++;
			}

			in>>line;
			if(stdout) std::cout<<line<<std::endl;
			
		}		

		cont = -1; 
		cont2=0;
		value_map.clear();
		while(cont<num_edge_attr){
			if(line.at(0)==escape_char || in.eof()){
				//add previous node if information exists
				if (!value_map.empty()){
					decoded_attributes.at("edge_attr").emplace(attr, value_map);
					value_map.clear();
				}
				cont2 = 0;
				attr = line.substr(1);
				cont++;
				if(cont==num_edge_attr) break;

			}
			else{
				value_map.emplace(std::to_string(cont2), line);
				cont2++;
			}

			in>>line;
			if(stdout) std::cout<<line<<std::endl;
			
		}		


		in.close();
	}
	else{
		std::cout<<"Unable to open info_file. Stopping decode execution"<<std::endl;
		exit(1);
	}

	/*
	std::cout<<"Print decoding structure"<<std::endl;
	for(auto const &header: decoded_attributes){
		std::cout<<header.first<<std::endl;
		for(auto const &attr: header.second){
			std::cout<<attr.first<<std::endl;
			for(auto const &line: attr.second){
				std::cout<<line.first<<" -> "<<line.second<<std::endl;
			}
		}
	}
	*/


	// Get children structure
	std::map<std::size_t, std::vector<std::size_t>> children;
	std::size_t num_nodes = arborescence.size()+1; 
	for(std::size_t i=0; i<num_nodes; i++){
		children.emplace(std::make_pair(i, std::vector<std::size_t> ()));
	}	

	for(std::size_t i=0; i<arborescence.size(); i++){
		children.at(arborescence.at(i)).emplace_back(i);			
	}

	// Create all graphs and load them to environement

	std::size_t root = arborescence.size();
	std::string graph_file;
	std::size_t num_subs=0;

	std::size_t size_read=0;
	
	unsigned short int type_read=0;

	ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> g1;
	ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> g2;
	ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> empty_env;
	ged::GEDGraph::GraphID empty_id = empty_env.add_graph("empty","");

	ged::NodeMap node_map_aux = ged::NodeMap(1,1);
	
	std::size_t graph_id;
	std::vector<std::size_t> v_d;
	std::vector<std::size_t> v_s;
	std::vector<std::size_t> v_rest;
	std::vector<std::size_t> v_i;
	std::vector<std::size_t> v_is;
	std::vector<std::size_t> e_ed;
	std::vector<std::size_t> e_is;
	std::vector<std::size_t> e_s;
	std::vector<std::pair<std::size_t,std::size_t>> e_ed_v;
	std::vector<std::pair<std::size_t,std::size_t>> e_is_v;
	std::vector<std::pair<std::size_t,std::size_t>> e_s_v;
	//std::vector<bool> marker; // edge indices
	std::map<std::pair<std::size_t, std::size_t>, bool> marker;
	std::map<std::pair<std::size_t, std::size_t>, UserEdgeLabel> marker_labels;


	std::size_t from;
	std::size_t to;
	std::size_t pos;
	
	std::vector<UserNodeID> aux_node_ids;
	
	std::vector<UserNodeLabel> aux_node_labels;

	std::vector<UserEdgeLabel> aux_edge_labels;
	std::vector<std::pair<std::size_t, std::size_t>> aux_edges;
	std::size_t num_edges;
	std::map<std::size_t, std::size_t> map_indices;

	std::size_t parent_num;
	std::list<std::size_t> to_do;
	to_do.emplace_front(root);
	std::map<std::size_t, std::size_t> pos_to_id;

	std::map<std::string, std::string> label;
	std::string first;
	std::string second;

	std::ifstream graph_in;
	std::pair<ged::GEDGraph::GraphID, ged::GEDGraph::GraphID> limits;
	std::vector<ged::GEDGraph::GraphID> all_ids;

	typename std::list<std::pair<std::pair<std::size_t, std::size_t>,UserEdgeLabel>>::iterator iter;
	std::pair<std::pair<std::size_t, std::size_t>,UserEdgeLabel> edge;


	std::vector<std::string> graph_files;
	info_file = path + "/" + file_name + ".collection";
	in.open(info_file.c_str());

	if(in.is_open()){
		while(in>>line){
			graph_files.emplace_back(line);
		}
	}
	else{
		std::cout<<"Unable to open collection file "<<info_file<<std::endl;
	}


	while(!to_do.empty()){
		if(stdout) std::cout<<"START WHILE"<<std::endl;
		parent_num = to_do.front();
		to_do.pop_front();
		if(stdout) std::cout<<"Parent: "<<parent_num<<std::endl;
		if(stdout) std::cout<<"children: "<<children.at(parent_num).size()<<std::endl;
		if(children.at(parent_num).size()<1){
			if(stdout) std::cout<<"skip"<<std::endl;
		} 
		if(stdout) std::cout<<"star for"<<std::endl;
		for(auto const child : children.at(parent_num)){
			//std::cout<<"START FOR: "<<child<<std::endl;
			to_do.emplace_front(child);
			
			if(stdout){
				std::cout<<".............Current pos to id map: ........."<<std::endl;
				for(auto const &entry: pos_to_id){
					std::cout<<entry.first<<" -> "<<entry.second<<std::endl;
				}
				std::cout<<std::endl;		
			}
			
					
			if(stdout) std::cout<<"Parent: "<<parent_num<<", Child: "<<child<<std::endl;
			graph_file = graph_files.at(child);
			if(stdout) std::cout<<"File: "<<graph_file<<std::endl;
			graph_in.open(graph_file.c_str());
			if(graph_in.is_open()){

				v_d.clear();
				v_s.clear();
				v_is.clear();
				e_ed.clear();
				e_is.clear();
				e_s.clear();
				v_i.clear();

				e_ed_v.clear();
				e_is_v.clear();
				e_s_v.clear();

				aux_node_ids.clear();
	
				aux_node_labels.clear();

				aux_edge_labels.clear();
				aux_edges.clear();
				v_rest.clear();

				map_indices.clear();

				marker.clear();
				marker_labels.clear();



				if(parent_num == root){
					//std::cout<<"From empty graph"<<std::endl;
					g1 = empty_env.get_graph(empty_id, true, true, true);
					g2 = empty_env.get_graph(empty_id, true, true, true);
				}
				else{
					//std::cout<<"From "<<parent_num<<std::endl;
					g1 = env.get_graph(pos_to_id.at(parent_num), true, true, true);
					g2 = env.get_graph(pos_to_id.at(parent_num), true, true, true);
				}
				if(stdout){
					std::cout<<" ____________ BEGIN SOURCE GRAPH ______________"<<std::endl;
					describe_graph(g1);
					std::cout<<" ____________ END SOURCE GRAPH ______________"<<std::endl;					
				}
				
				// V
				// Line 1: # of vertices
				graph_in >> line;
				if(stdout) std::cout<<"1. "<<line<<std::endl;
				num_nodes = std::stoi(line);

				//Line 2: # subs
				graph_in >> line;
				if(stdout) std::cout<<"2. "<<line<<std::endl;
				num_subs = std::stoi(line);

				// Determine the situation							
				if(g1.num_nodes >= num_nodes){
					if(g1.num_nodes - num_nodes < num_nodes - num_subs){
						type_read = 1;	// read v_d
					}
					else{
						type_read = 2;	// read v_is
					}
				}
				else{
					type_read = 3; // insertions
				}
				if(stdout) std::cout<<"Type read: "<<type_read<<std::endl;

				cont=0;
				// Read either v_d or v_is or insertions
				switch(type_read){
					case 1:
						// Read v_d
						for (std::size_t i = 0; i < g1.num_nodes - num_nodes; i++)
						{
							graph_in >> line;
							if(stdout) std::cout<<line<<std::endl;
							v_d.emplace_back(std::stoi(line));
							v_rest.emplace_back(std::stoi(line));
						}

						break;
					case 2:
						// Read v_is
						for (std::size_t i = 0; i < num_nodes - num_subs; i++)
						{
							graph_in >> line;
							if(stdout) std::cout<<line<<std::endl;
							v_is.emplace_back(std::stoi(line));
							v_rest.emplace_back(std::stoi(line));
						} 

						break;
					case 3:
						// Read insertions
					
						if(stdout) std::cout<<"Insertions "<<num_nodes - g1.num_nodes<<std::endl;
						if(stdout) std::cout<<g1.num_nodes<<" to "<< num_nodes<<std::endl;
						
						
						for (std::size_t i = 0; i < num_nodes - g1.num_nodes; i++)
						{
							
							graph_in >> line; // index
							if(stdout) std::cout<<line<<std::endl;
							second = line.substr(1);
							v_i.emplace_back(std::stoi(second));
							map_indices.emplace(std::stoi(second), aux_node_ids.size());
							aux_node_ids.emplace_back(second);
							label.clear();
							for(int a = 0; a < num_node_attr; a++){
								graph_in >> line;	
								if(stdout) std::cout<<line<<std::endl;
								pos = line.find(":");
								first = line.substr(0,pos);
								second = line.substr(pos+1);
								label.emplace(first, second);							
							}

							aux_node_labels.emplace_back(label);
							
						}

						break;
				}


				// Read substitutions
				if(stdout) std::cout<<"Substitutions "<<num_subs<<std::endl;
				
				for (std::size_t i = 0; i < num_subs; i++)
				{

					graph_in >> line; // index
					if(stdout) std::cout<<line<<std::endl;
					second = line.substr(1);
					v_s.emplace_back(std::stoi(second));
					v_rest.emplace_back(std::stoi(second));
					map_indices.emplace(std::stoi(second), aux_node_ids.size());
					aux_node_ids.emplace_back(second);

					label.clear();

					
					for(int a = 0; a < num_node_attr; a++){
						graph_in >> line;
						if(stdout) std::cout<<line<<std::endl;
						pos = line.find(":");
						first = line.substr(0,pos);
						second = line.substr(pos+1);
						label.emplace(first, second);
					}	


					aux_node_labels.emplace_back(label);

				}

				// End V

				if(stdout) std::cout<< "Deduce v_d os v_is"<<std::endl;
				// Now deduce v_d to create the auxiliary node map
				if(type_read!=3){ // V > V'
					if(type_read==2){
						for(std::size_t n = 0; n<g1.num_nodes; n++){
							if(std::find(v_rest.begin(), v_rest.end(), n) == v_rest.end()){
								v_d.emplace_back(n);
								if(stdout) std::cout<< n << " to v_d"<<std::endl;
							}
						}
					}
				}

				if(type_read!=2){
					for(std::size_t n = 0; n<g1.num_nodes; n++){
						if(std::find(v_rest.begin(), v_rest.end(), n) == v_rest.end()){
							v_is.emplace_back(n);
							if(stdout) std::cout<< n << " to v_is"<<std::endl;
						}
					}
				}

				if(stdout) std::cout<< "Get auxiliary node map"<<std::endl;
				node_map_aux = get_aux_node_map(g1.num_nodes, num_nodes, v_d, v_i);
				if(stdout) std::cout<< node_map_aux<<std::endl;


				if(stdout) std::cout<< "Index correction"<<std::endl;
				// Correct the indices
				g2.node_labels.clear();
				g2.original_node_ids.clear();
				g2.num_nodes=num_nodes;

				label.clear();
				g2.node_labels = std::vector<UserNodeLabel>(num_nodes, label);
				g2.original_node_ids = std::vector<UserNodeID>(num_nodes, "");


				if(stdout) std::cout<< " ==== V_S ==="<<std::endl;
				for(auto const &v : v_s){
					if(stdout) std::cout<< v << " -> "<< node_map_aux.image(v) << " - > " << map_indices.at(v) <<std::endl;
					// AQUI
					g2.node_labels.at(node_map_aux.image(v)) = aux_node_labels.at(map_indices.at(v));
					//g2.node_labels.at(v) = aux_node_labels.at(map_indices.at(v));
					g2.original_node_ids.at(node_map_aux.image(v)) = aux_node_ids.at(map_indices.at(v));
					//g2.original_node_ids.at(v) = aux_node_ids.at(map_indices.at(v));
				}

				if(stdout) std::cout<< " ==== V_IS ==="<<std::endl;
				for(auto const &v : v_is){
					if(stdout) std::cout<< v << " -> "<< node_map_aux.image(v) << " - > " << v <<std::endl;
					// AQUI
					g2.node_labels.at(node_map_aux.image(v)) = g1.node_labels.at(v);
					//g2.node_labels.at(v) = g1.node_labels.at(v);
					g2.original_node_ids.at(node_map_aux.image(v)) = g1.original_node_ids.at(v);
					//g2.original_node_ids.at(v) = g1.original_node_ids.at(node_map_aux.pre_image(v));
				}				

				if(stdout) std::cout<< " ==== V_I ==="<<std::endl;
				for(auto const &v : v_i){
					if(stdout) std::cout<< v << " -> "<< map_indices.at(v)<<std::endl;
					g2.node_labels.at(v) = aux_node_labels.at(map_indices.at(v));
					g2.original_node_ids.at(v) = aux_node_ids.at(map_indices.at(v));
				}

				// Fix problem with node ids not corresponding to position in the node list
				for(std::size_t n=0; n< g2.original_node_ids.size(); n++){
					g2.original_node_ids.at(n) = std::to_string(n);
				}

				

				
				map_indices.clear();
				aux_edges.clear();
				// -------------------------------------------------------------------
				// E
				if(stdout) std::cout<<"---------edges-------------------"<<std::endl;

				// init. Get current edges and labels
				if(stdout) std::cout<<"-------- Existing edges"<<std::endl;				
				for(iter = g1.edge_list.begin(); iter != g1.edge_list.end(); iter++){
					edge = (*iter);
					if(stdout) {
						std::cout<<"Edge (in g1 tilde): "<<edge.first.first << " - "<< edge.first.second<<std::endl;
						label.clear();
						for(auto const & ll : edge.second){
							std::cout<<"\t"<<ll.first<< " = "<<ll.second<<std::endl;
							label.emplace(ll.first, ll.second);
						}
					}
					if(stdout) std::cout<<"Edge (in g2 tilde): "<<node_map_aux.image(edge.first.first) << " - "<< node_map_aux.image(edge.first.second)<<std::endl;
					//from = node_map_aux.image(edge.first.first);
					//to = node_map_aux.image(edge.first.second);
					// marker in g1 tilde
					from = edge.first.first;
					to = edge.first.second;
					marker.emplace(std::make_pair(std::make_pair(from, to), true));
					marker_labels.emplace(std::make_pair(std::make_pair(from, to), label));
				}
				if(stdout) std::cout<<"-------- end Existing edges"<<std::endl;

				if(stdout){
					std::cout<<"-------- marker"<<std::endl;
					for( auto const & mm : marker){
						std::cout<<mm.first.first<< " -- "<< mm.first.second<<" = "<<mm.second<<std::endl;
						
					}	
					std::cout<<"-------- marker_labels"<<std::endl;
					for( auto const & mm : marker_labels){
						std::cout<<mm.first.first<< " -- "<< mm.first.second<<std::endl;
						for(auto const & ll : mm.second){
							std::cout<<"\t"<<ll.first<< " = "<<ll.second<<std::endl;
						}
					}					
				}
				
				// marker = std::vector<bool> (g1.num_edges, true);

				// Line 1: # edges
				graph_in >> line;
				if(stdout) std::cout<<line<<std::endl;
				num_edges = std::stoi(line);

				// Line 2: # subs
				graph_in >> line;
				if(stdout) std::cout<<line<<std::endl;
				num_subs = std::stoi(line);

				//Line 3: type
				graph_in >> line;
				if(stdout) std::cout<<line<<std::endl;
				type_read = std::stoi(line);

				//Line 4: size of the set to read (either e_ed or e_is)
				graph_in >> line;
				if(stdout) std::cout<<line<<std::endl;
				size_read = std::stoi(line);

				

				// Read either e_ed or e_is or insertions
				if(stdout) std::cout<<"Read e_ed or e_is"<<std::endl;
				for (std::size_t i = 0; i < size_read; i++)
				{
					graph_in >> line;
					if(stdout) std::cout<<line<<std::endl;

					// Reading in g2 tilde
					pos = line.find(",");
					first = line.substr(1,pos-1);
					second = line.substr(pos+1);	
					if(stdout) std::cout<<first<<" -> "<<second<<std::endl;	
					
					if(type_read == 0){
						//e_ed.emplace_back(std::stoi(line));	// edge index						
						e_ed_v.emplace_back(std::make_pair(std::stoi(first),std::stoi(second)));						
					}
					else{
						//e_is.emplace_back(std::stoi(line));  // edge index					
						e_is_v.emplace_back(std::make_pair(std::stoi(first),std::stoi(second)));
						e_is.emplace_back(0);  // For the limits in the other loops
					}
					from = node_map_aux.pre_image(std::stoi(first));
					to = node_map_aux.pre_image(std::stoi(second));

					// Must change for directed graphs
					if(marker.count(std::make_pair(from,to))>0){
						marker.at(std::make_pair(from,to))=false;
					}
					else{
						if(marker.count(std::make_pair(to,from))>0){
							marker.at(std::make_pair(to, from))=false;
						}
						else{
							std::cout<<"Edge not found: "<<from<<" -> "<<to<<std::endl;
							exit(1);
						}
					}

				}

				// Read substitutions
				if(stdout) std::cout<<"Substitutions "<<num_subs<<std::endl;

				for (std::size_t i = 0; i < num_subs; i++)
				{		
					graph_in >> line; //index			
					if(stdout) std::cout<<line<<std::endl;
					/*
					second = line.substr(1);	
					e_s.emplace_back(std::stoi(second));
					marker.at(std::stoi(second))=false;
					iter = g1.edge_list.begin();
					std::advance (iter, std::stoi(second));
					edge = (*iter);
					from = node_map_aux.image(edge.first.first);
					to = node_map_aux.image(edge.first.second);
					*/

					pos = line.find(",");
					first = line.substr(1,pos-1);
					second = line.substr(pos+1);	
					if(stdout) std::cout<<first<<" -> "<<second<<std::endl;	
					from = std::stoi(first);
					to = std::stoi(second);
					
					aux_edges.emplace_back(std::make_pair(from, to));

					from = node_map_aux.pre_image(std::stoi(first));
					to = node_map_aux.pre_image(std::stoi(second));

					// Must change for directed graphs
					if(marker.count(std::make_pair(from,to))>0){
						marker.at(std::make_pair(from,to))=false;
					}
					else{
						if(marker.count(std::make_pair(to,from))>0){
							marker.at(std::make_pair(to, from))=false;
						}
						else{
							std::cout<<"Edge not found: "<<from<<" -> "<<to<<std::endl;
							exit(1);
						}
					}
					

					label.clear();
					for(int a = 0; a < num_edge_attr; a++){
						graph_in >> line;	
						if(stdout) std::cout<<line<<std::endl;
						pos = line.find(":");
						first = line.substr(0,pos);
						second = line.substr(pos+1);

						label.emplace(first, second);							
					}
					aux_edge_labels.emplace_back(label);

				}

				// Delete edges by node deletion
				if(stdout) std::cout<<"Deletions by node deletion "<<std::endl;
				cont=0;
				for(iter = g1.edge_list.begin(); iter != g1.edge_list.end(); iter++){
					edge = (*iter);
					for(auto const & v: v_d){
						if(edge.first.first == v || edge.first.second == v){
							if(stdout) std::cout<<"Edge deleted in g1 tilde: "<<edge.first.first << " - "<< edge.first.second<<std::endl;
							e_ed.emplace_back(cont);
							//marker.at(cont)=false;
							from = node_map_aux.image(edge.first.first);
							to = node_map_aux.image(edge.first.second);
							
							e_ed_v.emplace_back(std::make_pair(from,to));
							from = edge.first.first;
							to = edge.first.second;
							
							// Must change for directed graphs
							if(marker.count(std::make_pair(from,to))>0){
								marker.at(std::make_pair(from,to))=false;
							}
							else{
								if(marker.count(std::make_pair(to,from))>0){
									marker.at(std::make_pair(to, from))=false;
								}
								else{
									std::cout<<"Edge not found: "<<from<<" -> "<<to<<std::endl;
									exit(1);
								}
							}
							//break;
						}
					}
					cont++;
				}


				// Deduce the set that was not read
				if(stdout) std::cout<<"Deduce other sets "<<std::endl;
				cont=0;
				if(type_read == 0){
					if(stdout) std::cout<<"deduce IS "<<marker.size()<<std::endl;
					// Now deduce IS
					/*
					for(iter = g1.edge_list.begin(); iter != g1.edge_list.end(); iter++){
						if(marker.at(cont)){
							edge = (*iter);
							if(stdout) std::cout<<"Edge in IS: "<<edge.first.first << " - "<< edge.first.second<<std::endl;
							e_is.emplace_back(cont);
							//map_indices.emplace_back(cont, aux_edges.size());
							from = node_map_aux.image(edge.first.first);
							to = node_map_aux.image(edge.first.second);
							aux_edges.emplace_back(std::make_pair(from, to));
							aux_edge_labels.emplace_back(edge.second);
							marker.at(cont) = false;
						}
						cont++;
					}	
					*/

					for(const auto & m : marker){
						
						if(m.second){
							if(stdout) std::cout<< "Edge in IS in g1 tilde: "<<m.first.first << " - " << m.first.second << std::endl;
							
							from = node_map_aux.image(m.first.first);
							to = node_map_aux.image(m.first.second);

							if(stdout) std::cout<<" - after node_map_aux: " << from << " - " << to << std::endl;	

							marker.at(std::make_pair(m.first.first, m.first.second))=false;
							
							e_is_v.emplace_back(std::make_pair(from,to));
							aux_edges.emplace_back(std::make_pair(from,to));

							
							aux_edge_labels.emplace_back(marker_labels.at(m.first));

							e_is.emplace_back(cont); // For the limits in the other loops
						}
					}

					if(stdout) std::cout<<"END ______ deduce IS "<<std::endl;



				}
				else{
					if(stdout) std::cout<<"deduce ED "<<std::endl;
					// Now deduce ED
					/*
					for(iter = g1.edge_list.begin(); iter != g1.edge_list.end(); iter++){
						if(marker.at(cont)){
							edge = (*iter);
							if(stdout) std::cout<<"Edge deleted: "<<edge.first.first << " - "<< edge.first.second<<std::endl;
							e_ed.emplace_back(cont);
							marker.at(cont) = false;
						}
						cont++;
					}	
					*/

					for(const auto & m : marker){
						
						if(m.second){
							if(stdout) std::cout<<"Edge in ED: "<<m.first.first << " - "<< m.first.second<<std::endl;
							from = node_map_aux.image(m.first.first);
							to = node_map_aux.image(m.first.second);
							marker.at(std::make_pair(m.first.first,m.first.second))=false;
							e_ed_v.emplace_back(std::make_pair(from,to));
							
						}
					}

					if(stdout) std::cout<<"END ______ deduce ED "<<std::endl;
					
				}

				// Read insertions
				if(num_edges > num_subs + e_is.size()){ 
					if(stdout) std::cout<<"Insertions: "<< num_edges - num_subs - e_is.size() <<std::endl;
					
					for (std::size_t i = 0; i < num_edges - num_subs - e_is.size(); i++)
					{		
						graph_in >> line; // index			
						if(stdout) std::cout<<line<<std::endl;
						pos = line.find(",");
						first = line.substr(1,pos-1);
						second = line.substr(pos+1);	
						if(stdout) std::cout<<first<<" -> "<<second<<std::endl;
						//map_indices.emplace_back(std::stoi(second), aux_edges.size());
						aux_edges.emplace_back(std::make_pair(std::stoi(first),std::stoi(second)));
						

						label.clear();
						
						for(int a = 0; a < num_edge_attr; a++){
							graph_in >> line;	
							if(stdout) std::cout<<line<<std::endl;
							pos = line.find(":");
							first = line.substr(0,pos);
							second = line.substr(pos+1);
							if(stdout) std::cout<<"Edge label: "<<first << " - "<< second <<"."<<std::endl;
							label.emplace(first, second);						
						}

						aux_edge_labels.emplace_back(label);

					}

				}
				

				g2.edge_list.clear();
				if (stdout) std::cout<<"BUILD EDGES: "<<aux_edges.size()<<", "<<aux_edge_labels.size()<<std::endl;
				for(std::size_t i =0; i<aux_edges.size(); i++){
					g2.edge_list.push_front(std::make_pair(aux_edges.at(i), aux_edge_labels.at(i)));
				}

				g2.num_edges = g2.edge_list.size();
				
				// End E
				if(stdout){
					std::cout<<" ____________ BEGIN RESULTING GRAPH ______________"<<std::endl;
					describe_graph(g2);
					std::cout<<" ____________ END RESULTING GRAPH ______________"<<std::endl;					
				}

				graph_id = env.load_exchange_graph(g2, ged::undefined(), ged::Options::ExchangeGraphType::EDGE_LIST, decoded_attributes.at("graphs").at("name").at(std::to_string(child)), "class");
				if (stdout) std::cout<<"Done: "<<env.num_graphs()<<std::endl;
				all_ids.emplace_back(graph_id);
				pos_to_id.emplace(child, graph_id);
				graph_in.close();
				//std::cout<<"End "<<std::endl;
			}
			else{
				std::cout<<"Unable to open graph file "<<graph_file<<std::endl;
				exit(1);
			}
		//std::cout<<"END FOR"<<std::endl;
		}

	if(stdout) std::cout<<"END WHILE"<<std::endl;
	}

	// Decode environement
	
	if(stdout) 
	{
		std::cout<<"Decode environement"<<std::endl;
	
		std::cout<<"Print decoding structure"<<std::endl;
		for(auto const &header: decoded_attributes){
			std::cout<<header.first<<std::endl;
			for(auto const &attr: header.second){
				std::cout<<attr.first<<std::endl;
				for(auto const &l: attr.second){
					std::cout<<l.first<<" -> "<<l.second<<std::endl;
				}
			}
		}
	}
	
	ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> env_decoded;
	
	if(stdout) std::cout<<"Num graphs: "<<env.num_graphs()<<std::endl;
	ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> g;
	std::pair<std::pair<std::size_t, std::size_t>, ged::GXLLabel> edge_gxl;
	limits = env.graph_ids();



	for(ged::GEDGraph::GraphID i = limits.first; i<limits.second; i++){
		if(i == root) continue;
		if(stdout) std::cout<<"Graph: "<<i<<std::endl;
		g = env.get_graph(i, true, true, true);
		
		if(stdout){
			std::cout<<"describe .............. "<<std::endl;
			describe_graph(g);
			std::cout<<"end ............. "<<std::endl;
		}
		for(std::size_t n=0; n< g.node_labels.size(); n++){
			for(auto const l: g.node_labels.at(n)){
				if(stdout) std::cout<<l.first<<" -> "<<l.second<<std::endl;
				if(stdout) std::cout<<decoded_attributes.at("node_attr").at(l.first).at(l.second)<<std::endl;
				g.node_labels.at(n).at(l.first) = decoded_attributes.at("node_attr").at(l.first).at(l.second);
			}
		}

		typename std::list<std::pair<std::pair<std::size_t, std::size_t>, ged::GXLLabel>>::iterator iter_gxl;
		
		for(iter_gxl = g.edge_list.begin(); iter_gxl != g.edge_list.end(); iter_gxl++){
			edge_gxl = (*iter_gxl);
			for(auto const l: edge_gxl.second){
				if(stdout) std::cout<<l.first<<" -> "<<l.second<<std::endl;
				if(stdout) std::cout<<decoded_attributes.at("edge_attr").at(l.first).at(l.second)<<std::endl;
				(*iter_gxl).second.at(l.first) = decoded_attributes.at("edge_attr").at(l.first).at(l.second);
			}
		}

		env_decoded.load_exchange_graph(g, ged::undefined(), ged::Options::ExchangeGraphType::EDGE_LIST, env.get_graph_name(i), env.get_graph_class(i));

		
	}


	//end
	env = env_decoded;
	if(stdout) std::cout<<"END DECODE"<<std::endl;

}



//TODO: add template
void 
get_compression_data(
	std::vector<std::string> &headers,
	std::vector<std::string> &values,
	std::map<std::string, std::string> &args
	){

	bool stdout=false;
	if(args.count("stdout")>0){
		if(args.at("stdout")=="true"){
			stdout = true;
		}
	}

	if(stdout) std::cout<<"--------------START and INPUTS-----------------\n";
	ged::Seconds runtime;
	auto start = std::chrono::high_resolution_clock::now();

	std::string collection_file;
	std::string graph_dir;

	if(args.count("collection_file")>0){
		collection_file = args.at("collection_file");
	}
	else{
		std::cout<<"No field collection_file in args. Stopping execution"<<std::endl;
		return;
	}

	if(args.count("graph_dir_file")>0){
		graph_dir = args.at("graph_dir_file");
	}
	else{
		std::cout<<"No field graph_dir in args. Stopping execution"<<std::endl;
		return;
	}

	headers.emplace_back("collection_file");
	values.emplace_back(collection_file);

	headers.emplace_back("graph_dir");
	values.emplace_back(graph_dir);

	if(stdout) std::cout<<"Collection file: "<<collection_file<<std::endl;
	if(stdout) std::cout<<"Graph directory: "<<graph_dir<<std::endl;


	if(stdout) std::cout<<"--------------LOAD GRAPHS, GET GRAPH STRUCTURE-----------------"<<std::endl;	

	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;

	std::vector<ged::GEDGraph::GraphID> graph_ids(env.load_gxl_graphs(graph_dir, collection_file,
			ged::Options::GXLNodeEdgeType::LABELED, ged::Options::GXLNodeEdgeType::LABELED));

	// Add empty graph at the end
	ged::GEDGraph::GraphID empty_id = env.add_graph("empty","");

	if(stdout) std::cout<<"Number of graphs: "<<env.num_graphs()<<std::endl;
	std::map<std::string, std::map<std::string, std::vector<std::string>>> distribution;
	std::map<std::string, std::map<std::string, std::set<std::string>>> alphabets;
	

	double b_ni;
	double b_na;
	double b_ei; 
	double b_ea;

	get_graphs_structure<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>(env, distribution, alphabets, b_ni, b_na, b_ei, b_ea);

	if(stdout){
		std::cout<<"Structure:"<<std::endl;
		std::cout<<"Nodes: "<<std::endl;
		for(auto const& a: alphabets.at("node_attr")){
			std::cout<<"\t"<<a.first<<": "<<a.second.size()<<" values"<<std::endl;	
		}
		std::cout<<"Edges: "<<std::endl;
		for(auto const& a: alphabets.at("edge_attr")){
			std::cout<<"\t"<<a.first<<": "<<a.second.size()<<" values"<<std::endl;	
		}
		
		std::cout<<"b_ni: "<<b_ni<<", b_na: "<<b_na<<std::endl;
		std::cout<<"b_ei: "<<b_ei<<", b_ea: "<<b_ea<<std::endl;

	}


	std::string output_root_path = "";
	std::string dataset_name = "";
	std::map<std::string, std::map<std::string, std::map<std::string,std::string>>> encoded_attributes;
	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env_coded;
	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env_uncoded;
	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env_dummy;

	if(args.count("encode")>0){
		if(args.at("encode")=="true"){
			
			if(stdout) std::cout<<"--------------ENCODE COLLECTION-----------------"<<std::endl;


			if(args.count("output_root_file")>0){
				output_root_path = args.at("output_root_file");
			}
			else{
				std::cout<<"No field output_root_file in args. Stopping execution"<<std::endl;
				return;
			}

			if(args.count("dataset_file")>0){
				dataset_name = args.at("dataset_file");
			}
			else{
				std::cout<<"No field dataset_file in args. Stopping execution"<<std::endl;
				return;
			}

			if(stdout){
				std::cout<<"############   Collection BEFORE encoding   ################"<<std::endl;
				for(auto i : graph_ids){
					std::cout<<std::endl;
					std::cout<<env.get_graph_name(i)<<std::endl;
					describe_graph(env.get_graph(i, true, true, true));
					std::cout<<std::endl;
				}
				std::cout<<"############   Collection BEFORE encoding - end   ################"<<std::endl;
			}

			std::vector<std::size_t> dummy;
			// First encoding just to generate the coded environment. 
			// In the second encoding we will create the actual files
			encode_environment<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>(
			env_coded, encoded_attributes, output_root_path + "/" + dataset_name, dataset_name, env, alphabets, dummy, empty_id);

			env_uncoded = env;
			env = env_coded;

			if(stdout && env.num_graphs()<10){
				std::cout<<"Uncoded collection"<<std::endl;
				for(auto i : graph_ids){
					std::cout<<env_uncoded.get_graph_name(i)<<std::endl;
					describe_graph(env_uncoded.get_graph(i, true, true, true));
				}

				std::cout<<"Coded collection"<<std::endl;
				for(auto i : graph_ids){
					std::cout<<env.get_graph_name(i)<<std::endl;
					describe_graph(env.get_graph(i, true, true, true));
				}
			}
			
		}
	}


	if(stdout) std::cout<<"--------------SET COMPRESSION EDIT COST-----------------"<<std::endl;

	headers.emplace_back("b_ni");
	values.emplace_back(std::to_string(b_ni));

	headers.emplace_back("b_na");
	values.emplace_back(std::to_string(b_na));

	headers.emplace_back("b_ei");
	values.emplace_back(std::to_string(b_ei));

	headers.emplace_back("b_ea");
	values.emplace_back(std::to_string(b_ea));

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

	if(stdout) std::cout<<"--------------ENV INIT-----------------"<<std::endl;
	env.set_edit_costs(ged::Options::EditCosts::COMPRESSION, comp_costs);
	env.init(ged::Options::InitType::LAZY_WITHOUT_SHUFFLED_COPIES);


	// Set method
	std::string ged_method_options = "";
	if(args.count("ged_method_options")>0){
		ged_method_options = args.at("ged_method_options");
	}

	if(args.count("ged_method")>0){
		if(stdout) std::cout<<"GED method: "<<args.at("ged_method")<<std::endl;
		if(args.at("ged_method") == "branch_uniform"){
			env.set_method(ged::Options::GEDMethod::BRANCH_UNIFORM, "--threads " + ged_method_options);
		}
		else{
			if(args.at("ged_method") == "branch_fast"){
				env.set_method(ged::Options::GEDMethod::BRANCH_FAST, "--threads " + ged_method_options);
			}
			else{
				if(args.at("ged_method") == "ipfp"){
					env.set_method(ged::Options::GEDMethod::IPFP, "--threads " + ged_method_options);
				}
				else{
					std::cout<<"No valid ged_method in args. Setting it to BRANCH_UNIFORM"<<std::endl;
					env.set_method(ged::Options::GEDMethod::BRANCH_UNIFORM, "--threads " + ged_method_options);
				}
			}
		}
	}
	else{
		std::cout<<"No field ged_method in args. Setting it to BRANCH_UNIFORM"<<std::endl;
		env.set_method(ged::Options::GEDMethod::BRANCH_UNIFORM, "--threads " + ged_method_options);
		
	}
	
	

	if(stdout) std::cout<<"--------------RUN METHOD-----------------"<<std::endl;
	std::pair<ged::GEDGraph::GraphID, ged::GEDGraph::GraphID> limits = env.graph_ids();

	if(stdout) std::cout<<"Limits (no of graphs): "<<limits.first<<", "<<limits.second-1<<std::endl;
	if(stdout) std::cout<<"Empty graph id: "<<empty_id<<std::endl;
	
	std::vector<std::vector<double>> upper_bounds;
	std::vector<std::vector<double>> upper_bounds_refined;
	std::vector<double> aux_line;
	ged::GEDGraph::GraphID i_par;
	ged::GEDGraph::GraphID j_par;
	ged::GEDGraph::GraphID k_par;

	// Compute GED upper bounds
	for(k_par = limits.first; k_par<(limits.second)*(limits.second); k_par++){
		env.run_method(k_par/limits.second,k_par % limits.second);
	}


	for(i_par = limits.first; i_par<limits.second; i_par++){
		aux_line.clear();
		for(j_par = limits.first; j_par<limits.second; j_par++){
			if(j_par!=empty_id){
				if(i_par==j_par){
					aux_line.emplace_back(std::numeric_limits<std::size_t>::max());
				}
				else{
					// Correction factor of 3*b_ni+2*b_ei
					aux_line.emplace_back(env.get_upper_bound(i_par,j_par)+ 3*b_ni+2*b_ei);					
				}				
			}			
		}
		upper_bounds.emplace_back(aux_line);
		upper_bounds_refined.emplace_back(aux_line);
	}


	if(stdout && upper_bounds.size()<10){
		std::cout<<"Cost matrix: "<<upper_bounds.size()<<" x "<<upper_bounds.at(0).size()<<std::endl;
		for(std::size_t i=0; i<upper_bounds.size(); i++){
			for(std::size_t j=0; j<upper_bounds.at(i).size(); j++){
				std::cout<<upper_bounds.at(i).at(j)<<", ";
			}
			std::cout<<std::endl;
		}
	}

	std::vector<std::pair<std::size_t,std::size_t>>  collection_graph;
	for(ged::GEDGraph::GraphID i{limits.first}; i<limits.second; i++){
		for(ged::GEDGraph::GraphID j{limits.first}; j<limits.second; j++){
			if(i==j || j==empty_id){
				continue;
			}
			else{
				collection_graph.emplace_back(std::make_pair(i,j));
			}

		}
	}
	auto start_arb = std::chrono::high_resolution_clock::now();
	runtime = start_arb - start;
	double gedlib_time = runtime.count();
	

	// Free up memory before passing to the next section
	// Just before, get the "total" cost for reference
	double base_compression_cost = calculate_total_compression_cost<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>(env, empty_id, b_ni, b_na, b_ei, b_ea);
	//env = ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>();


	if(stdout) std::cout<<"--------------SPANNING ARBORESCENCE OF MINIMUM WEIGHT-----------------"<<std::endl;
	std::size_t root = upper_bounds.size()-1; // empty id
	if(stdout) std::cout<<"Root: "<<root<<std::endl;
	

	std::vector<std::size_t> arborescence;
	double cost_arb = 0;

	spanning_arborescence_of_minimum_weight(arborescence, cost_arb, collection_graph, upper_bounds,root, false);
	
	if(stdout) std::cout<<"Cost of arborescence: "<<cost_arb<<std::endl;
	if(stdout){
		std::cout<<"Arborescence: "<<std::endl;
		for(std::size_t i =0; i<arborescence.size(); i++){
			std::cout<<arborescence.at(i)<<", ";
		}
		std::cout<<std::endl;
	} 

	std::map<std::size_t, std::vector<std::size_t>> depth_degrees;

	// Manually inserted values
	headers.emplace_back("cost_arborescence");
	values.emplace_back(std::to_string(cost_arb));

	headers.emplace_back("base_compression_cost");
	values.emplace_back(std::to_string(base_compression_cost));

	get_arborescence_info(headers, values, depth_degrees, arborescence, root);

	auto start_ref = std::chrono::high_resolution_clock::now();
	runtime = start_ref - start_arb;
	double arb_time = runtime.count();
	

	if(stdout) std::cout<<"--------------REFINEMENT-----------------"<<std::endl;


	// Set method
	std::string ged_method_refinement_options = "";
	if(args.count("ged_method_refinement_options")>0){
		ged_method_refinement_options = args.at("ged_method_refinement_options");
	}

	if(args.count("ged_method_refinement")>0){
		if(stdout) std::cout<<"GED method (refinement): "<<args.at("ged_method_refinement")<<std::endl;
		if(args.at("ged_method_refinement") == "branch_uniform"){
			env.set_method(ged::Options::GEDMethod::BRANCH_UNIFORM, "--threads " + ged_method_refinement_options);
		}
		else{
			if(args.at("ged_method_refinement") == "branch_fast"){
				env.set_method(ged::Options::GEDMethod::BRANCH_FAST, "--threads " + ged_method_refinement_options);
			}
			else{
				if(args.at("ged_method_refinement") == "branch_tight"){
					env.set_method(ged::Options::GEDMethod::BRANCH_TIGHT, "--threads " + ged_method_refinement_options);
				}
				else{
					if(args.at("ged_method_refinement") == "ipfp"){
						env.set_method(ged::Options::GEDMethod::IPFP, "--threads " + ged_method_refinement_options + " --initialization-method BRANCH_UNIFORM");
					}
					else{
						if(args.at("ged_method_refinement") == "ring"){

							std::string ring_train_path="";
							if(args.count("ring_train_path")>0){
								ring_train_path = args.at("ged_method_train_path");
							}

							std::string ring_led_method="";
							if(args.count("ring_led_method")>0){
								ring_led_method = args.at("ring_led_method");
							}

							std::string ring_suffix="";
							if(args.count("ring_suffix")>0){
								ring_suffix = args.at("ring_suffix");
							}

							// Train data must be in the ring_train_path folder. Its name is ring_method_suffix
							env.set_method(ged::Options::GEDMethod::RING, init_options(ring_train_path, ring_led_method + "_" + ring_suffix, false, true, std::stoi(ged_method_refinement_options)) +  " --led-method " + ring_led_method);
							env.init_method();

						}
						else{
							std::cout<<"No valid ged_method in args. Setting it to BRANCH_UNIFORM"<<std::endl;
							env.set_method(ged::Options::GEDMethod::BRANCH_UNIFORM, "--threads " + ged_method_refinement_options);	
						}
					}
				}
			}
		}
	}
	else{
		std::cout<<"No field ged_refine_method in args. Setting it to BRANCH_UNIFORM"<<std::endl;
		env.set_method(ged::Options::GEDMethod::BRANCH_UNIFORM, "--threads " + ged_method_refinement_options);
		
	}

	//env.set_edit_costs(ged::Options::EditCosts::COMPRESSION, comp_costs);
	//env.init(ged::Options::InitType::LAZY_WITHOUT_SHUFFLED_COPIES);


	std::size_t refinement_size = 0;
	if(args.count("refinement_size")>0){
		refinement_size = std::stoi(args.at("refinement_size"));
	}
	else{
		std::cout<<"No field refinement_size in args. Setting it to 0"<<std::endl;
	}	

	if(stdout) std::cout<<"Refinement size: "<<refinement_size<<std::endl;

	std::size_t step;
	std::size_t node;
	for(std::size_t n =0; n<arborescence.size(); n++){
		step=0;
		node = n;
		while(step<refinement_size && node!=root){
			step++;
			env.run_method(arborescence.at(node), node);
			upper_bounds_refined.at(arborescence.at(node)).at(node) = min(env.get_upper_bound(arborescence.at(node), node) + 3*b_ni+2*b_ei, upper_bounds.at(arborescence.at(node)).at(node));			
			node = arborescence.at(node);
		}
	}


	if(stdout && upper_bounds_refined.size()<15){
		std::cout<<"Cost matrix (refined): "<<upper_bounds_refined.size()<<" x "<<upper_bounds_refined.at(0).size()<<std::endl;
		for(std::size_t i=0; i<upper_bounds_refined.size(); i++){
			for(std::size_t j=0; j<upper_bounds_refined.at(i).size(); j++){
				std::cout<<upper_bounds_refined.at(i).at(j)<<", ";
			}
			std::cout<<std::endl;
		}
	}

	auto end_ref = std::chrono::high_resolution_clock::now();
	runtime = end_ref - start_ref;
	double ref_time = runtime.count();


	std::vector<std::size_t> arborescence_ref;
	double cost_arb_ref = 0;

	spanning_arborescence_of_minimum_weight(arborescence_ref, cost_arb_ref, collection_graph, upper_bounds_refined,root, false);

	if(stdout) std::cout<<"Cost of arborescence (refined): "<<cost_arb_ref<<std::endl;
	if(stdout && arborescence_ref.size()<15){
		std::cout<<"Arborescence (refined): "<<std::endl;
		for(std::size_t i =0; i<arborescence_ref.size(); i++){
			std::cout<<arborescence_ref.at(i)<<", ";
		}
		std::cout<<std::endl;
	} 

	// Manually inserted values
	headers.emplace_back("cost_arborescence_refined");
	values.emplace_back(std::to_string(cost_arb_ref));

	get_arborescence_info(headers, values, depth_degrees, arborescence_ref, root);

	auto end_ref_arb = std::chrono::high_resolution_clock::now();
	runtime = end_ref_arb - end_ref;
	double ref_arb_time = runtime.count();


	if(stdout) std::cout<<"--------------TO COMPARE-----------------"<<std::endl;

	double sum=0;
	for(ged::GEDGraph::GraphID i{limits.first}; i<limits.second; i++){
		if(i != empty_id) sum+=upper_bounds.back().at(i);
	}
	
	headers.emplace_back("base_compression_cost_2");
	values.emplace_back(std::to_string(sum));

	double compression_ratio = cost_arb/base_compression_cost;
	
	headers.emplace_back("compression_ratio");
	values.emplace_back(std::to_string(compression_ratio));

	double compression_ratio_ref = cost_arb_ref/base_compression_cost;
	
	headers.emplace_back("compression_ratio_refined");
	values.emplace_back(std::to_string(compression_ratio_ref));

	if(stdout) std::cout<<"Total cost 1: "<<base_compression_cost<<std::endl;
	if(stdout) std::cout<<"Total cost 2 (GED from empty graph): "<<sum<<std::endl;
	if(stdout) std::cout<<"Compression ratio: "<<compression_ratio<<std::endl;
	if(stdout) std::cout<<"Compression ratio (refined): "<<compression_ratio_ref<<std::endl;
	
	headers.emplace_back("gedlib_runtime_initial");
	values.emplace_back(std::to_string(gedlib_time));

	headers.emplace_back("gedlib_runtime_refinement");
	values.emplace_back(std::to_string(ref_time));

	if(stdout) std::cout<<"GEDLIB time (initial): "<<gedlib_time<<std::endl;
	if(stdout) std::cout<<"GEDLIB time (refinement): "<<ref_time<<std::endl;
	
	headers.emplace_back("spanning_arb_runtime");
	values.emplace_back(std::to_string(arb_time));
	if(stdout) std::cout<<"Spanning arborescence time: "<<arb_time<<std::endl;

	headers.emplace_back("refine_arb_runtime");
	values.emplace_back(std::to_string(ref_arb_time));
	if(stdout) std::cout<<"Spanning arborescence time (refinement): "<<ref_arb_time<<std::endl;



	if(args.count("encode")>0){
		if(args.at("encode")=="true"){
			
			if(stdout) std::cout<<"--------------ENCODE COLLECTION-----------------"<<std::endl;

			// Re encode to get the arborescence right
			encode_environment<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>(
			env_dummy, encoded_attributes, output_root_path + "/" + dataset_name, dataset_name, env_uncoded, alphabets, arborescence_ref, empty_id);
			
			encode_arborescence<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>(output_root_path + "/" + dataset_name, dataset_name,
			 env, encoded_attributes, arborescence_ref, root, stdout);
		}
	}
}

void treat_dataset(std::map<std::string, std::string> args){

	std::vector<std::string> headers;
	std::vector<std::string> values;
	std::vector<std::string> gxl_file_names;
	std::vector<std::string> graph_classes;

	std::cout<<"**********    GET COMPRESSION DATA AND ENCODE   ************"<<std::endl;
	get_compression_data(headers, values, args);	


	std::cout<<"************    WRITE RESULTS FILE    **************"<<std::endl;

	std::ofstream output_file;
	output_file.open(args.at("output_results_file").c_str(), ios::out | ios::app);
	
	if(output_file.is_open()){
		if(args.at("first_iteration")=="true"){
			write_to_file(output_file, headers);			
		}
		write_to_file(output_file, values);		
		output_file.close();	
	}
	else{
		std::cout<<"Error when opening output file"<<std::endl;
		exit(1);
	}


	if(args.count("decode") && args.at("decode")=="true"){

		std::cout<<"**********    DECODE    ************"<<std::endl;	

		ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env_decoded;
		decode_collection<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>(env_decoded, args);	

		std::pair<ged::GEDGraph::GraphID, ged::GEDGraph::GraphID> graph_ids;
		graph_ids = env_decoded.graph_ids();

		if(args.count("stdout")>0 && args.at("stdout")=="true"){
			std::cout<<"############   Collection AFTER decoding   ################"<<std::endl;
			for(std::size_t i=graph_ids.first; i<graph_ids.second; i++){
				std::cout<<std::endl;
				std::cout<<env_decoded.get_graph_name(i)<<std::endl;
				describe_graph(env_decoded.get_graph(i, true, true, true));
				std::cout<<std::endl;
			}

			std::cout<<"############   Collection AFTER decoding - end   ################"<<std::endl;
		}

		
		if(args.count("write_decoded")>0 && args.at("write_decoded")=="true"){
			gxl_file_names.clear();
			graph_classes.clear();
			for(std::size_t i=graph_ids.first; i<graph_ids.second; i++){
				gxl_file_names.emplace_back(env_decoded.get_graph_name(i));
				graph_classes.emplace_back(env_decoded.get_graph_class(i));
				env_decoded.save_as_gxl_graph(i,  args.at("output_root_file") + "/" + args.at("dataset_file") + "/decoded/" + env_decoded.get_graph_name(i));
				
			}	
			env_decoded.save_graph_collection(args.at("output_root_file") + "/" + args.at("dataset_file") + "/decoded/" + args.at("dataset_file") + ".xml",  gxl_file_names,  graph_classes);
		}

	}
	
}


int main(int argc, char* argv[]){

	std::vector<std::string> headers;
	std::vector<std::string> values;
	std::map<std::string, std::string> args;
	
	// LLenar args y lanzar

	
	args.emplace(std::make_pair("stdout",argv[1]));

	args.emplace(std::make_pair("collection_file",argv[2]));
	args.emplace(std::make_pair("graph_dir_file",argv[3]));
	args.emplace(std::make_pair("output_root_file",argv[4]));
	args.emplace(std::make_pair("dataset_file",argv[5]));

	args.emplace(std::make_pair("output_results_file",argv[6]));

	args.emplace(std::make_pair("ged_method",argv[7]));
	args.emplace(std::make_pair("ged_method_options",argv[8]));
	args.emplace(std::make_pair("ged_method_refinement",argv[9]));
	args.emplace(std::make_pair("ged_method_refinement_options",argv[10]));
	args.emplace(std::make_pair("refinement_size",argv[11]));
	args.emplace(std::make_pair("encode",argv[12]));
	args.emplace(std::make_pair("decode",argv[13]));
	args.emplace(std::make_pair("write_decoded",argv[14]));

	args.emplace(std::make_pair("train_set",argv[15]));
	args.emplace(std::make_pair("train_path",argv[16]));
	args.emplace(std::make_pair("ring_method",argv[17]));

	


	std::ifstream in_file_collections(args.at("collection_file").c_str());
	std::ifstream in_file_graphs(args.at("graph_dir_file").c_str());
	std::ifstream in_file_output_root(args.at("output_root_file").c_str());
	std::ifstream in_file_dataset(args.at("dataset_file").c_str());
	


    if (!in_file_collections || !in_file_graphs || !in_file_output_root || !in_file_dataset) {
        std::cout << "Unable to open files";
        exit(1); // terminate with error
    }

    std::string input_collection_file;
    std::string input_graph_dir;
    std::string input_output_root;
    std::string input_dataset;

    args.emplace(std::make_pair("first_iteration", "true"));

    while (in_file_collections >> input_collection_file) {
        in_file_graphs >> input_graph_dir;
        in_file_output_root >> input_output_root;
        in_file_dataset >> input_dataset;


        args.at("collection_file") = input_collection_file;
        args.at("graph_dir_file") = input_graph_dir;
        args.at("output_root_file") = input_output_root;
        args.at("dataset_file") = input_dataset;

		treat_dataset(args); 
		args.at("first_iteration") = "false";
		
    }
   
	return 0;
	
}


