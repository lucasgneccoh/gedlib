#ifndef COMPRESS_SRC_ABC_HPP_
#define COMPRESS_SRC_ABC_HPP_

#include <iostream>
#include <random>
#include <string>
#include <cmath>
#include <algorithm>
#include <limits>
#include <climits> // for CHAR_BITS
#include <list>
#include "../ext/frangio68/Minimal-Spanning-Arborescence-master/MSArbor.C"

#define COMPRESS_EDIT_COST
#include "../../src/env/ged_env.hpp"
#undef COMPRESS_EDIT_COST


/*!
 * @namespace ged
 * @brief Global namespace for GEDLIB.
 */
namespace ged {

class GED_ABC {
public:

	/*!
	 * @brief Destructor.
	 */
	~GED_ABC();

	/*!
	 * @brief Constructor.
	 */
	GED_ABC();



	void set_omega(double d);

	/*!
	 * @brief Adds a new uninitialized graph to the environment. Call init() after calling this method.
	 * @param[in] graph_name The name of the added graph. Empty if not specified.
	 * @param[in] graph_class The class of the added graph. Empty if not specified.
	 * @return The ID of the newly added graph.
	 */
	template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
	void describe_graph(ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> g);


	template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
	void get_graphs_structure(ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env,
		std::map<std::string, std::map<std::string, std::vector<std::string>>> &distribution,
		std::map<std::string, std::map<std::string, std::set<std::string>>> &alphabets,
		std::map<std::string, std::map<std::string, std::size_t>> &attr_sizes,
		std::size_t &b_ni, std::size_t &b_na, std::size_t &b_ei, std::size_t &b_ea, 
		bool &fast_node_translate, bool &fast_edge_translate);

	
	std::pair<double, double> simplified_compression_size_from_empty_graph(
		std::size_t num_nodes, std::size_t num_edges, 
		std::size_t &b_ni, std::size_t &b_na,
		std::size_t &b_ei, std::size_t &b_ea);

	std::size_t compute_cost_edge_pairs(std::size_t V, std::size_t E, std::size_t b_ni, std::size_t b_na, std::size_t b_ei, std::size_t b_ea);

	std::size_t compute_cost_triangular_matrix(std::size_t V, std::size_t  E, std::size_t b_ni, std::size_t b_na, std::size_t b_ei, std::size_t b_ea);

	std::size_t compute_cost_abc(std::size_t V, std::size_t  E, std::size_t b_ni, std::size_t b_na, std::size_t b_ei, std::size_t b_ea);

	template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
	std::size_t base_compr_cost_edge_pairs(ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env, std::size_t &b_ni, std::size_t &b_na, std::size_t &b_ei, std::size_t &b_ea, std::size_t constant=0);


	template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
	std::size_t base_compr_cost_triangular_matrix(ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env, std::size_t &b_ni, std::size_t &b_na, std::size_t &b_ei, std::size_t &b_ea, std::size_t constant=0);

	template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
	std::size_t base_compr_cost_abc(ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env, std::size_t &b_ni, std::size_t &b_na, std::size_t &b_ei, std::size_t &b_ea, std::size_t constant=0);

	/*
		Args can have the following fields:
		edit_cost_type: Used to decide between traditional or reformulated costs
			mod: used reformulated compression costs (see overleaf). They are used to get tighter bounds, but can lead to execution issues.
			otherwise: use traditional compression costs. Can lead to worse bounds, but will have less execution errors.

		ged_method_options: For the moment its only the number of threads, as a string

		ged_method: GED method to use. Be careful because not all methods work with the compression edit costs, specially the the reformulated ones
			branch_uniform
			branch_fast
			ipfp

		graph_sample_size: The percentage of graphs to use for the calculation of the GED. Must be a positive number. For 50% use 50, not 0.5. If greater than 100, then all the GED's are calculated

		write_ged_matrix: Used to decide if the GED matrix should be written in csv format. It is valid for the initial and refined one.
			true: Write the matrix in the output_root using the file_preffix
			otherwise: Don't write the matrix

		write_arb: Used to decide if the vector representing the arborescence should be written. The format is the following: For a node x, then its parent is the node arb.at(x). The root of the arborescence is the node with value arb.size()
			true: Write the vector in the output_root using the file_preffix
			otherwise: Don't write the vector

		relaxed_compression: Wether to use the relaxed coding (allows insertions and suppresions of nodes at the same time) or the strict one. If strict mode is selected, then the traditional edit costs can lead to execution errors due to the prescence of insertions and suppresions simultaneously
			true: Use the relaxed compression format
			otherwise: Use the strict compression format (may lead to errors if traditional costs are used)
		
		write_results: Wheter to fill string vectors with the results of the compression process
			true: Do it
			otherwise: Don't

	*/
	template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
	void compress_collection(ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env, 
		std::string output_root,
		std::string file_preffix,
		std::string folder_for_encoded,
		std::map<std::string, std::string> &args,
		int stdout = 0,
		std::vector<std::string> &headers = nullptr,
		std::vector<std::string> &values = nullptr
		);

	template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
	void compress_collection_from_empty(ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env, 
		std::string output_root,
		std::string file_preffix,
		std::string folder_for_encoded,
		std::map<std::string, std::string> &args,
		int stdout = 0
		);


	template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
	void decode_collection( ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env,
		std::string output_root, std::string file_preffix, std::map<std::string, std::string> &args,
		int stdout=0);

	template<class T>
	void write_to_file(std::ofstream &file, std::vector<T> &values);

	template<class T>
	void write_to_file(std::string path, std::vector<T> &values);



private:

	double omega_prop_;

	template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
	void decode_single_graph(std::string path, std::string graph_name, std::string graph_class,
		ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env,
		std::size_t parent_num, std::size_t child,
		std::size_t root,
		std::map<std::size_t, std::size_t> &pos_to_id,
		std::size_t b_ni, std::size_t b_ei,
		std::map<std::string, std::vector<std::string>> &ordered_attributes,
		std::map<std::string, std::map<std::string, std::size_t>> &attr_sizes,
		int stdout = 0);


	template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
	void decode_single_graph_relaxed(std::string path, std::string graph_name, std::string graph_class,
		ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env,
		std::size_t parent_num, std::size_t child, std::size_t root,
		std::map<std::size_t, std::size_t> &pos_to_id,
		std::size_t b_ni, std::size_t b_ei,
		std::map<std::string, std::vector<std::string>> &ordered_attributes,
		std::map<std::string, std::map<std::string, std::size_t>> &attr_sizes,
		int stdout=0);

	template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
	void print_compression_sets(
		vector<ged::GEDGraph::NodeID> &v_d,
		vector<ged::GEDGraph::NodeID> &v_i,
		vector<ged::GEDGraph::NodeID> &v_s,
		vector<ged::GEDGraph::NodeID> &v_is,
		vector<std::pair<ged::GEDGraph::NodeID,ged::GEDGraph::NodeID>> &e_nd,
		vector<std::pair<ged::GEDGraph::NodeID,ged::GEDGraph::NodeID>> &e_ed,
		vector<std::pair<ged::GEDGraph::NodeID,ged::GEDGraph::NodeID>> &e_ni,
		vector<std::pair<ged::GEDGraph::NodeID,ged::GEDGraph::NodeID>> &e_ei,	
		vector<std::pair<ged::GEDGraph::NodeID,ged::GEDGraph::NodeID>> &e_s,
		vector<std::pair<ged::GEDGraph::NodeID,ged::GEDGraph::NodeID>> &e_is,
		ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &g1,
		ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &g2
		);


	void spanning_arborescence_of_minimum_weight(
		std::vector<std::size_t> &tree,
		std::size_t &cost,
		std::vector<std::vector<std::size_t>> &w,
		std::size_t &root,
		std::size_t max_arc_cost
		);

	void get_arborescence_info(
		std::vector<std::string> &headers,
		std::vector<std::string> &values,
		std::map<std::size_t, std::vector<std::size_t>> &depth_degrees,
		std::vector<std::size_t> &arborescence,
		std::size_t root);


	template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
	void translate_env(std::map<std::string, std::map<std::string, std::map<std::string,std::string>>> &encoded_attributes,
		ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env_orig,
		ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env_coded,
		bool fast_node_translate = false,
		bool fast_edge_translate = false,
		int stdout = 0);



	template<class T>
	std::map<std::string, std::vector<std::string>> get_ordered_attributes(std::map<std::string, std::map<std::string, T>> &container);


	ged::NodeMap get_aux_node_map(std::size_t v1, std::size_t v2, vector<ged::GEDGraph::NodeID> &v_d, vector<ged::GEDGraph::NodeID> &v_i, int stdout=0);



	ged::NodeMap get_id_node_map(std::size_t  num_nodes, ged::NodeMap node_map, ged::NodeMap node_map_aux, ged::NodeMap node_map_id, int stdout=0);


	void permute_nodes(ged::NodeMap permutation, vector<ged::GEDGraph::NodeID> &v);

	std::size_t get_size_t_from_bytes(std::ifstream &file, std::size_t num);

	std::vector<std::size_t> random_sample(std::vector<std::size_t> &population, std::size_t k, std::size_t curr_graph);


	void modified_costs(std::vector<double> & comp_costs, std::size_t v1, std::size_t v2, double &c_nd, double &c_ni, double &c_ns, double &c_ed, double &c_ei, double &c_es, double &c_es_id, std::size_t  &omega, std::size_t &b_ni,
		std::size_t & b_na,	std::size_t & b_ei, std::size_t & b_ea);


	void write_matrix(std::string path, std::vector<std::vector<std::size_t>> upper_bounds);

	std::map<std::string, std::map<std::string, std::map<std::string,std::string>>>
	get_attribute_encoding(std::map<std::string, std::map<std::string, std::set<std::string>>> &alphabets);


	template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
	void create_info_file(std::string path, std::string file_preffix, std::string folder_for_encoded,
		ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env, 
		std::map<std::string, std::map<std::string, std::set<std::string>>> &alphabets,
		std::map<std::string, std::map<std::string, std::size_t>> &attr_sizes,
		std::size_t b_ni, std::size_t b_ei,
		std::vector<std::size_t> &arborescence,
		std::size_t root,
		bool fast_node_translate, bool fast_edge_translate,
		int stdout = 0,	char separator = '\n'
		);

	template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
	void encode_single_graph(std::string path,
		ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env_coded, 
		std::size_t parent_num, std::size_t child,
		std::map<ged::GEDGraph::GraphID, ged::NodeMap> &graph_permutations,
		std::map<std::string, std::vector<std::string>> &ordered_attributes,
		std::map<std::string, std::map<std::string, std::size_t>> &attr_sizes,
		std::size_t b_ni, std::size_t b_ei,
		int stdout = 0
		);

	template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
	void encode_single_graph_relaxed(std::string path,
		ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env_coded, 
		std::size_t parent_num, std::size_t child,
		std::map<ged::GEDGraph::GraphID, ged::NodeMap> &graph_permutations,
		std::map<std::string, std::vector<std::string>> &ordered_attributes,
		std::map<std::string, std::map<std::string, std::size_t>> &attr_sizes,
		std::size_t b_ni, std::size_t b_ei,
		int stdout = 0
		);

	template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
	void encode_arborescence(std::string output_root, std::string file_preffix,	std::string folder_for_encoded,
		ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env_coded,
		std::map<std::string, std::map<std::string, std::size_t>> &attr_sizes,
		std::size_t b_ni, std::size_t b_ei,
		std::vector<std::size_t> &arborescence,
		std::size_t root,
		bool relaxed = false,
		int stdout=0
		);

	template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
	std::size_t compute_induced_compression_cost(
		ged::NodeMap node_map, 
		ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &g1,
		ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &g2,
		std::size_t &b_ni, std::size_t &b_na,
		std::size_t &b_ei, std::size_t &b_ea);


	bool compare_label(std::map<std::string, std::string> label_1, std::map<std::string, std::string> label_2);


	template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
	void get_all_compression_sets(ged::NodeMap node_map,
		vector<ged::GEDGraph::NodeID> &v_d,
		vector<ged::GEDGraph::NodeID> &v_i,
		vector<UserNodeLabel> &varphi_i,
		vector<ged::GEDGraph::NodeID> &v_s,
		vector<ged::GEDGraph::NodeID> &v_is,
		vector<UserNodeLabel> &varphi_s,
		vector<std::pair<ged::GEDGraph::NodeID,ged::GEDGraph::NodeID>> &e_nd,
		vector<std::pair<ged::GEDGraph::NodeID,ged::GEDGraph::NodeID>> &e_ed,
		vector<std::pair<ged::GEDGraph::NodeID,ged::GEDGraph::NodeID>> &e_ni,
		vector<std::pair<ged::GEDGraph::NodeID,ged::GEDGraph::NodeID>> &e_ei,
		vector<UserEdgeLabel> &phi_ni,
		vector<UserEdgeLabel> &phi_ei,
		vector<std::pair<ged::GEDGraph::NodeID,ged::GEDGraph::NodeID>> &e_s,
		vector<std::pair<ged::GEDGraph::NodeID,ged::GEDGraph::NodeID>> &e_is,
		vector<UserEdgeLabel> &phi_s,
		ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &g1,
		ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &g2);

	
	unsigned char* to_binary(std::size_t value, std::size_t num_chars);


	std::size_t interpret(unsigned char* oData, std::size_t start, std::size_t num);

	unsigned char* read_chars(std::ifstream &file, std::size_t num);

	void write_chars(std::ofstream &file, std::size_t value, std::size_t num_chars);

};
}

#include "arborescence_based_compression.cpp"


#endif /* COMPRESS_SRC_ABC_HPP_ */

