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
#include <dirent.h>
// Must include this file for the MSA calculation.
// Currently defining the types for the GED matrix in MSArbor.h to use the same types
#include "../ext/frangio68/Minimal-Spanning-Arborescence-master/MSArbor.C"

#define COMPRESS_EDIT_COST // For redefining behaviors in ged_env and ged_data
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

	/*!
	 * @brief sets the proportion of std::numeric_limits<std::size_t> to use as omega for the modified edit costs.
	 * @param[in] d Floating pint number between 0 and 1. The omega used in the end will be std::size_t cast of std::numeric_limits<std::size_t>::max()*d.
	 */
	void set_omega(double d);

	/*!
	 * @brief Prints a description of the graph and its attributes.
	 * @tparam UserNodeID Class of user-specific node IDs.
 	 * @tparam UserNodeLabel Class of user-specific node labels. If nodes are unlabeled, use ged::NoLabel or define your own dummy label class.
 	 * @tparam UserEdgeLabel Class of user-specific edge labels. If edges are unlabeled, use ged::NoLabel or define your own dummy label class.
	 * @param[in] g ExchangeGraph representation of the graph.
	 */
	template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
	void describe_graph(ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> g);

	/*!
	 * @brief Creates structures with the alphabets and distributions of the node and edge attributes. Calculates also the number of bytes needed for encoding graphs nodes, edges and attributes.
	 * @tparam UserNodeID Class of user-specific node IDs.
 	 * @tparam UserNodeLabel Class of user-specific node labels. If nodes are unlabeled, use ged::NoLabel or define your own dummy label class.
 	 * @tparam UserEdgeLabel Class of user-specific edge labels. If edges are unlabeled, use ged::NoLabel or define your own dummy label class.
	 * @param[in] env ged_env containing the graph collection.
	 * @param[in] distribution structure to fill with the different attribute distributions. Used as output.
	 * @param[in] alphabets structure to fill with the different attribute alphabets. Used as output.
	 * @param[in] attr_sizes structure to fill with the different attribute sizes in bytes for encoding and decoding. Used as output.
	 * @param[in] b_ni number of bytes needed to encode a node index. Used as output.
	 * @param[in] b_na number of bytes needed to encode a node attribute dictionary. Used as output.
	 * @param[in] b_ei number of bytes needed to encode an edge index. Used as output.
	 * @param[in] b_ea number of bytes needed to encode an edge attribute dictionary. Used as output.
	 * @param[in] fast_node_translate boolean indicating if the node attributes change over nodes, forcing to introduce a dummy value for default.
	 * @param[in] fast_edge_translate boolean indicating if the edge attributes change over nodes, forcing to introduce a dummy value for default.
	 */
	template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
	void get_graphs_structure(ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env,
		std::map<std::string, std::map<std::string, std::vector<std::string>>> &distribution,
		std::map<std::string, std::map<std::string, std::set<std::string>>> &alphabets,
		std::map<std::string, std::map<std::string, std::size_t>> &attr_sizes,
		std::size_t &b_ni, std::size_t &b_na, std::size_t &b_ei, std::size_t &b_ea, 
		bool &fast_node_translate, bool &fast_edge_translate);

	/*!
	 * @brief Calculates the number of bytes needed to encode a graph with V vertices and E edges if
	 * edges are saved as node pairs.
	 * @param[in] V number of nodes.
	 * @param[in] E number of edges.
	 * @param[in] b_ni number of bytes needed to encode a node index.
	 * @param[in] b_na number of bytes needed to encode a node attribute dictionary.
	 * @param[in] b_ei number of bytes needed to encode an edge index.
	 * @param[in] b_ea number of bytes needed to encode an edge attribute dictionary.
	 * @return The number of bytes needed to encode the graph in the edge pair format.
	 */
	std::size_t compute_cost_edge_pairs(std::size_t V, std::size_t E, std::size_t b_ni, std::size_t b_na, std::size_t b_ei, std::size_t b_ea);


	/*!
	 * @brief Calculates the number of bytes needed to encode a graph with V vertices and E edges if
	 * edges are saved as the upper half of a triangular adjacency matrix.
	 * @param[in] V number of nodes.
	 * @param[in] E number of edges.
	 * @param[in] b_ni number of bytes needed to encode a node index.
	 * @param[in] b_na number of bytes needed to encode a node attribute dictionary.
	 * @param[in] b_ei number of bytes needed to encode an edge index.
	 * @param[in] b_ea number of bytes needed to encode an edge attribute dictionary.
	 * @return The number of bytes needed to encode the graph in the triangular matrix format.
	 */
	std::size_t compute_cost_triangular_matrix(std::size_t V, std::size_t  E, std::size_t b_ni, std::size_t b_na, std::size_t b_ei, std::size_t b_ea);


	/*!
	 * @brief Calculates the number of bytes needed to encode a graph with V vertices and E edges if
	 * the ABC optimal format is used. See paper for details
	 * @param[in] V number of nodes.
	 * @param[in] E number of edges.
	 * @param[in] b_ni number of bytes needed to encode a node index.
	 * @param[in] b_na number of bytes needed to encode a node attribute dictionary.
	 * @param[in] b_ei number of bytes needed to encode an edge index.
	 * @param[in] b_ea number of bytes needed to encode an edge attribute dictionary.
	 * @return The number of bytes needed to encode the graph in the ABC format.
	 */
	std::size_t compute_cost_abc(std::size_t V, std::size_t  E, std::size_t b_ni, std::size_t b_na, std::size_t b_ei, std::size_t b_ea);


	/*!
	 * @brief Calculates the total cost of saving a collection of files using the edge pair format
	 * @tparam UserNodeID Class of user-specific node IDs.
 	 * @tparam UserNodeLabel Class of user-specific node labels. If nodes are unlabeled, use ged::NoLabel or define your own dummy label class.
 	 * @tparam UserEdgeLabel Class of user-specific edge labels. If edges are unlabeled, use ged::NoLabel or define your own dummy label class.
	 * @param[in] env ged_env containing the graph collection.
	 * @param[in] b_ni number of bytes needed to encode a node index.
	 * @param[in] b_na number of bytes needed to encode a node attribute dictionary.
	 * @param[in] b_ei number of bytes needed to encode an edge index.
	 * @param[in] b_ea number of bytes needed to encode an edge attribute dictionary.
	 * @param[in] constant Constant to add to the cost of each graph.
	 * @param[in] ignore_false whether to include the cost of the last graph. This is used because
	 * in the method, the empty graph is inserted as the last graph. 
	 * @return Total encoding cost for the collection using the edge pair format
	 */
	template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
	std::size_t base_compr_cost_edge_pairs(ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env, std::size_t &b_ni, std::size_t &b_na, std::size_t &b_ei, std::size_t &b_ea, std::size_t constant=0, bool ignore_last=false);

	/*!
	 * @brief Calculates the total cost of saving a collection of files using the triangular format
	 * @tparam UserNodeID Class of user-specific node IDs.
 	 * @tparam UserNodeLabel Class of user-specific node labels. If nodes are unlabeled, use ged::NoLabel or define your own dummy label class.
 	 * @tparam UserEdgeLabel Class of user-specific edge labels. If edges are unlabeled, use ged::NoLabel or define your own dummy label class.
	 * @param[in] env ged_env containing the graph collection.
	 * @param[in] b_ni number of bytes needed to encode a node index.
	 * @param[in] b_na number of bytes needed to encode a node attribute dictionary.
	 * @param[in] b_ei number of bytes needed to encode an edge index.
	 * @param[in] b_ea number of bytes needed to encode an edge attribute dictionary.
	 * @param[in] constant Constant to add to the cost of each graph.
	 * @param[in] ignore_false whether to include the cost of the last graph. This is used because
	 * in the method, the empty graph is inserted as the last graph. 
	 * @return Total encoding cost for the collection using the triangular matrix format
	 */
	template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
	std::size_t base_compr_cost_triangular_matrix(ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env, std::size_t &b_ni, std::size_t &b_na, std::size_t &b_ei, std::size_t &b_ea, std::size_t constant=0, bool ignore_last=false);

	/*!
	 * @brief Calculates the total cost of saving a collection of files using the ABC format. See paper for details.
	 * @tparam UserNodeID Class of user-specific node IDs.
 	 * @tparam UserNodeLabel Class of user-specific node labels. If nodes are unlabeled, use ged::NoLabel or define your own dummy label class.
 	 * @tparam UserEdgeLabel Class of user-specific edge labels. If edges are unlabeled, use ged::NoLabel or define your own dummy label class.
	 * @param[in] env ged_env containing the graph collection.
	 * @param[in] b_ni number of bytes needed to encode a node index.
	 * @param[in] b_na number of bytes needed to encode a node attribute dictionary.
	 * @param[in] b_ei number of bytes needed to encode an edge index.
	 * @param[in] b_ea number of bytes needed to encode an edge attribute dictionary.
	 * @param[in] constant Constant to add to the cost of each graph.
	 * @param[in] ignore_false whether to include the cost of the last graph. This is used because
	 * in the method, the empty graph is inserted as the last graph. 
	 * @return Total encoding cost for the collection using the ABC format
	 */
	template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
	std::size_t base_compr_cost_abc(ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env, std::size_t &b_ni, std::size_t &b_na, std::size_t &b_ei, std::size_t &b_ea, std::size_t constant=0, bool ignore_last=false);

	/*! NEEDS MODIFICATION
	 * @brief Compresses a collection of graphs using the ABC method.
	 * @tparam UserNodeID Class of user-specific node IDs.
 	 * @tparam UserNodeLabel Class of user-specific node labels. If nodes are unlabeled, use ged::NoLabel or define your own dummy label class.
 	 * @tparam UserEdgeLabel Class of user-specific edge labels. If edges are unlabeled, use ged::NoLabel or define your own dummy label class.
	 * @param[in] env ged_env containing the graph collection.
	 * @param[in] output_root the root for all the output, including compressed files and metadata files.
	 * @param[in] file_preffix preffix to use when naming output files.
	 * @param[in] folder_for_encoded folder inside output_root where to store the encoded files (compressed graphs).
	 * @param[in] args map containing different name-value pairs for changing the behaviour of the method 
	 *	Args can have the following fields:
	 *	edit_cost_type: Used to decide between traditional or reformulated costs
	 *		mod: used reformulated compression costs (see overleaf). They are used to get tighter bounds, but can lead to execution issues.
	 *		otherwise: use traditional compression costs. Can lead to worse bounds, but will have less execution errors.
	 *	ged_method_options: For the moment its only the number of threads, as a string
	 * 	ged_method: GED method to use. Be careful because not all methods work with the compression edit costs, specially the the reformulated ones
	 *		branch_uniform
	 *		branch_fast
	 *		ipfp
	 *	graph_sample_size: The percentage of graphs to use for the calculation of the GED. Must be a positive number. For 50% use 50, not 0.5. If greater than 100, then all the GED's are calculate
	 *  write_ged_matrix: Used to decide if the GED matrix should be written in csv format. It is valid for the initial and refined one.
	 *		true: Write the matrix in the output_root using the file_preffix
	 *		otherwise: Don't write the matrix
	 *	write_arb: Used to decide if the vector representing the arborescence should be written. The format is the following: For a node x, then its parent is the node arb.at(x). The root of the arborescence is the node with value arb.size()
	 *		true: Write the vector in the output_root using the file_preffix
	 *		otherwise: Don't write the vector
	 *	relaxed_compression: Whether to use the relaxed coding (allows insertions and suppresions of nodes at the same time) or the strict one. If strict mode is selected, then the traditional edit costs can lead to execution errors due to the prescence of insertions and suppresions simultaneously
	 *		true: Use the relaxed compression format
	 *		otherwise: Use the strict compression format (may lead to errors if traditional costs are used)	
	 *	write_results: Wheter to fill string vectors with the results of the compression process
	 *		true: Do it
	 *		otherwise: Don't
	 * @param[in] stdout integer used to control the output in the terminal.
	 * @param[in] headers string vector to fill the headers of the output results file. Used as output.
	 * @param[in] values string vector to fill results. Used as output.
	*/
	template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
	void compress_collection(bool binary, bool relaxed, std::string graph_dir, std::string xml_file,
		Options::GXLNodeEdgeType node_type, Options::GXLNodeEdgeType edge_type,
		const std::unordered_set<std::string> & irrelevant_node_attributes, 
		const std::unordered_set<std::string> & irrelevant_edge_attributes,
		std::map<std::string, std::map<std::string, char>> & attr_types,
		std::string output_encoded,
		std::string output_other,
		std::string file_preffix,
		std::map<std::string, std::string> &args,
		int stdout = 0,
		std::vector<std::string> &headers = nullptr,
		std::vector<std::string> &values = nullptr,
		char separator_1 = "\n",
		char separator_2 = "\0"
		);

	/*!
	 * @brief Loads a collection of graphs from their comressed ABC format to a ged_env object.
	 * @tparam UserNodeID Class of user-specific node IDs.
 	 * @tparam UserNodeLabel Class of user-specific node labels. If nodes are unlabeled, use ged::NoLabel or define your own dummy label class.
 	 * @tparam UserEdgeLabel Class of user-specific edge labels. If edges are unlabeled, use ged::NoLabel or define your own dummy label class.
	 * @param[in] env ged_env to store the graph collection.
	 * @param[in] output_root the root for all the output, including compressed files and metadata files.
	 * @param[in] file_preffix preffix to use when naming output files.
	 * @param[in] folder_for_encoded folder inside output_root where to store the encoded files (compressed graphs).
	 * @param[in] args map containing different name-value pairs for changing the behaviour of the method 
	 *	Args can have the following fields:
	 *	relaxed_compression: Whether to use the relaxed coding (allows insertions and suppresions of nodes at the same time) or the strict one. If strict mode is selected, then the traditional edit costs can lead to execution errors due to the prescence of insertions and suppresions simultaneously
	 *		true: Use the relaxed compression format
	 *		otherwise: Use the strict compression format (may lead to errors if traditional costs are used)	
	 * @param[in] stdout integer used to control the output in the terminal.
	*/
	template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
	void decode_collection(bool binary, bool relaxed, ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env,
		std::string path, std::string file_preffix, std::map<std::string, std::string> &args,
		int stdout=0, 	char separator_1 = "\n",	char separator_2 = "\0");

	/*!
	 * @brief Writes a vector to a file using comma as separator and \n as a line breaker
	 * @tparam T type inside the container.
	 * @param[in] file ofstream object.
	 * @param[in] values vector to write to the file.
	*/
	template<class T>
	void write_to_file(std::ofstream &file, std::vector<T> &values);

	/*!
	 * @brief Writes a vector to a file using comma as separator and \n as a line breaker
	 * @tparam T type inside the container.
	 * @param[in] path path to the file.
	 * @param[in] values vector to write to the file.
	*/
	template<class T>
	void write_to_file(std::string path, std::vector<T> &values);

	// Take the xml collection file as the meta data for the compression. Fills the list graphs with the <name, class> of each graph
	void read_xml_graph_collection(const std::string & file, std::list<std::pair<std::string, std::string>> & graphs);

	// Returns the size in bytes of a file
	std::size_t get_file_size(std::string path);

	// Used to determine an order in the attribute reading and writing. 
	// This ensures that at decompression, the attributes are read in the same order they were written when comrpessing
	template<class T>
	std::map<std::string, std::vector<std::string>> get_ordered_attributes(std::map<std::string, std::map<std::string, T>> &container);

		// Creates a new env with the attributes translated according to encoded_attributes
	template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
	void translate_env(std::map<std::string, std::map<std::string, std::map<std::string,std::string>>> &encoded_attributes,
		ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env_orig,
		ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env_coded,
		bool fast_node_translate = false,
		bool fast_edge_translate = false,
		int stdout = 0);

	// Given the alphabets, generates a map encoding the attributes with their index in the alphabet.
	std::map<std::string, std::map<std::string, std::map<std::string,std::string>>>
	get_attribute_encoding(std::map<std::string, std::map<std::string, std::set<std::string>>> &alphabets);


	// Returns a sorted vector (lexicographic order) of the file names inside path whose names end with suffix
	std::vector<std::string> get_sorted_file_names(std::string path, std::string suffix);

private:
	// Number between 0 and 1. Determines the omega used in the modified edit costs as
	// std::size_t omega = std::numeric_limits<std::size_t>::max() * omega_prop_;
	double omega_prop_;


	template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
	ged::NodeMap
	get_trivial_node_map(std::vector<std::pair<std::size_t, std::string>> &att_1, std::vector<std::pair<std::size_t, std::string>> &att_2);

	// Loads the graph contained in path to the environment env.
	template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
	void decode_single_graph(bool binary, bool relaxed, std::string path, std::string graph_name, std::string graph_class,
		ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env,
		std::size_t parent_num, std::size_t child,
		std::size_t root,
		std::map<std::size_t, std::size_t> &pos_to_id,
		std::size_t b_ni, std::size_t b_ei,
		std::map<std::string, std::vector<std::string>> &ordered_attributes,
		std::map<std::string, std::map<std::string, std::size_t>> &attr_sizes,
		int stdout = 0, char separator_1 = '\n');

	// Loads the graph contained in path to the environment env using the relaxed format.
	template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
	void decode_single_graph_relaxed(bool binary, std::string path, std::string graph_name, std::string graph_class,
		ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env,
		std::size_t parent_num, std::size_t child, std::size_t root,
		std::map<std::size_t, std::size_t> &pos_to_id,
		std::size_t b_ni, std::size_t b_ei,
		std::map<std::string, std::vector<std::string>> &ordered_attributes,
		std::map<std::string, std::map<std::string, std::size_t>> &attr_sizes,
		int stdout=0);

	// Print in an ordered way the compression sets
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


	// Uses the implementation of frangio68 to calculate the MSA from the cost matrix csts
	void spanning_arborescence_of_minimum_weight(
		std::vector<std::size_t> &tree,
		std::size_t &cost,
		MSA_di_unipi_it::MSArbor::CRow csts,
		std::size_t &root,
		std::size_t max_arc_cost
		);

	// Used to get some statistics from the arborescence
	void get_arborescence_info(
		std::vector<std::string> &headers,
		std::vector<std::string> &values,
		std::map<std::size_t, std::vector<std::size_t>> &depth_degrees,
		std::vector<std::size_t> &arborescence,
		std::size_t root);




	
	// Creates the auxiliar node map used for the decompression path.
	ged::NodeMap get_aux_node_map(std::size_t v1, std::size_t v2, vector<ged::GEDGraph::NodeID> &v_d, vector<ged::GEDGraph::NodeID> &v_i, int stdout=0);


	// Uses the original node map and the auxiliar node map to create the node map between 
	// a graph and its decompressed version.
	ged::NodeMap get_id_node_map(std::size_t  num_nodes, ged::NodeMap node_map, ged::NodeMap node_map_aux, ged::NodeMap node_map_id, int stdout=0);


	// Permutes the nodes of a graph according to a node map. 
	// User must be sure that the node map permutation has an image for every element in v.
	// If not, it will result in an error
	void permute_nodes(ged::NodeMap permutation, vector<ged::GEDGraph::NodeID> &v);

	// Reads num bytes from file and interpret this bytes as an unsigned integer
	std::size_t get_size_t_from_bytes(std::ifstream &file, std::size_t num);

	// Creates a random uniform sample of size k from population without including curr_graph in the result. 
	std::vector<std::size_t> random_sample(std::vector<std::size_t> &population, std::size_t k, std::size_t curr_graph);

	// Calculates the edit cost constants depending on the case (v1 < v2 or v2 <= v1). Only applies for transformed costs
	void modified_costs(std::vector<double> & comp_costs, std::size_t v1, std::size_t v2, double &c_nd, double &c_ni, double &c_ns, double &c_ed, double &c_ei, double &c_es, double &c_es_id, std::size_t  &omega, std::size_t &b_ni,
		std::size_t & b_na,	std::size_t & b_ei, std::size_t & b_ea);

	// Write the GED matrix into a file using comma as separator and \n as line break.
	// The matrix is stored in a MSA_di_unipi_it::MSArbor::CRow, which is a pointer to the type defined in
	// MSArbor.h as Index (currently unsigned int).
	void write_matrix(std::string path, MSA_di_unipi_it::MSArbor::CRow upper_bounds,std::size_t lines);

	

	// Writes the info_file containing all the metadata required for decompression.
	template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
	void create_info_file(bool binary, std::string path, std::string file_preffix,
		ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env, 
		std::map<std::string, std::map<std::string, std::set<std::string>>> &alphabets,
		std::map<std::string, std::map<std::string, std::size_t>> &attr_sizes,
		std::map<std::string, std::map<std::string, char>> &attr_types,
		std::size_t b_ni, std::size_t b_ei,
		std::vector<std::size_t> &arborescence,
		std::size_t root,
		bool fast_node_translate, bool fast_edge_translate,
		int stdout,	char separator_1, char separator_2
		);

	// Creates the graph compressed file.
	template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
	void encode_single_graph(bool binary,bool relaxed, std::string path,
		ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env_coded, 
		std::size_t parent_num, std::size_t child,
		std::map<ged::GEDGraph::GraphID, ged::NodeMap> &graph_permutations,
		std::map<std::string, std::vector<std::string>> &ordered_attributes,
		std::map<std::string, std::map<std::string, std::size_t>> &attr_sizes,
		std::size_t b_ni, std::size_t b_ei,
		int stdout = 0, char separator = '\n'
		);

	// Creates the graph compressed file using the relaxed format.
	template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
	void encode_single_graph_relaxed(bool binary, std::string path,
		ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env_coded, 
		std::size_t parent_num, std::size_t child,
		std::map<ged::GEDGraph::GraphID, ged::NodeMap> &graph_permutations,
		std::map<std::string, std::vector<std::string>> &ordered_attributes,
		std::map<std::string, std::map<std::string, std::size_t>> &attr_sizes,
		std::size_t b_ni, std::size_t b_ei,
		int stdout = 0, char separator = '\n'
		);

	// Encodes the collection of graphs using the structure given in arborescence
	template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
	void encode_arborescence(bool binary, bool relaxed, std::string path, std::string file_preffix,
		ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env_coded,
		std::map<std::string, std::map<std::string, std::size_t>> &attr_sizes,
		std::size_t b_ni, std::size_t b_ei,
		std::vector<std::size_t> &arborescence,
		std::size_t root,
		int stdout=0, char separator='\n'
		);


	// Returns the compression cost induced by node_map between g1 and g2.
	template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
	std::size_t compute_induced_compression_cost(
		ged::NodeMap node_map, 
		ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &g1,
		ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &g2,
		std::size_t &b_ni, std::size_t &b_na,
		std::size_t &b_ei, std::size_t &b_ea);


	// Comparator for labels
	bool compare_label(std::map<std::string, std::string> label_1, std::map<std::string, std::string> label_2);

	// Comparator for labels
	template<class A, class B>
	bool compare_pair_second(std::pair<std::size_t, std::string> p1,std::pair<std::size_t, std::string> p2);



	// Calculates the compression sets from the node map and the graphs.
	// All the input vectors are filled and used to return output.
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

	
	// Writes value in binary format using num_chars bytes (if possible), adding trailing 0 if necessary
	unsigned char* size_t_to_binary(std::size_t value, std::size_t num_chars);


	// Reads a series of characters and translates their value to a unsigned integer.
	std::size_t interpret_binary_size_t(unsigned char* oData, std::size_t start, std::size_t num);

	// Reads num chars one by one from file and returns an array with them
	unsigned char* read_chars(std::ifstream &file, std::size_t num);

	// Writes value in binary using num_chars bytes.
	void write_size_t_to_chars(std::ofstream &file, std::size_t value, std::size_t num_chars);

	// To order the <name,class> pairs using the name attribute. 
	static bool compare_graph_names(std::pair<std::string, std::string> p1, std::pair<std::string, std::string> p2);


	void write_word_binary(std::ofstream & out_file, std::string word);

	std::string	read_word_binary(std::ifstream & in_file, char sep);

};
}

#include "arborescence_based_compression.cpp"


#endif /* COMPRESS_SRC_ABC_HPP_ */

