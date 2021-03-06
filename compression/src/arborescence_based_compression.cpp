/***************************************************************************
*                                                                          *
*   Copyright (C) 2020 by Lucas Gnecco                                     *
*                                                                          *
*   This file is part of GEDLIB-compression.                               *
*                                                                          *
*   GEDLIB is free software: you can redistribute it and/or modify it      *
*   under the terms of the GNU Lesser General Public License as published  *
*   by the Free Software Foundation, either version 3 of the License, or   *
*   (at your option) any later version.                                    *
*                                                                          *
*   GEDLIB is distributed in the hope that it will be useful,              *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of         *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the           *
*   GNU Lesser General Public License for more details.                    *
*                                                                          *
*   You should have received a copy of the GNU Lesser General Public       *
*   License along with GEDLIB. If not, see <http://www.gnu.org/licenses/>. *
*                                                                          *
***************************************************************************/

/*!
 * @file  arborescence_based_compression.cpp
 * @brief ged::GED_ABC class definition with the compression methods.
 */

#ifndef COMPRESS_SRC_ABC_IPP_
#define COMPRESS_SRC_ABC_IPP_

#define bytes(num) ceil(log2(num+1)/CHAR_BIT)



class compression_exception : public std::exception{
	public:
	    compression_exception(const std::string& msg) : m_msg(msg)
	    {}

	   ~compression_exception()
	   {
	        std::cout << "compression_exception::~compression_exception" << std::endl;
	   }

	   virtual const char* what() const throw () 
	   {
	        std::cout << "compression_exception - what:" << std::endl;
	        return m_msg.c_str();
	   }

	   const std::string m_msg;
};

namespace ged{


GED_ABC::
~GED_ABC() {
	
}

GED_ABC::
GED_ABC():
omega_prop_{0.2}
{}

void
GED_ABC::
set_omega(double d){
	if(d>1){
		std::cout<<"Wrong value. Omega must be less than 1"<<std::endl;
	}
	else{
		omega_prop_ = d;	
	}
}

std::vector<std::string>
GED_ABC::
get_sorted_file_names(std::string path, std::string suffix){
	DIR* rep = NULL;
    rep = opendir(path.c_str());
    if (rep == NULL) 
        throw(compression_exception("Could not open path to get sorted file names"));

    struct dirent* file = NULL;
    std::vector<std::string> filenames;
    std::string name;     
    std::size_t pos;
    while ((file = readdir(rep)) != NULL){
    	name = file->d_name;
    	if(suffix.empty()){
    		filenames.emplace_back(name);
    	}
    	else{
    		pos = (name.length()>=suffix.length()) ? name.length()-suffix.length() : 0;
	    	if(name.substr(pos) == suffix){
	    		filenames.emplace_back(name);	
	    	}
    	}
    }

    if (closedir(rep) == -1) 
        throw(compression_exception("Could not close path to get sorted file names"));

    std::sort(filenames.begin(), filenames.end());
    return filenames;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void 
GED_ABC::
describe_graph(ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> g){
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
ged::NodeMap
GED_ABC::
get_trivial_node_map(std::vector<std::pair<std::size_t, std::string>> &att_1, std::vector<std::pair<std::size_t, std::string>> &att_2){

	ged::NodeMap node_map = ged::NodeMap(att_1.size(), att_2.size());
	std::vector<bool> missing(att_2.size(), true);

	// Sorting can help find things easier, but then the real positions
	// of nodes are lost
	std::vector<std::string> keys_2;
	for(std::size_t i = 0; i<att_2.size(); i++){
		keys_2.emplace_back(att_2.at(i).second);
	} 

	typename std::vector<std::string>::iterator it;
	typename std::vector<std::string>::iterator beg;
	bool cont = false;
	//std::cout<<"Start trivial Noe map"<<std::endl;
	for(std::size_t i=0; i<att_1.size(); i++){
		cont = false;
		beg = keys_2.begin();
		while(!cont){
			it = std::find(beg, keys_2.end(), att_1.at(i).second);
			if(it == keys_2.end()){
				// Erase node
				node_map.add_assignment(att_1.at(i).first, ged::GEDGraph::dummy_node());
				cont=true;
			}
			else{				
				// Just in case, check this node has not been found
				if(! missing.at(it-keys_2.begin())){
					// This node is already in the nodemap
					beg = it+1;					
				}
				else{
					// Add assignment
					node_map.add_assignment(att_1.at(i).first, att_2.at(it-keys_2.begin()).first);	
					missing.at(it-keys_2.begin())	= false;
					cont=true;
				}				
			}
		}
	}
	// If nodes in g2 are missing, we have to add them
	for(std::size_t i = 0; i<missing.size(); i++){
		if(missing.at(i)){
			node_map.add_assignment(ged::GEDGraph::dummy_node(), att_2.at(i).first);
		}
	}

	return node_map;

}


unsigned char*
GED_ABC::
size_t_to_binary(std::size_t value, std::size_t num_chars){
	std::size_t aux = value;
	std::size_t mod = 1;
	std::size_t base = pow(2, CHAR_BIT);
	std::size_t bytes = bytes(value);
	if(bytes > num_chars){
		throw(compression_exception("to_binary: overflow"));
	}
	unsigned char * line = new unsigned char[num_chars];
	
	for(std::size_t i =0; i< num_chars; i++){
		if(i< bytes){
			mod = aux % base;
			line[i] = static_cast<unsigned char>(mod);		
			aux = aux / base;	
		}
		else{
			line[i] = static_cast<unsigned char>(0);
		}
		
	}

	return line;
}


std::size_t 
GED_ABC::
interpret_binary_size_t(unsigned char* oData, std::size_t start, std::size_t num){
	std::size_t base = 1;
	std::size_t sum=0;
	std::size_t index=0;
	//index = start + num - 1;
	index = start;

	for ( std::size_t i = 0; i< num ; i++ )
	{		
		sum = sum + base * static_cast<std::size_t> (oData[index]);
		base = base * pow(2,CHAR_BIT);
		index++;	
	}
	return sum;
}

std::size_t 
GED_ABC::
get_size_t_from_bytes(std::ifstream &file, std::size_t num){
	unsigned char* line = read_chars(file, num);
	return interpret_binary_size_t(line, 0, num);
}


unsigned char* 
GED_ABC::
read_chars(std::ifstream &file, std::size_t num){
	unsigned char* line = new unsigned char[num];
	unsigned char c;
	for(std::size_t i = 0; i<num; i++){
		file.read(reinterpret_cast<char*> (&c), 1);
		line[i] = c;
	}
	return line;
}

void 
GED_ABC::
write_size_t_to_chars(std::ofstream &file, std::size_t value, std::size_t num_chars){
	unsigned char* res;
	res = size_t_to_binary(value, num_chars); // use the minimal number of chars needed
	file.write(reinterpret_cast<char*> (&res[0]), num_chars);
}

void
GED_ABC::
write_word_binary(std::ofstream & out_file, std::string word){
	for(auto s_iter = word.begin(); s_iter != word.end(); s_iter++){
        out_file.put((*s_iter));
    }
}

std::string
GED_ABC::
read_word_binary(std::ifstream & in_file, char sep){
	std::string ans = "";
	char c;
	in_file.get(c);
	while(c != sep){
		ans.push_back(c);
		in_file.get(c);
	}
	return ans;
}


template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void 
GED_ABC::
get_graphs_structure(ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env,
	std::map<std::string, std::map<std::string, std::vector<std::string>>> &distribution,
	std::map<std::string, std::map<std::string, std::set<std::string>>> &alphabets,
	std::map<std::string, std::map<std::string, std::size_t>> &attr_sizes,
	std::size_t &b_ni, std::size_t &b_na, std::size_t &b_ei, std::size_t &b_ea, 
	bool &fast_node_translate, bool &fast_edge_translate){

	std::map<std::string, std::vector<std::string>> node_attr;
	std::map<std::string, std::vector<std::string>> edge_attr;
	std::map<std::string, std::set<std::string>> node_attr_set;
	std::map<std::string, std::set<std::string>> edge_attr_set;
	std::map<std::string, std::size_t> node_attr_size;
	std::map<std::string, std::size_t> edge_attr_size;

	std::size_t max_nodes = 0;
	std::size_t max_edges = 0;
	std::size_t first_label_found_pos_node = 0;
	std::size_t first_label_found_graph_node = 0;
	bool first_label_found_node = false;
	std::size_t first_label_found_pos_edge = 0;
	std::size_t first_label_found_graph_edge = 0;
	bool first_label_found_edge = false;

	fast_node_translate = true;
	fast_edge_translate = true;
	std::size_t cont;

	ged::ExchangeGraph<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> g;
	std::pair<ged::GEDGraph::GraphID, ged::GEDGraph::GraphID> graph_ids = env.graph_ids();
	
	for (ged::GEDGraph::GraphID g_id{graph_ids.first}; g_id<graph_ids.second; g_id++){
		g = env.get_graph(g_id, false, false, true);
		max_nodes = max(max_nodes, g.num_nodes);
		max_edges = max(max_edges, g.num_edges);
		for (std::size_t i{0}; i < g.num_nodes; i++) {
			for(auto l : g.node_labels.at(i)){
				if(node_attr.count(l.first)==0){
					if(!first_label_found_node){
						first_label_found_node = true;
						first_label_found_pos_node = i;
						first_label_found_graph_node = g_id;
					}
					if(g_id != first_label_found_graph_node || i > first_label_found_pos_node){
					// new label name after having already found a label
					// this means nodes have different attribute names
						fast_node_translate = false;
					}
					node_attr.emplace(std::make_pair(l.first, std::vector<std::string>()));
					node_attr_set.emplace(std::make_pair(l.first, std::set<std::string>()));
				}
				node_attr.at(l.first).emplace_back(l.second);
				node_attr_set.at(l.first).insert(l.second);
			}
		}

		std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel> edge;
		typename std::list<std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel>>::iterator iter;
		cont = 0;
		for(iter = g.edge_list.begin(); iter != g.edge_list.end(); iter ++ ){
			edge = (*iter);
			cont++;
			for(auto e : edge.second){
				if(edge_attr.count(e.first)==0){
					if(!first_label_found_edge){
						first_label_found_edge = true;
						first_label_found_pos_edge = cont;
						first_label_found_graph_edge = g_id;
					}
					if(g_id != first_label_found_graph_edge || cont > first_label_found_pos_edge){
					// new label name after having already found a label
					// this means edges have different attribute names
						fast_edge_translate = false;
					}
					edge_attr.emplace(std::make_pair(e.first, std::vector<std::string>()));
					edge_attr_set.emplace(std::make_pair(e.first, std::set<std::string>()));
				}
				edge_attr.at(e.first).emplace_back(e.second);
				edge_attr_set.at(e.first).insert(e.second);
			}
		}	
	}

	distribution.clear();
	distribution.emplace(std::make_pair("node_attr", node_attr));
	distribution.emplace(std::make_pair("edge_attr", edge_attr));

	alphabets.clear();
	alphabets.emplace(std::make_pair("node_attr", node_attr_set));
	alphabets.emplace(std::make_pair("edge_attr", edge_attr_set));

	b_ni = static_cast<std::size_t> (bytes(max_nodes));
	b_ei = static_cast<std::size_t> (bytes(max_edges));
	b_na = 0;
	b_ea = 0;
	std::size_t aux;
	for(auto const& a: alphabets.at("node_attr")){
		aux = a.second.size();
		if(!fast_node_translate){
			aux += 1;
		}
		node_attr_size.emplace(std::make_pair(a.first, bytes(aux)));
		b_na += static_cast<std::size_t> (bytes(aux));
	}
	for(auto const& a: alphabets.at("edge_attr")){
		aux = a.second.size();
		if(!fast_node_translate){
			aux += 1;
		}
		edge_attr_size.emplace(std::make_pair(a.first, bytes(aux)));
		b_ea += static_cast<std::size_t> (bytes(aux));
	}

	attr_sizes.clear();
	attr_sizes.emplace(std::make_pair("node_attr", node_attr_size));
	attr_sizes.emplace(std::make_pair("edge_attr", edge_attr_size));

}

bool
GED_ABC::
compare_label(std::map<std::string, std::string> label_1, std::map<std::string, std::string> label_2){
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


struct {
    bool operator()(std::pair<std::size_t, std::string> a, std::pair<std::size_t, std::string> b) const
    {   
        return a.second < b.second;
    }   
} pair_sort;

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void 
GED_ABC::
get_all_compression_sets(ged::NodeMap node_map,
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
	ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &g2){
	
	std::size_t index;

	//Nodes
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
	e_nd.clear();
	e_ed.clear();
	e_ni.clear();
	e_ei.clear();
	phi_ni.clear();
	phi_ei.clear();
	e_s.clear();
	e_is.clear();
	phi_s.clear();

	std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel> edge;
	typename std::list<std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel>>::iterator iter;
	index = 0;
	for(iter = g1.edge_list.begin(); iter != g1.edge_list.end(); iter ++ ){
		
		edge = (*iter);
		if(node_map.image(edge.first.first)==ged::GEDGraph::dummy_node() || node_map.image(edge.first.second)==ged::GEDGraph::dummy_node()){
			e_nd.emplace_back(edge.first); // in G1.edge_list
		}
		else{
			if(g2.adj_matrix[node_map.image(edge.first.first)][node_map.image(edge.first.second)]==0 && g2.adj_matrix[node_map.image(edge.first.second)][node_map.image(edge.first.first)]==0){
				e_ed.emplace_back(edge.first); // in G1.edge_list
			}
			else{
				if(g2.adj_matrix[node_map.image(edge.first.first)][node_map.image(edge.first.second)]==1){

					if(compare_label(edge.second,g2.edge_labels.at(std::make_pair(node_map.image(edge.first.first),node_map.image(edge.first.second))))){
						e_is.emplace_back(edge.first); // in G1.edge_list
					}
					else{
						e_s.emplace_back(edge.first); // in G1.edge_list
						phi_s.emplace_back(g2.edge_labels.at(std::make_pair(node_map.image(edge.first.first),node_map.image(edge.first.second))));
					}
				}
				else{
					if(compare_label(edge.second,g2.edge_labels.at(std::make_pair(node_map.image(edge.first.second),node_map.image(edge.first.first))))){
						e_is.emplace_back(edge.first); // in G1.edge_list
					}
					
					else{
						e_s.emplace_back(edge.first); // in G1.edge_list
						phi_s.emplace_back(g2.edge_labels.at(std::make_pair(node_map.image(edge.first.second),node_map.image(edge.first.first))));
					}
				}
			}		
		}
		index++;
	}
	index = 0;
	std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel> edge2;

	// Uses the adjacency matrix. A variant can be searching for edges in the edge_list. 
	// This would allow for stronger comparisson using labels. For now it works for graphs with no multiple edges
	for(iter = g2.edge_list.begin(); iter != g2.edge_list.end(); iter ++ ){
		edge2 = (*iter);
		if(node_map.pre_image(edge2.first.first)==ged::GEDGraph::dummy_node() || node_map.pre_image(edge2.first.second)==ged::GEDGraph::dummy_node()){
			e_ni.emplace_back(edge2.first); // in G2.edge_list
			phi_ni.emplace_back(edge2.second);			
		}
		else{
			if(g1.adj_matrix[node_map.pre_image(edge2.first.first)][node_map.pre_image(edge2.first.second)]==0 && g1.adj_matrix[node_map.pre_image(edge2.first.second)][node_map.pre_image(edge2.first.first)]==0){
				e_ei.emplace_back(edge2.first); // in G2.edge_list
				phi_ei.emplace_back(edge2.second);
			}
		}
		index++;
	}
	
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
std::size_t
GED_ABC::
compute_induced_compression_cost(
	ged::NodeMap node_map, 
	ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &g1,
	ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> &g2,
	std::size_t &b_ni, std::size_t &b_na,
	std::size_t &b_ei, std::size_t &b_ea){
	
	std::size_t result=0;
	std::size_t for_v_d=0;
	std::size_t for_v_is=0;
	std::size_t for_e_d=0;
	std::size_t for_e_is=0;

	//Nodes
	std::vector<ged::NodeMap::Assignment> relation;
	node_map.as_relation(relation);

	for(auto const assignment : relation){
		if(assignment.second == ged::GEDGraph::dummy_node()){
			for_v_d += b_ni; // deletion
		}
		if(assignment.first == ged::GEDGraph::dummy_node()){
			result += b_na;	// insertion
		}
		if(assignment.first != ged::GEDGraph::dummy_node() &&  assignment.second != ged::GEDGraph::dummy_node()){
			if(compare_label(g1.node_labels.at(assignment.first),g2.node_labels.at(assignment.second))){
				for_v_is += b_ni; // identical substitution
			}			
			else{
				result += b_ni + b_na; // substitution
			}			
		}
	}

	if(for_v_d < for_v_is){
		result += for_v_d;
	}
	else{
		result += for_v_is;
	}
	

	std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel> edge;
	typename std::list<std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel>>::iterator iter;	
	for(iter = g1.edge_list.begin(); iter != g1.edge_list.end(); iter ++ ){
		
		edge = (*iter);
		if(node_map.image(edge.first.first)==ged::GEDGraph::dummy_node() || node_map.image(edge.first.second)==ged::GEDGraph::dummy_node()){
			// Do nothing, edge deletion by node deletion.
		}
		else{
			if(g2.adj_matrix[node_map.image(edge.first.first)][node_map.image(edge.first.second)]==0 && g2.adj_matrix[node_map.image(edge.first.second)][node_map.image(edge.first.first)]==0){
				for_e_d += 2*b_ni; // deletion
			}
			else{
				if(g2.adj_matrix[node_map.image(edge.first.first)][node_map.image(edge.first.second)]==1){
					if(compare_label(edge.second,g2.edge_labels.at(std::make_pair(node_map.image(edge.first.first),node_map.image(edge.first.second))))){
						for_e_is += 2*b_ni; // identical substitution
					}
					else{
						result += 2*b_ni + b_ea; // substitution
					}
				}
				else{
					if(compare_label(edge.second,g2.edge_labels.at(std::make_pair(node_map.image(edge.first.second),node_map.image(edge.first.first))))){
						for_e_is += 2*b_ni; // identical substitution
					}
					else{
						result += 2*b_ni + b_ea; // substitution
					}
				}
			}		
		}
		
	}

	if(for_e_d < for_e_is){
		result += for_e_d;
	}
	else{
		result += for_e_is;
	}
	std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel> edge2;

	// Uses the adjacency matrix. A variant can be searching for edges in the edge_list. 
	// This would allow for stronger comparisson using labels. For now it works for graphs with no multiple edges
	for(iter = g2.edge_list.begin(); iter != g2.edge_list.end(); iter ++ ){
		edge2 = (*iter);
		if(node_map.pre_image(edge2.first.first)==ged::GEDGraph::dummy_node() || node_map.pre_image(edge2.first.second)==ged::GEDGraph::dummy_node()){
			result += 2*b_ni + b_ea;		
		}
		else{
			if(g1.adj_matrix[node_map.pre_image(edge2.first.first)][node_map.pre_image(edge2.first.second)]==0 && g1.adj_matrix[node_map.pre_image(edge2.first.second)][node_map.pre_image(edge2.first.first)]==0){
				result += 2*b_ni + b_ea;
			}
		}
		
	}

	return result;
}

std::size_t 
GED_ABC::
compute_cost_edge_pairs(std::size_t V, std::size_t  E, std::size_t b_ni, std::size_t b_na, std::size_t b_ei, std::size_t b_ea){return b_ni + V * b_na + E*(2*b_ni + b_ea);}

std::size_t 
GED_ABC::
compute_cost_triangular_matrix(std::size_t  V, std::size_t  E, std::size_t b_ni, std::size_t b_na, std::size_t b_ei, std::size_t b_ea){return b_ni + V*(V-1)/2 + V*b_na + E*b_ea;}

std::size_t 
GED_ABC::
compute_cost_abc(std::size_t V, std::size_t  E, std::size_t b_ni, std::size_t b_na, std::size_t b_ei, std::size_t b_ea){return 2*b_ni + V * b_na + 1 + 3*b_ei + E*(2*b_ni + b_ea);}


template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
std::size_t
GED_ABC::
base_compr_cost_edge_pairs(ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env, std::size_t &b_ni, std::size_t &b_na, std::size_t &b_ei, std::size_t &b_ea, std::size_t constant, bool ignore_last){
	std::size_t total_cost_compression = 0;
	std::pair<ged::GEDGraph::GraphID, ged::GEDGraph::GraphID> limits = env.graph_ids();
	std::size_t last = ignore_last ? limits.second-1 : limits.second;
	for(ged::GEDGraph::GraphID i{limits.first}; i<last; i++){
		total_cost_compression += compute_cost_edge_pairs(env.get_num_nodes(i), env.get_num_edges(i), b_ni, b_na, b_ei, b_ea);
		total_cost_compression += constant;
	}
	return total_cost_compression;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
std::size_t 
GED_ABC::
base_compr_cost_triangular_matrix(ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env, std::size_t &b_ni, std::size_t &b_na, std::size_t &b_ei, std::size_t &b_ea, std::size_t constant, bool ignore_last){
	std::size_t total_cost_compression = 0;
	std::pair<ged::GEDGraph::GraphID, ged::GEDGraph::GraphID> limits = env.graph_ids();
	std::size_t last = ignore_last ? limits.second-1 : limits.second;
	for(ged::GEDGraph::GraphID i{limits.first}; i<last; i++){
		total_cost_compression += compute_cost_triangular_matrix(env.get_num_nodes(i), env.get_num_edges(i), b_ni, b_na, b_ei, b_ea);
		total_cost_compression += constant;
	}
	return total_cost_compression;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
std::size_t 
GED_ABC::
base_compr_cost_abc(ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env, std::size_t &b_ni, std::size_t &b_na, std::size_t &b_ei, std::size_t &b_ea, std::size_t constant, bool ignore_last){
	std::size_t total_cost_compression = 0;
	std::pair<ged::GEDGraph::GraphID, ged::GEDGraph::GraphID> limits = env.graph_ids();
	std::size_t last = ignore_last ? limits.second-1 : limits.second;
	for(ged::GEDGraph::GraphID i{limits.first}; i<last; i++){
		total_cost_compression += compute_cost_abc(env.get_num_nodes(i), env.get_num_edges(i), b_ni, b_na, b_ei, b_ea);
		total_cost_compression += constant;
	}
	return total_cost_compression;
}


template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void 
GED_ABC::
print_compression_sets(
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
	){
	
	std::cout<<"-------------------NODES----------------------"<<std::endl;
	std::cout<<"v_d: deleted nodes"<<std::endl;
	for(const auto & n : v_d){
		std::cout<<n<<", ";
	}
	std::cout<<std::endl;
	std::cout<<"v_i: inserted nodes"<<std::endl;
	for(const auto & n : v_i){
		std::cout<<n<<", ";
	}
	std::cout<<std::endl;
	std::cout<<"v_is: identically substituted nodes"<<std::endl;
	for(const auto & n : v_is){
		std::cout<<n<<", ";
	}
	std::cout<<std::endl;
	std::cout<<"v_s: non-identically substituted nodes"<<std::endl;
	for(const auto & n : v_s){
		std::cout<<n<<", ";
	}
	std::cout<<std::endl;
	std::cout<<"-------------------EDGES----------------------"<<std::endl;


	std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel> edge;
	
	std::cout<<"e_nd: deleted edges by node deletion"<<std::endl;
	for(const auto & e : e_nd){
		std::cout<<e.first<<" -- "<<e.second<<", ";
	}
	std::cout<<std::endl;

	std::cout<<"e_ed: other deleted edges"<<std::endl;
	for(const auto & e : e_ed){
		std::cout<<e.first<<" -- "<<e.second<<", ";
	}
	std::cout<<std::endl;

	std::cout<<"e_is: identically substituted edges"<<std::endl;
	for(const auto & e : e_is){
		std::cout<<e.first<<" -- "<<e.second<<", ";
	}
	std::cout<<std::endl;

	std::cout<<"e_s: non-indentically substituted edges"<<std::endl;
	for(const auto & e : e_s){
		std::cout<<e.first<<" -- "<<e.second<<", ";
	}
	std::cout<<std::endl;

	std::cout<<"e_ni: inserted edges on new nodes"<<std::endl;
	for(const auto & e : e_ni){
		std::cout<<e.first<<" -- "<<e.second<<", ";
	}
	std::cout<<std::endl;

	std::cout<<"e_ei: other inserted edges"<<std::endl;
	for(const auto & e : e_ei){
		std::cout<<e.first<<" -- "<<e.second<<", ";
	}
	std::cout<<std::endl;

}

void 
GED_ABC::
spanning_arborescence_of_minimum_weight(
	std::vector<std::size_t> &tree,
	std::size_t &cost,
	MSA_di_unipi_it::MSArbor::CRow csts,
	std::size_t &root,
	std::size_t max_arc_cost
	){

	// FROM FISCHETTI 1993
	//https://github.com/frangio68/Minimal-Spanning-Arborescence
	// Depending on the dataset, it may be necessary to adjust the public types of 
	//	MSA_di_unipi_it::MSArbor::Index 
	//	MSA_di_unipi_it::MSArbor::CNumber
	// To do so, go to the original MSArbor.h file 
	
	MSA_di_unipi_it::MSArbor::Index n = root+1;
	MSA_di_unipi_it::MSArbor MSA( n );
	cost = MSA.Solve( csts );
	tree.clear();
	for(MSA_di_unipi_it::MSArbor::Index i = 0 ; i < n - 1 ; i++ ){
	 	tree.emplace_back(MSA.ReadPred()[ i ]);
	}

	return;
} 

void 
GED_ABC::
get_arborescence_info(
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

}

template<class T>
void 
GED_ABC::
write_to_file(std::ofstream &file, std::vector<T> &values){
	
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

template<class T>
void 
GED_ABC::
write_to_file(std::string path, std::vector<T> &values){
	std::ofstream file(path.c_str(), ios::app);
	for(std::size_t j=0; j<values.size(); j++){
		file << values.at(j);
		if(j < values.size()-1){
			file << ","; 	
		} 
		else{
			file << "\n"; 	
		}
	}
	file.close();
}  

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void 
GED_ABC::
translate_env(std::map<std::string, std::map<std::string, std::map<std::string,std::string>>> &encoded_attributes,
	ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env_orig,
	ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env_coded,
	bool fast_node_translate,
	bool fast_edge_translate,
	int stdout){

	std::pair<ged::GEDGraph::GraphID, ged::GEDGraph::GraphID> limits;
	limits = env_orig.graph_ids();	
	ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> g;
	std::pair<std::pair<std::size_t, std::size_t>, ged::GXLLabel> edge;
	ged::ProgressBar progress(limits.second);

	std::map<std::size_t, std::vector<std::string>> to_erase;
	std::size_t index;

	if (stdout > 0) std::cout << "\rTranslating graphs: " << progress << std::flush;
	for(ged::GEDGraph::GraphID i{limits.first}; i<limits.second; i++){
		g = env_orig.get_graph(i, false, false, true); // Edge list
		to_erase.clear();
		for(std::size_t n=0; n< g.node_labels.size(); n++){	
			to_erase.emplace(std::make_pair(n, std::vector<std::string>()));
			if(fast_node_translate){
				for(const auto dict : g.node_labels.at(n)){
					if(stdout>3) std::cout<<"Node :"<<n<<" --- "<<dict.first<<" , "<<dict.second<<std::endl;
					if(stdout>3) std::cout<<encoded_attributes.at("node_attr").at(dict.first).at(dict.second)<<std::endl;
					if(encoded_attributes.at("node_attr").at(dict.first).at(dict.second)=="dummy"){
						//g.node_labels.at(n).erase(dict.first);
						to_erase.at(n).emplace_back(dict.first);
					}
					else{
						g.node_labels.at(n).at(dict.first) = encoded_attributes.at("node_attr").at(dict.first).at(dict.second);
					}
				}
			}
			else{
				for(const auto dict : encoded_attributes.at("node_attr")){
					if(g.node_labels.at(n).count(dict.first)>0){
						if (stdout >3) std::cout<<"Node "<<n<<": "<<dict.first<<" -> "<<dict.second.at(g.node_labels.at(n).at(dict.first))<<std::endl;
						g.node_labels.at(n).at(dict.first) = dict.second.at(g.node_labels.at(n).at(dict.first));
					}
					else{
						// Add default value
						if (stdout >3) std::cout<<"Node "<<n<<": "<<dict.first<<" -> default"<<std::endl;
						g.node_labels.at(n).emplace(std::make_pair(dict.first, 
							dict.second.at("dummy")));
					}
				}	
			}
		}
		// Erase the dummy labels
		for(const auto & l : to_erase){
			for(const auto & a : l.second){
				if (stdout >3) std::cout<<"Erasing label: Node "<<l.first<<", "<<a<<std::endl;
				g.node_labels.at(l.first).erase(a);	
			}
			
		}	


		typename std::list<std::pair<std::pair<std::size_t, std::size_t>, ged::GXLLabel>>::iterator iter;
		index=0;
		to_erase.clear();
		for(iter = g.edge_list.begin(); iter != g.edge_list.end(); iter++){
			edge = (*iter);		
			to_erase.emplace(std::make_pair(index, std::vector<std::string>()));
			if(fast_edge_translate){
				for(auto const l: edge.second){
					if (stdout>3) std::cout<<"Edge "<<(*iter).first.first<<","<<(*iter).first.second<<": "<<std::endl;
					if (stdout>3) std::cout<<encoded_attributes.at("edge_attr").at(l.first).at(l.second)<<std::endl;
					if(encoded_attributes.at("edge_attr").at(l.first).at(l.second)=="dummy"){						
						//(*iter).second.erase(l.first);
						to_erase.at(index).emplace_back(l.first);
					}
					else{
						(*iter).second.at(l.first) = encoded_attributes.at("edge_attr").at(l.first).at(l.second);	
					}
					
				}
			}
			else{
				for(const auto dict : encoded_attributes.at("edge_attr")){					
					if((*iter).second.count(dict.first)>0){
						if (stdout>3) std::cout<<"Edge "<<(*iter).first.first<<","<<(*iter).first.second<<": "<<dict.first<<" -> "<<dict.second.at((*iter).second.at(dict.first))<<std::endl;
						(*iter).second.at(dict.first) = dict.second.at((*iter).second.at(dict.first));
					}
					else{
						// Add default value
						if (stdout >3) std::cout<<"Edge "<<(*iter).first.first<<","<<(*iter).first.second<<": "<<dict.first<<" -> default"<<std::endl;
						(*iter).second.emplace(std::make_pair(dict.first,
							dict.second.at("dummy")));
					}
				}
			}
			index++;

		}

		for(const auto & l : to_erase){
			iter = g.edge_list.begin();
			std::advance(iter, l.first);
			for(const auto & a : l.second){
				if (stdout >3) std::cout<<"Erasing label: Edge "<<(*iter).first.first<<" - "<<(*iter).first.second<<", "<<a<<std::endl;
				(*iter).second.erase(a);
			}
		}	

		env_coded.load_exchange_graph(g, ged::undefined(), ged::Options::ExchangeGraphType::EDGE_LIST, env_orig.get_graph_name(i), env_orig.get_graph_class(i));
		progress.increment();
		if (stdout >0) std::cout << "\rTranslating graphs: " << progress << std::flush;
	}
	if (stdout >0) std::cout << std::endl;
}


template<class T>
std::map<std::string, std::vector<std::string>> 
GED_ABC::
get_ordered_attributes(std::map<std::string, std::map<std::string, T>> &container){
	std::vector<std::string> node_attr;
	std::vector<std::string> edge_attr;
	for(const auto attr : container.at("node_attr")){
		node_attr.emplace_back(attr.first);
	}
	for(const auto attr : container.at("edge_attr")){
		edge_attr.emplace_back(attr.first);
	}
	std::sort(node_attr.begin(), node_attr.end());
	std::sort(edge_attr.begin(), edge_attr.end());

	std::map<std::string, std::vector<std::string>> res;
	res.emplace(std::make_pair("node_attr", node_attr));
	res.emplace(std::make_pair("edge_attr", edge_attr));
	return res;
}


ged::NodeMap 
GED_ABC::
get_aux_node_map(std::size_t v1, std::size_t v2, vector<ged::GEDGraph::NodeID> &v_d, vector<ged::GEDGraph::NodeID> &v_i, int stdout){
	
	if(stdout>1){
		std::cout<<"get_aux_node_map: "<<v1 <<" to "<<v2<<std::endl;
		std::cout<<"v_d: "<<std::endl;
		for(std::size_t i=0; i<v_d.size(); i++){
			std::cout<<v_d.at(i)<<", ";
		}
		std::cout<<std::endl;

		std::cout<<"v_i: "<<std::endl;
		for(std::size_t i=0; i<v_i.size(); i++){
			std::cout<<v_i.at(i)<<", ";
		}
		std::cout<<std::endl;
	}
	

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
		np.add_assignment(ged::GEDGraph::dummy_node(), i + v1 - v_d.size());
		//np.add_assignment(ged::GEDGraph::dummy_node(), i + v1 );
	}

	return np;

}

ged::NodeMap 
GED_ABC::
get_id_node_map(std::size_t  num_nodes, ged::NodeMap node_map, ged::NodeMap node_map_aux, ged::NodeMap node_map_id, int stdout){
	
	if(stdout>1){
		std::cout<<"get_id_node_map"<<std::endl;
		std::cout<<"node_map"<<std::endl;
		std::cout<<node_map<<std::endl;
		std::cout<<"node_map_aux"<<std::endl;
		std::cout<<node_map_aux<<std::endl;
		std::cout<<"node_map_id"<<std::endl;
		std::cout<<node_map_id<<std::endl;
	}
	
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
	// v_i must be sorted (?)
	for(std::size_t  n = 0; n<v_i.size(); n++){
		res.add_assignment(v_i.at(n), cont);
		cont++;

	}
	
	if(stdout>1){
		std::cout<<"res"<<std::endl;
		std::cout<<res<<std::endl;
	}
	return res;
}

void 
GED_ABC::
permute_nodes(ged::NodeMap permutation, vector<ged::GEDGraph::NodeID> &v){
	for(std::size_t i = 0; i < v.size();i ++){
		v.at(i) = permutation.image(v.at(i));
	}
}


std::vector<std::size_t> 
GED_ABC::
random_sample(std::vector<std::size_t> &population, std::size_t k, std::size_t curr_graph){
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<int> distr(0, population.size()-1);
	std::vector<std::size_t> sample;
	std::size_t x;
	while(sample.size()<k){
		x = distr(gen);
		if(x != curr_graph && std::find(sample.begin(), sample.end(), population.at(x)) == sample.end()){
			sample.emplace_back(population.at(x));
		}
	}
	std::sort(sample.begin(), sample.end());
	return sample;
}

void 
GED_ABC::
modified_costs(std::vector<double> & comp_costs, std::size_t v1, std::size_t v2, double &c_nd, double &c_ni, double &c_ns, double &c_ed, double &c_ei, double &c_es, double &c_es_id, std::size_t  &omega, std::size_t &b_ni,
	std::size_t & b_na,
	std::size_t & b_ei, 
	std::size_t & b_ea){
	// Compression costs: second version (eq 34 - 38)
	comp_costs.clear();
	if(v1 > v2){
		c_nd = 0;		
	}
	else{
		c_nd = omega;
	}
	if(v1 < v2){
		c_ni = 0;
		c_ei = 0;
		c_es = 0;
		c_es_id = -2*b_ni - b_ea;

	}
	else{
		c_ni = omega;
		c_ei = 2*b_ni + b_ea;
		c_es = 2*b_ni + b_ea;	
		c_es_id = 0;	
	}

	c_ns = b_ni + b_na;
	c_ed = 2*b_ni;

	comp_costs.emplace_back(c_ni);
	comp_costs.emplace_back(c_nd);
	comp_costs.emplace_back(c_ns);
	comp_costs.emplace_back(c_ei);
	comp_costs.emplace_back(c_ed);
	comp_costs.emplace_back(c_es);
	comp_costs.emplace_back(c_es_id);
								
}


void 
GED_ABC::
write_matrix(std::string path, MSA_di_unipi_it::MSArbor::CRow upper_bounds, std::size_t lines){
	std::ofstream file(path.c_str());
	if(file.is_open()){
		for(std::size_t i=0; i< lines; i++){
			for(std::size_t j=0; j< lines-1; j++){
				file<<upper_bounds[i + lines*j];
				if(j == lines-2){
					file<<"\n";
				}
				else{
					file<<",";
				}
			}
		}
		file.close();
	}
	else{
		std::string msg;
		msg = "Unable to open " + path + " for matrix writing";
		throw(compression_exception(msg));
	}
}

std::map<std::string, std::map<std::string, std::map<std::string,std::string>>>
GED_ABC::
get_attribute_encoding(std::map<std::string, std::map<std::string, std::set<std::string>>> &alphabets){

	std::map<std::string, std::map<std::string, std::map<std::string,std::string>>> encoded_attributes;
	encoded_attributes.emplace("node_attr", std::map<std::string, std::map<std::string,std::string>> ());
	encoded_attributes.emplace("edge_attr", std::map<std::string, std::map<std::string,std::string>> ());


	std::size_t cont;
	for(const auto & attr: alphabets.at("node_attr")){
		encoded_attributes.at("node_attr").emplace(std::make_pair(attr.first, std::map<std::string,std::string> ()));
		cont=0;
		for(const auto val : attr.second){			
			encoded_attributes.at("node_attr").at(attr.first).emplace(val, std::to_string(cont));
			cont++;
		}			
		encoded_attributes.at("node_attr").at(attr.first).emplace("dummy", std::to_string(cont));
	}

	for(const auto & attr: alphabets.at("edge_attr")){
		encoded_attributes.at("edge_attr").emplace(std::make_pair(attr.first, std::map<std::string,std::string> ()));
		cont=0;
		for(const auto val : attr.second){			
			encoded_attributes.at("edge_attr").at(attr.first).emplace(val, std::to_string(cont));
			cont++;
		}
		encoded_attributes.at("edge_attr").at(attr.first).emplace("dummy", std::to_string(cont));
	}

	return encoded_attributes;
}

std::size_t
GED_ABC::
get_file_size(std::string path){
	// Taken from StackOverflow. Not working
    struct stat stat_buf;
    int rc = stat(path.c_str(), &stat_buf);
    return rc == 0 ? stat_buf.st_size : 0;
}


void
GED_ABC::
read_xml_graph_collection(const std::string & file, std::list<std::pair<std::string, std::string>> & graphs){

	graphs.clear();

	boost::property_tree::ptree root;
	try {
		read_xml(file, root);
	}
	catch (const boost::property_tree::xml_parser_error & error) {
		throw Error(std::string("Error reading file ") + file + ": " + error.message() + ".");
	}
	// first sanity checks
	if (root.count("GraphCollection") == 0) {
		throw Error("The file " + file + " has the wrong format: no xml-element <GraphCollection>.");
	}
	if (root.count("GraphCollection") >= 2) {
		throw Error("The file " + file + " has the wrong format: more than one xml-element <GraphCollection>.");
	}
	root = root.get_child("GraphCollection");


	// Read the listed .gxl-files into the environment.
	std::vector<GEDGraph::GraphID> graph_ids;
	std::string gxl_file("");
	std::string graph_class("");
	for (const boost::property_tree::ptree::value_type & val : root) {
		if (val.first == "graph") {
			try {
				gxl_file = val.second.get<std::string>("<xmlattr>.file");
			}
			catch (const boost::property_tree::ptree_bad_path & error) {
				throw Error("The file " + file + " has the wrong format: missing xml-attribute \"file\" in element <GraphCollection>.<graph>");
			}
			catch (const boost::property_tree::ptree_bad_data & error) {
				throw Error("The file " + file + " has the wrong format: corrupted content in xml-attribute \"file\" of element <GraphCollection>.<graph>");
			}
			try {
				graph_class = val.second.get<std::string>("<xmlattr>.class");
			}
			catch (const boost::property_tree::ptree_bad_path & error) {
				throw Error("The file " + file + " has the wrong format: missing xml-attribute \"class\" in element <GraphCollection>.<graph>");
			}
			catch (const boost::property_tree::ptree_bad_data & error) {
				throw Error("The file " + file + " has the wrong format: corrupted content in xml-attribute \"class\" of element <GraphCollection>.<graph>");
			}
			graphs.push_back(std::make_pair(gxl_file, graph_class));	
		}
		else if (val.first != "<xmlattr>") {
			throw Error("The file " + file + " has the wrong format: unexpected element <GraphCollection>.<" + val.first + ">.");
		}
		
	}
}

bool
GED_ABC::
compare_graph_names(std::pair<std::string, std::string> p1, std::pair<std::string, std::string> p2){
	return(p1.first<=p2.first);
}


template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void 
GED_ABC::
create_info_file(bool binary, std::string path, std::string file_preffix,
	ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env, 
	std::map<std::string, std::map<std::string, std::set<std::string>>> &alphabets,
	std::map<std::string, std::map<std::string, std::size_t>> &attr_sizes,
	std::map<std::string, std::map<std::string, char>> &attr_types,
	std::size_t b_ni, std::size_t b_ei,
	std::vector<std::size_t> &arborescence,
	std::size_t root,
	bool fast_node_translate, bool fast_edge_translate,
	int stdout,	char separator_1, char separator_2
	){


	if(stdout>1) std::cout<<"create_info_file: start"<<std::endl;
	std::ofstream ofile;
	if(binary){
		ofile.open(path + "/_000_" + file_preffix + ".met", ios::out | ios::binary);	
	}
	else{
		ofile.open(path + "/_000_" + file_preffix + ".met", ios::out);		
	}
	
	
	std::string aux_string;
	int aux_int, aux_int_2;
	float aux_float;
	double aux_double;
	char type;
	std::vector<std::string> aux_map_labels = {"node_attr", "edge_attr"};

	std::pair<ged::GEDGraph::GraphID, ged::GEDGraph::GraphID> limits;
	limits = env.graph_ids();

	// info_file
	if(ofile.is_open()){

		// graph classes 
		if(stdout>1) std::cout<<"create_info_file: graph_classes"<<std::endl;
		type = attr_types.at("graph_attr").at("class");
		if(binary){
			ofile.put(type);
		}
		for(ged::GEDGraph::GraphID i{limits.first}; i<limits.second; i++){
			if(i == root) continue;	
			if(binary){				
				aux_string = env.get_graph_class(i);
				switch(type){ 
			        case 'i': // int
			        	aux_int = std::stoi(aux_string);
			            ofile.write(reinterpret_cast<char*>(&aux_int), sizeof(int));
			            break;
			        case 'f': // float
			        	aux_float = std::stof(aux_string);
			            ofile.write(reinterpret_cast<char*>(&aux_float), sizeof(float));
			            break;
			        case 'd': // double
			        	aux_double = std::stod(aux_string);
			            ofile.write(reinterpret_cast<char*>(&aux_double), sizeof(double));
			            break;
			        case 'c': // char   
			            ofile.put(aux_string[0]);
			            break;
			        case 's': // string   
			            write_word_binary(ofile, aux_string);
			            ofile.put(separator_1);
			            break;
			    }

			}	
			else{
				ofile<<env.get_graph_class(i)<<separator_1;	
			}	
		}

		// arborescence
		if(stdout>1) std::cout<<"create_info_file: Arborescence"<<std::endl;
		if(binary){
			// Number of bytes to read each time			
			aux_int = bytes(arborescence.size());
			if(stdout>1) std::cout<<aux_int<<std::endl;
			write_size_t_to_chars(ofile, aux_int, 1);						
			aux_int = bytes(arborescence.size());
		}
		for(std::size_t k=0; k<arborescence.size();k++){
			if(binary){
				write_size_t_to_chars(ofile, arborescence.at(k), aux_int);
			}
			else{
				ofile<<arborescence.at(k)<<separator_1;
			}
			
		}

		// # of attributes
		// for nodes and edges
		// b_ni and b_ei
		aux_int = alphabets.at("node_attr").size();
		aux_int_2 = alphabets.at("edge_attr").size();
		if(binary){			
			write_size_t_to_chars(ofile, aux_int, 1);
			write_size_t_to_chars(ofile, aux_int_2, 1);
			write_size_t_to_chars(ofile, b_ni, 1);
			write_size_t_to_chars(ofile, b_ei, 1);
		}
		else{
			ofile<<aux_int<<separator_1;
			ofile<<aux_int_2<<separator_1;

			ofile<<b_ni<<separator_1;
			ofile<<b_ei<<separator_1;
		}
		

		std::map<std::string, std::vector<std::string>> ordered_attributes = get_ordered_attributes<std::set<std::string>>(alphabets); 
		// Size in bytes of attributes
		for(const auto & node_edge : aux_map_labels){
			for(const auto attr: ordered_attributes.at(node_edge)){
				if(binary){
					write_size_t_to_chars(ofile, attr_sizes.at(node_edge).at(attr), 1);		
				}
				else{
					ofile<< attr_sizes.at(node_edge).at(attr)<<separator_1;	
				}
			}
		}

		// Alphabets

		for(const auto & node_edge : aux_map_labels){
			if(stdout>1) std::cout<<"create_info_file: alphabets, "<<node_edge<<std::endl;
			for(const auto attr: ordered_attributes.at(node_edge)){
				if(binary){
					// The name followed by a separator
					if(stdout>3) std::cout<<attr<<std::endl;
					write_word_binary(ofile, attr);
					ofile.put(separator_1);
					// type of attributes
					if(attr_types.at(node_edge).count(attr)==0) {
						aux_string = "The attribute " + attr + " is not specified in attr_types." + node_edge; 
						throw compression_exception( aux_string );
					}
					type = attr_types.at(node_edge).at(attr);
					if(stdout>3) std::cout<<type<<std::endl;
					ofile.put(type);
					// Number of values to read in 4 bytes (max 4 294 967 296 values)
					write_size_t_to_chars(ofile, alphabets.at(node_edge).at(attr).size(), 4);

					if(type == 's'){
						// Worst case, must use separators
						for(auto const val : alphabets.at(node_edge).at(attr)){							
							write_word_binary(ofile, val);
							if(stdout>3) std::cout<<val<<std::endl;
							ofile.put(separator_1);
						}
					}
					else{
						// the values using fixed size
						for(auto const val : alphabets.at(node_edge).at(attr)){
							if(stdout>3) std::cout<<val<<std::endl;
							switch(type){ 
					        case 'i': // int
					        	aux_int = std::stoi(val);
					            ofile.write(reinterpret_cast<char*>(&aux_int), sizeof(int));
					            break;
					        case 'f': // float
					        	aux_float = std::stof(val);
					            ofile.write(reinterpret_cast<char*>(&aux_float), sizeof(float));
					            break;
					        case 'd': // double
					        	aux_double = std::stod(val);
					            ofile.write(reinterpret_cast<char*>(&aux_double), sizeof(double));
					            break;
					        case 'c': // char   
					            ofile.put(val[0]);
					            break;				        
					    	}				
						}	
					}
					// End of attribute
				}
				else{
					// attr name
					ofile<<attr<<separator_1;
					// attr size (# of elements)
					ofile<<alphabets.at(node_edge).at(attr).size()<<separator_1;					
					for(auto const val : alphabets.at(node_edge).at(attr)){
						ofile<<val<<separator_1;				
					}
				}
							
			}

		}		

		ofile.close();	
	}
	else{
		throw compression_exception( "Error when opening meta data file" );

	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void 
GED_ABC::
encode_single_graph(bool binary, bool relaxed, std::string path,
	ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env_coded, 
	std::size_t parent_num, std::size_t child,
	std::map<ged::GEDGraph::GraphID, ged::NodeMap> &graph_permutations,
	std::map<std::string, std::vector<std::string>> &ordered_attributes,
	std::map<std::string, std::map<std::string, std::size_t>> &attr_sizes,
	std::size_t b_ni, std::size_t b_ei,
	int stdout, char separator
	){
	std::ofstream ofile;

	ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> g1, g2;
	std::vector<ged::GEDGraph::NodeID> v_d, v_i, v_s, v_is, v_i_aux;
	std::vector<UserNodeLabel> varphi_i, varphi_s, varphi_i_aux;
	std::vector<std::pair<ged::GEDGraph::NodeID,ged::GEDGraph::NodeID>> e_nd, e_ed, e_ni, e_ei, e_s, e_is;
	std::vector<UserEdgeLabel> phi_ni, phi_ei, phi_s;
	std::map<std::size_t,std::size_t> v_i_map;

	std::size_t from, to;

	ged::NodeMap node_map_aux = ged::NodeMap(1,1);
	ged::NodeMap node_map_id = ged::NodeMap(1,1);
	ged::NodeMap node_map_id_before = ged::NodeMap(1,1);

	if(stdout>1) std::cout<<"Get NodeMap"<<std::endl;
	ged::NodeMap node_map = env_coded.get_node_map(parent_num, child);
	if(stdout>1) std::cout<<"Get g1"<<std::endl;
	g1 = env_coded.get_graph(parent_num, true, false, true); // edge list and adj matrix
	if(stdout>1) std::cout<<"Get g2"<<std::endl;
	g2 = env_coded.get_graph(child, true, false, true); // edge list and adj matrix


	if (stdout>1) {
		std::cout<<std::endl<<"START ENCODING OF GRAPHS-------------------------"<<std::endl;
		std::cout<<parent_num<<" -> "<<child<<std::endl<< node_map<<std::endl<<std::endl;
	}

	if(stdout>1) std::cout<<"Get compression sets Firts time"<<std::endl;

	get_all_compression_sets<UserNodeID, UserNodeLabel, UserEdgeLabel>(node_map,v_d,v_i,varphi_i,v_s,v_is,varphi_s,e_nd,e_ed,e_ni,e_ei,
		phi_ni,phi_ei,e_s,e_is,phi_s,g1,g2);

	if (stdout>3) print_compression_sets(v_d,v_i,v_s,v_is,e_nd,e_ed,e_ni,e_ei,e_s,e_is, g1, g2);

	node_map_id_before = graph_permutations.at(parent_num);
	node_map_aux.clear();

	if(stdout>1) std::cout<<"Get aux node map"<<std::endl;
	permute_nodes(node_map_id_before, v_d); // Now in the compressed version (g1 tilda)
	if(stdout>1) std::cout<<"*** v_d transformed"<<std::endl;
	permute_nodes(node_map_id_before, v_s); // Now in the compressed version (g1 tilda)
	if(stdout>1) std::cout<<"*** v_s transformed"<<std::endl;
	permute_nodes(node_map_id_before, v_is); // Now in the compressed version (g1 tilda)
	if(stdout>1) std::cout<<"*** v_is transformed"<<std::endl;

	if(stdout>1) std::cout<<"With modified vertex"<<std::endl;
	if (stdout>3) print_compression_sets(v_d,v_i,v_s,v_is,e_nd,e_ed,e_ni,e_ei,e_s,e_is, g1, g2);

	if(stdout>1) std::cout<<"Get aux node map"<<std::endl;
	// NodeMap between the compressed versions
	node_map_aux = get_aux_node_map(g1.num_nodes, g2.num_nodes, v_d, v_i, stdout); 

	if(stdout>1) std::cout<<"Get id node map"<<std::endl;
	node_map_id = get_id_node_map(g2.num_nodes, node_map, node_map_aux, node_map_id_before, stdout);
	graph_permutations.emplace(std::make_pair(child, node_map_id));

	
	// Check if deletions and insertions are mutually exclusive
	if(stdout>0 && v_d.size()>0 && v_i.size()>0){
		std::cout<<"WARNING: Vertex deletions and insertions at the same time"<<std::endl;	
	}
	
	std::string graph_file_name = path + "/" + env_coded.get_graph_name(child) + ".cmp";
	if(stdout>1) std::cout<< "File: " << graph_file_name << std::endl;
	

	if(binary){
		ofile.open(graph_file_name, ios::out | ios::binary);
	}
	else{
		ofile.open(graph_file_name, ios::out);	
	}

	if(! ofile.is_open()){
		
		std::string msg = "Unable to open file " + graph_file_name;
		throw compression_exception(msg);
		
	}	
	//Rec V
	if(stdout>1) std::cout<<"REC V"<<std::endl;
	if(binary){
		write_size_t_to_chars(ofile, g2.num_nodes, b_ni);
		write_size_t_to_chars(ofile, v_s.size(), b_ni);
		if(relaxed) write_size_t_to_chars(ofile, v_is.size(), b_ni);
	}
	else{
		ofile<<g2.num_nodes<<separator;
		ofile<<v_s.size()<<separator;
		if(relaxed) ofile<<v_is.size()<<separator;
	}
		
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
	if(!relaxed && (g1.num_nodes < g2.num_nodes)){
						
		// Node Insertions
		for(std::size_t i=0; i<v_i.size(); i++){
			// Already in g2 tilda			
			if(binary){
				for(const auto attr: ordered_attributes.at("node_attr")){
					write_size_t_to_chars(ofile, std::stoi(varphi_i_aux.at(i).at(attr)), attr_sizes.at("node_attr").at(attr));
				}
			}
			else{
				for(const auto attr: ordered_attributes.at("node_attr")){
					ofile<<varphi_i_aux.at(i).at(attr)<<separator;
				}	
			}											
		}	
	}
	else{
		if(v_d.size() < v_is.size()){
			// Node Deletions				
			if(binary){
				for(std::size_t i=0; i<v_d.size(); i++){					
					// Nodes in g1 tilda
					write_size_t_to_chars(ofile, v_d.at(i), b_ni);
				}	
			}
			else{
				for(std::size_t i=0; i<v_d.size(); i++){					
					// Nodes in g1 tilda
					ofile<<v_d.at(i)<<separator;
				}		
			}
			
		}
		else{
			// Node Identical Substitutions
			if(binary){
				for(std::size_t i=0; i<v_is.size(); i++){
					// Nodes in g1 tilda
					write_size_t_to_chars(ofile, v_is.at(i), b_ni);				
				}		
			}
			else{
				for(std::size_t i=0; i<v_is.size(); i++){
					// Nodes in g1 tilda
					ofile<< v_is.at(i)<<separator;
				}	
			}

		}
	}

	// Node Substitutions				
	for(std::size_t i=0; i<v_s.size(); i++){
		// Nodes already in g1 tilda
		if(binary){
			write_size_t_to_chars(ofile, v_s.at(i), b_ni);
			for(const auto attr: ordered_attributes.at("node_attr")){
				write_size_t_to_chars(ofile, std::stoi(varphi_s.at(i).at(attr)), attr_sizes.at("node_attr").at(attr));
			}
		}
		else{
			ofile << v_s.at(i) << separator;
			for(const auto attr: ordered_attributes.at("node_attr")){
				ofile << varphi_s.at(i).at(attr) << separator;
			}
		}
	}

	if(relaxed){
		// Node insertions
		if(binary){
			for(std::size_t i=0; i<v_i.size(); i++){
				for(const auto attr: ordered_attributes.at("node_attr")){
					write_size_t_to_chars(ofile, std::stoi(varphi_i_aux.at(i).at(attr)), attr_sizes.at("node_attr").at(attr));
				}
			}
		}
		else{
			for(std::size_t i=0; i<v_i.size(); i++){
				for(const auto attr: ordered_attributes.at("node_attr")){
					ofile<<varphi_i_aux.at(i).at(attr)<<separator;
				}	
			}
		}
	}

	//Rec E
							
	std::pair<ged::GEDGraph::NodeID,ged::GEDGraph::NodeID> e;
	
	if(binary){
		write_size_t_to_chars(ofile, g2.num_edges, b_ei);
		write_size_t_to_chars(ofile, e_s.size(), b_ei);
	}
	else{
		ofile << g2.num_edges << separator;
		ofile << e_s.size() << separator;
	}
		


	if(e_ed.size() <= e_is.size()){
		// Write e_d (edge deletions)
		if(binary){
			write_size_t_to_chars(ofile, 0, 1);
			write_size_t_to_chars(ofile, e_ed.size(), b_ei);			
			for(auto e : e_ed){
				// Nodes in g2 tilda
				from = node_map_aux.image(node_map_id_before.image(e.first));
				to = node_map_aux.image(node_map_id_before.image(e.second));
				write_size_t_to_chars(ofile, from, b_ni);
				write_size_t_to_chars(ofile, to, b_ni);
			}	
		}
		else{
			ofile << 0 << separator;
			ofile << e_ed.size() << separator;
			for(auto e : e_ed){
				// Nodes in g2 tilda
				from = node_map_aux.image(node_map_id_before.image(e.first));
				to = node_map_aux.image(node_map_id_before.image(e.second));
				ofile << from << separator;
				ofile << to << separator;
			}
		}
		
	}
	else{
		// Write e_is (edge identical substitutions)
		if(binary){
			write_size_t_to_chars(ofile, 1, 1);
			write_size_t_to_chars(ofile, e_is.size(), b_ei);
			for(auto e : e_is){
				// Nodes in g2 tilda
				from = node_map_aux.image(node_map_id_before.image(e.first));
				to = node_map_aux.image(node_map_id_before.image(e.second));
				write_size_t_to_chars(ofile, from, b_ni);
				write_size_t_to_chars(ofile, to, b_ni);
			}	
		}
		else{
			ofile << 1 << separator;
			ofile << e_is.size() << separator;
			for(auto e : e_is){
				// Nodes in g2 tilda
				from = node_map_aux.image(node_map_id_before.image(e.first));
				to = node_map_aux.image(node_map_id_before.image(e.second));
				ofile << from << separator;
				ofile << to << separator;
			}
		}
	}

	
	// Edge Substitutions (Non identical)
	if(binary){
		for(std::size_t i = 0; i<e_s.size(); i++){
			e = e_s.at(i);
			// Nodes in g2 tilda
			from = node_map_aux.image(node_map_id_before.image(e.first));
			to = node_map_aux.image(node_map_id_before.image(e.second));
			write_size_t_to_chars(ofile, from, b_ni);
			write_size_t_to_chars(ofile, to, b_ni);
			for(const auto attr: ordered_attributes.at("edge_attr")){
				write_size_t_to_chars(ofile, std::stoi(phi_s.at(i).at(attr)), attr_sizes.at("edge_attr").at(attr));
			}					
		}	
	}
	else{
		for(std::size_t i = 0; i<e_s.size(); i++){
			e = e_s.at(i);
			// Nodes in g2 tilda
			from = node_map_aux.image(node_map_id_before.image(e.first));
			to = node_map_aux.image(node_map_id_before.image(e.second));
			ofile << from << separator;
			ofile << to << separator;
			for(const auto attr: ordered_attributes.at("edge_attr")){
				ofile <<phi_s.at(i).at(attr) << separator;
			}					
		}
	}	
	

	// Edge Insertions
	if(binary){
		for(std::size_t i = 0; i<e_ni.size(); i++){
			e = e_ni.at(i);
			// Nodes in g2 tilda
			from = node_map_id.image(e.first);
			to = node_map_id.image(e.second);				

			write_size_t_to_chars(ofile, from, b_ni);
			write_size_t_to_chars(ofile, to, b_ni);
			for(const auto attr: ordered_attributes.at("edge_attr")){
				write_size_t_to_chars(ofile, std::stoi(phi_ni.at(i).at(attr)), attr_sizes.at("edge_attr").at(attr));
			}	
			
		}

		for(std::size_t i = 0; i<e_ei.size(); i++){
			e = e_ei.at(i);		
			// Nodes in g2 tilda
			from = node_map_id.image(e.first);
			to = node_map_id.image(e.second);
			
			write_size_t_to_chars(ofile, from, b_ni);
			write_size_t_to_chars(ofile, to, b_ni);
			for(const auto attr: ordered_attributes.at("edge_attr")){
				write_size_t_to_chars(ofile, std::stoi(phi_ei.at(i).at(attr)), attr_sizes.at("edge_attr").at(attr));
			}	
		}
	}
	else{
		for(std::size_t i = 0; i<e_ni.size(); i++){
			e = e_ni.at(i);
			// Nodes in g2 tilda
			from = node_map_id.image(e.first);
			to = node_map_id.image(e.second);				

			ofile << from << separator;
			ofile << to << separator;
			for(const auto attr: ordered_attributes.at("edge_attr")){
				ofile << phi_ni.at(i).at(attr)  << separator;
			}	
			
		}

		for(std::size_t i = 0; i<e_ei.size(); i++){
			e = e_ei.at(i);		
			// Nodes in g2 tilda
			from = node_map_id.image(e.first);
			to = node_map_id.image(e.second);
			
			ofile << from << separator;
			ofile << to << separator;
			for(const auto attr: ordered_attributes.at("edge_attr")){
				ofile << phi_ei.at(i).at(attr)  << separator;
			}	
		}
	}


	ofile.close();
	if(stdout>2) std::cout<<"Wrote: "<<graph_file_name<<std::endl;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void 
GED_ABC::
encode_single_graph_relaxed(bool binary, std::string path,
	ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env_coded, 
	std::size_t parent_num, std::size_t child,
	std::map<ged::GEDGraph::GraphID, ged::NodeMap> &graph_permutations,
	std::map<std::string, std::vector<std::string>> &ordered_attributes,
	std::map<std::string, std::map<std::string, std::size_t>> &attr_sizes,
	std::size_t b_ni, std::size_t b_ei,
	int stdout, char separator
	){
	std::ofstream ofile;

	ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> g1, g2;
	std::vector<ged::GEDGraph::NodeID> v_d, v_i, v_s, v_is, v_i_aux;
	std::vector<UserNodeLabel> varphi_i, varphi_s, varphi_i_aux;
	std::vector<std::pair<ged::GEDGraph::NodeID,ged::GEDGraph::NodeID>> e_nd, e_ed, e_ni, e_ei, e_s, e_is;
	std::vector<UserEdgeLabel> phi_ni, phi_ei, phi_s;
	std::map<std::size_t,std::size_t> v_i_map;

	std::size_t from, to;

	ged::NodeMap node_map_aux = ged::NodeMap(1,1);
	ged::NodeMap node_map_id = ged::NodeMap(1,1);
	ged::NodeMap node_map_id_before = ged::NodeMap(1,1);

	if(stdout>1) std::cout<<"Get NodeMap"<<std::endl;
	ged::NodeMap node_map = env_coded.get_node_map(parent_num, child);
	if(stdout>1) std::cout<<"Get g1"<<std::endl;
	g1 = env_coded.get_graph(parent_num, true, false, true); // edge list and adj matrix
	if(stdout>1) std::cout<<"Get g2"<<std::endl;
	g2 = env_coded.get_graph(child, true, false, true); // edge list and adj matrix


	if (stdout>1) {
		std::cout<<std::endl<<"START ENCODING OF GRAPHS-------------------------"<<std::endl;
		std::cout<<parent_num<<" -> "<<child<<std::endl<< node_map<<std::endl<<std::endl;
	}

	if(stdout>1) std::cout<<"Get compression sets Firts time"<<std::endl;

	get_all_compression_sets<UserNodeID, UserNodeLabel, UserEdgeLabel>(node_map,v_d,v_i,varphi_i,v_s,v_is,varphi_s,e_nd,e_ed,e_ni,e_ei,
		phi_ni,phi_ei,e_s,e_is,phi_s,g1,g2);
	if (stdout>3) print_compression_sets(v_d,v_i,v_s,v_is,e_nd,e_ed,e_ni,e_ei,e_s,e_is, g1, g2);

	node_map_id_before = graph_permutations.at(parent_num);
	node_map_aux.clear();

	if(stdout>1) std::cout<<"Get aux node map"<<std::endl;
	permute_nodes(node_map_id_before, v_d); // Now in the compressed version (g1 tilda)
	if(stdout>1) std::cout<<"*** v_d transformed"<<std::endl;
	permute_nodes(node_map_id_before, v_s); // Now in the compressed version (g1 tilda)
	if(stdout>1) std::cout<<"*** v_s transformed"<<std::endl;
	permute_nodes(node_map_id_before, v_is); // Now in the compressed version (g1 tilda)
	if(stdout>1) std::cout<<"*** v_is transformed"<<std::endl;

	if(stdout>1) std::cout<<"With modified vertex"<<std::endl;
	if (stdout>3) print_compression_sets(v_d,v_i,v_s,v_is,e_nd,e_ed,e_ni,e_ei,e_s,e_is, g1, g2);

	if(stdout>1) std::cout<<"Get aux node map"<<std::endl;
	// NodeMap between the compressed versions
	node_map_aux = get_aux_node_map(g1.num_nodes, g2.num_nodes, v_d, v_i, stdout); 

	if(stdout>1) std::cout<<"Get id node map"<<std::endl;
	node_map_id = get_id_node_map(g2.num_nodes, node_map, node_map_aux, node_map_id_before, stdout);
	graph_permutations.emplace(std::make_pair(child, node_map_id));

	
	// Check if deletions and insertions are mutually exclusive
	if(stdout>0 && v_d.size()>0 && v_i.size()>0) std::cout<<"WARNING: Vertex deletions and insertions at the same time"<<std::endl;
	if(stdout>1 && (e_nd.size()+e_ed.size())>0 && (e_ni.size()+e_ei.size())>0) std::cout<<"WARNING: Edge deletions and insertions at the same time"<<std::endl;

	std::string graph_file_name = path + "/" + env_coded.get_graph_name(child) + ".cmp";
	if(stdout>1) std::cout<< "File: " << graph_file_name << std::endl;

	ofile.open(graph_file_name, ios::out | ios::binary);
	if(! ofile.is_open()){
		
		std::string msg = "Unable to open file " + graph_file_name;
		throw compression_exception(msg);
		
	}
	//Eq 46 and 47
	//Rec V
	if(stdout>1) std::cout<<"REC V"<<std::endl;			
	write_size_t_to_chars(ofile, g2.num_nodes, b_ni);
	write_size_t_to_chars(ofile, v_s.size(), b_ni);
	// (new in relaxed model)
	write_size_t_to_chars(ofile, v_is.size(), b_ni);

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
		v_i_aux.emplace_back(k + g1.num_nodes -v_d.size());
		varphi_i_aux.emplace_back(varphi_i.at(k));
		v_i_map.emplace(v_i.at(k), k + g1.num_nodes -v_d.size());
	}

	if(v_d.size() < v_is.size()){
		// Node Deletions (type_read 0)								
		for(std::size_t i=0; i<v_d.size(); i++){	
			// Nodes in g1 tilda				
			write_size_t_to_chars(ofile, v_d.at(i), b_ni);				
		}
	}
	else{
		// Node Identical Substitutions (type_read 1)						
		for(std::size_t i=0; i<v_is.size(); i++){
			// Nodes in g1 tilda
			write_size_t_to_chars(ofile, v_is.at(i), b_ni);
			
		}					
	}

	// Node Substitutions	
	for(std::size_t i=0; i<v_s.size(); i++){
		// Nodes in g1 tilda
		write_size_t_to_chars(ofile, v_s.at(i), b_ni);
		for(const auto attr: ordered_attributes.at("node_attr")){
			write_size_t_to_chars(ofile, std::stoi(varphi_s.at(i).at(attr)), attr_sizes.at("node_attr").at(attr));
		}	
	}

	// Node Insertions
	for(std::size_t i=0; i<v_i.size(); i++){
		// Nodes in g2 tilda
		for(const auto attr: ordered_attributes.at("node_attr")){
				write_size_t_to_chars(ofile, std::stoi(varphi_i_aux.at(i).at(attr)), attr_sizes.at("node_attr").at(attr));
			}	
	}	

	//Rec E		
	std::pair<ged::GEDGraph::NodeID,ged::GEDGraph::NodeID> e;


	// Sorting is important as we are using the index in the edge list according to this order
	write_size_t_to_chars(ofile, g2.num_edges, b_ei);
	write_size_t_to_chars(ofile, e_s.size(), b_ei);

	
	if(e_ed.size() <= e_is.size()){
		// Write e_d (edge deletions)
		write_size_t_to_chars(ofile, 0, 1);
		write_size_t_to_chars(ofile, e_ed.size(), b_ei);	
		for(auto e : e_ed){

			// Nodes in g2 tilda
			from = node_map_aux.image(node_map_id_before.image(e.first));
			to = node_map_aux.image(node_map_id_before.image(e.second));
			write_size_t_to_chars(ofile, from, b_ni);
			write_size_t_to_chars(ofile, to, b_ni);
		}
	}
	else{
		// Write e_is (edge identical substitutions)
		write_size_t_to_chars(ofile, 1, 1);
		write_size_t_to_chars(ofile, e_is.size(), b_ei);
		for(auto e : e_is){

			// Nodes in g2 tilda
			from = node_map_aux.image(node_map_id_before.image(e.first));
			to = node_map_aux.image(node_map_id_before.image(e.second));
			write_size_t_to_chars(ofile, from, b_ni);
			write_size_t_to_chars(ofile, to, b_ni);
		}
	}

	// Edge Substitutions (Non identical)
	for(std::size_t i = 0; i<e_s.size(); i++){
		e = e_s.at(i);	
		// Nodes in g2 tilda
		from = node_map_aux.image(node_map_id_before.image(e.first));
		to = node_map_aux.image(node_map_id_before.image(e.second));
		write_size_t_to_chars(ofile, from, b_ni);
		write_size_t_to_chars(ofile, to, b_ni);
		for(const auto attr: ordered_attributes.at("edge_attr")){
			write_size_t_to_chars(ofile, std::stoi(phi_s.at(i).at(attr)), attr_sizes.at("edge_attr").at(attr));
		}	
		
	}

	// Edge Insertions
	for(std::size_t i = 0; i<e_ni.size(); i++){
		e = e_ni.at(i);	
		// Nodes in g2 tilda
		from = node_map_id.image(e.first);
		to = node_map_id.image(e.second);	
		write_size_t_to_chars(ofile, from, b_ni);
		write_size_t_to_chars(ofile, to, b_ni);
		for(const auto attr: ordered_attributes.at("edge_attr")){
			write_size_t_to_chars(ofile, std::stoi(phi_ni.at(i).at(attr)), attr_sizes.at("edge_attr").at(attr));
		}	
		
	}

	for(std::size_t i = 0; i<e_ei.size(); i++){
		e = e_ei.at(i);	
		// Nodes in g2 tilda
		from = node_map_id.image(e.first);
		to = node_map_id.image(e.second);				
		write_size_t_to_chars(ofile, from, b_ni);
		write_size_t_to_chars(ofile, to, b_ni);
		for(const auto attr: ordered_attributes.at("edge_attr")){
			write_size_t_to_chars(ofile, std::stoi(phi_ei.at(i).at(attr)), attr_sizes.at("edge_attr").at(attr));
		}	
	}
	ofile.close();
	if(stdout>2) std::cout<<"Wrote: "<<graph_file_name<<std::endl;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void 
GED_ABC::
encode_arborescence(bool binary, bool relaxed, std::string path, std::string file_preffix,
	ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env_coded,
	std::map<std::string, std::map<std::string, std::size_t>> &attr_sizes,
	std::size_t b_ni, std::size_t b_ei,
	std::vector<std::size_t> &arborescence,
	std::size_t root,
	int stdout,
	char separator
	){

	if(stdout>1) std::cout<<"encode_arborescence: START"<<std::endl;
	
	
	if(stdout>1) std::cout<<"encode_arborescence: Create children"<<std::endl;
	std::map<std::size_t, std::vector<std::size_t>> children;
	std::size_t num_nodes = arborescence.size()+1; 
	for(std::size_t i=0; i<num_nodes; i++){
		children.emplace(std::make_pair(i, std::vector<std::size_t> ()));
	}	

	for(std::size_t i=0; i<arborescence.size(); i++){
		children.at(arborescence.at(i)).emplace_back(i);			
	}

	std::map<std::string, std::vector<std::string>> ordered_attributes = get_ordered_attributes<std::size_t>(attr_sizes);

	ged::NodeMap node_map_id_before = ged::NodeMap(1,1);
	std::map<ged::GEDGraph::GraphID, ged::NodeMap> graph_permutations;
	node_map_id_before.add_assignment(0,0);
	graph_permutations.emplace(std::make_pair(root, node_map_id_before));

	// Encode each graph
	if(stdout>1) std::cout<<"Start encode arborescence"<<std::endl;
	std::list<std::size_t> to_do;
	std::size_t parent_num;
	to_do.emplace_front(root);
	while(!to_do.empty()){
		
		parent_num = to_do.front();
		to_do.pop_front();
		
		for(auto const child : children.at(parent_num)){			
			to_do.emplace_front(child);
			// Encode graph
			encode_single_graph(binary, relaxed, path,
					env_coded, parent_num, child,
					graph_permutations,ordered_attributes, attr_sizes, b_ni, b_ei, stdout, separator);				
		}
	}
}


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
void 
GED_ABC::
compress_collection(bool binary, bool relaxed, std::string graph_dir, std::string xml_file,
	Options::GXLNodeEdgeType node_type, Options::GXLNodeEdgeType edge_type,
	const std::unordered_set<std::string> & irrelevant_node_attributes, 
	const std::unordered_set<std::string> & irrelevant_edge_attributes,
	std::map<std::string, std::map<std::string, char>> & attr_types,
	std::string output_encoded,
	std::string output_other,
	std::string file_preffix,
	std::map<std::string, std::string> &args,
	int stdout,
	std::vector<std::string> &headers,
	std::vector<std::string> &values,
	char separator_1,
	char separator_2
	){

	ged::Seconds runtime;
	auto start = std::chrono::high_resolution_clock::now();

	std::list<std::pair<std::string, std::string>> graph_names_classes;

	read_xml_graph_collection(xml_file, graph_names_classes);

	graph_names_classes.sort(compare_graph_names);
	ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> env;
	
	ged::ProgressBar progress_load( graph_names_classes.size() );
	if (stdout >0) std::cout << "\rLoading graphs: " << progress_load << std::flush;
	typename std::list<std::pair<std::string, std::string>>::iterator iter;
	for(iter = graph_names_classes.begin(); iter != graph_names_classes.end(); iter ++){
		env.load_gxl_graph(graph_dir, (*iter).first, node_type, edge_type, 
			irrelevant_node_attributes, irrelevant_edge_attributes, ged::undefined(), (*iter).second);
		progress_load.increment();
		if (stdout >0) std::cout << "\rLoading graphs: " << progress_load << std::flush;
	}


	auto start_init = std::chrono::high_resolution_clock::now();
	runtime = start_init - start;
	double load_time = runtime.count();

	std::map<std::string, std::map<std::string, std::vector<std::string>>> distribution;
	std::map<std::string, std::map<std::string, std::set<std::string>>> alphabets;
	std::map<std::string, std::map<std::string, std::size_t>> attr_sizes;
	std::size_t b_ni, b_na, b_ei, b_ea;
	bool fast_node_translate=false, fast_edge_translate=false;

	std::string aux_string;
	
	// 1. Get structure
	get_graphs_structure<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>(env, distribution, alphabets, attr_sizes, b_ni, b_na, b_ei, b_ea, fast_node_translate, fast_edge_translate);	

	if(stdout>1){
		std::cout<<"Structure:"<<std::endl;
		std::cout<<"b_ni: "<<b_ni<<", b_na: "<<b_na<<std::endl;
		std::cout<<"b_ei: "<<b_ei<<", b_ea: "<<b_ea<<std::endl;

		std::cout<<"Nodes: fast_translate = "<<fast_node_translate<<std::endl;
		for(auto const& a: alphabets.at("node_attr")){
			std::cout<<"\t"<<a.first<<": "<<a.second.size()<<" values"<<std::endl;	
		}
		std::cout<<"Edges: fast_translate = "<<fast_edge_translate<<std::endl;
		for(auto const& a: alphabets.at("edge_attr")){
			std::cout<<"\t"<<a.first<<": "<<a.second.size()<<" values"<<std::endl;	
		}
	}

	std::map<std::string, std::map<std::string, std::map<std::string,std::string>>> encoded_attributes;
	encoded_attributes = get_attribute_encoding(alphabets);

	// 2. Translate
	ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> env_coded;
	ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> g1;
	ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> g2;
	

	if(stdout>1) std::cout<<"translate_env: "<<fast_node_translate<<", "<<fast_edge_translate<<std::endl;
	translate_env<UserNodeID,UserNodeLabel,UserEdgeLabel>(encoded_attributes, env, env_coded, fast_node_translate, fast_edge_translate, stdout);
	if(stdout>1) std::cout<<"end___translate_env"<<std::endl;

	// deallocate first env to save space (is this really deallocating?) The initial env is then "lost"
	env = ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>();

	ged::GEDGraph::GraphID empty_id = env_coded.add_graph("empty","");

	// 3. GED calculation

	// 3.1 Set edit costs
	std::vector<double> comp_costs;
	double c_nd, c_ni, c_ns, c_ed, c_ei, c_es, c_es_id;
	
	bool edit_modified = false;
	if(args.count("edit_cost_type")>0 && args.at("edit_cost_type")=="mod"){
		edit_modified = true;
		if(stdout>1) std::cout<<" EDIT COSTS: Modified"<<std::endl;
	}
	else{
		if(stdout>1) std::cout<<" EDIT COSTS: Simple"<<std::endl;	
	}
	// initialization (used in the not modified case)
	c_ni = b_na;
	c_nd = b_ni;
	c_ns = b_ni + b_na;
	c_ei = 2*b_ni + b_ea;
	c_ed = 2*b_ni;
	c_es = 2*b_ni + b_ea;
	c_es_id = 0;
	comp_costs.emplace_back(c_ni);
	comp_costs.emplace_back(c_nd);
	comp_costs.emplace_back(c_ns);
	comp_costs.emplace_back(c_ei);
	comp_costs.emplace_back(c_ed);
	comp_costs.emplace_back(c_es);
	comp_costs.emplace_back(c_es_id);

	env_coded.set_edit_costs(ged::Options::EditCosts::COMPRESSION, comp_costs);
	env_coded.init(ged::Options::InitType::LAZY_WITHOUT_SHUFFLED_COPIES);

	// 3.2 Set method
	std::string ged_method_options = "";
	if(args.count("ged_method_options")>0){
		ged_method_options = args.at("ged_method_options");
	}
	else{
		// Default number of threads
		ged_method_options = "16";	
	}

	// Currently working with BRANCH_UNIFORM. All other methods have not been thoroughly tested
	if(args.count("ged_method")>0){
		if(stdout) std::cout<<"GED method: "<<args.at("ged_method")<<std::endl;
		if(args.at("ged_method") == "branch_uniform"){
			env_coded.set_method(ged::Options::GEDMethod::BRANCH_UNIFORM, "--threads " + ged_method_options);
		}
		else{
			if(args.at("ged_method") == "branch_fast"){
				env_coded.set_method(ged::Options::GEDMethod::BRANCH_FAST, "--threads " + ged_method_options);
			}
			else{
				if(args.at("ged_method") == "branch"){
					env_coded.set_method(ged::Options::GEDMethod::BRANCH, "--threads " + ged_method_options);
				}
				else{
					if(args.at("ged_method") == "ipfp"){
						env_coded.set_method(ged::Options::GEDMethod::IPFP, "--threads " + ged_method_options);
					}
					else{
						std::cout<<"No valid ged_method in args. Setting it to BRANCH_UNIFORM"<<std::endl;
						env_coded.set_method(ged::Options::GEDMethod::BRANCH_UNIFORM, "--threads " + ged_method_options);
					}
				}
			}
		}
	}
	else{
		std::cout<<"No field ged_method in args. Setting it to BRANCH_UNIFORM"<<std::endl;
		env_coded.set_method(ged::Options::GEDMethod::BRANCH_UNIFORM, "--threads " + ged_method_options);
		
	}

	// 3.3 Get upper bounds

	std::pair<ged::GEDGraph::GraphID, ged::GEDGraph::GraphID> limits = env_coded.graph_ids();

	if(stdout>0) std::cout<<"Limits (no of graphs): "<<limits.first<<", "<<limits.second-1<<std::endl;
	if(stdout>0) std::cout<<"Empty graph id: "<<empty_id<<std::endl;
	
	// Allocate memory. Dimensions are limits.second x (limits.second-1)
	typedef MSA_di_unipi_it::MSArbor::CNumber ged_type;
	MSA_di_unipi_it::MSArbor::CRow upper_bounds = new MSA_di_unipi_it::MSArbor::CNumber[ (limits.second) * (limits.second - 1) ];	
	if(stdout>0) std::cout<<"GED matrix memory allocated"<<std::endl;

	ged::GEDGraph::GraphID i_par, j_par;
	std::vector<ged::GEDGraph::GraphID> subset;
	std::vector<ged::GEDGraph::GraphID> population(limits.second-1,0);
	// Allocate memory

	// Do not consider empty graph (last graph)
	for(i_par = limits.first; i_par<limits.second-1; i_par++){
		population.at(i_par) = i_par;
	}
	std::size_t graph_sample_size, aux_graph_sample_size;
	bool complete = true;
	double factor=1;
	if(args.count("graph_sample_size")>0){
		if(args.count("graph_sample_type")>0 && args.at("graph_sample_type")=="%"){
			factor = static_cast<double>(limits.second)/100 ;
		}
		aux_graph_sample_size = std::stoi(args.at("graph_sample_size")) * factor;
		graph_sample_size = (aux_graph_sample_size < limits.second-2) ? aux_graph_sample_size: limits.second-2 ;
		if (graph_sample_size <limits.second-2) complete=false;
    }
    else{
    	graph_sample_size = limits.second-2; // Max number of sample
    }

    if (stdout >0) std::cout <<"graph_sample_size: " << graph_sample_size << std::endl;
    ged::ProgressBar progress( graph_sample_size * (limits.second-1) + limits.second-1);
	if (stdout >0) std::cout << "\rComputing GED: " << progress << std::flush;

	// The type used in the std::numeric_limits should be the same used in MSArbor.h for CNumber
	ged_type max_arc_cost = MSA_di_unipi_it::MSArbor::C_INF - 1;
	// This part needs to be calibrated so that there is no overflow
	std::size_t omega = std::numeric_limits<std::size_t>::max()*omega_prop_;
	short int last_iter=0, this_iter=0;
	
	std::size_t cost_extra_constant = 2*b_ni+3*b_ei+1;

	std::size_t to_add = cost_extra_constant;
	if(relaxed) to_add += b_ni;

	std::size_t V1=0, V2=0;


	auto start_gedlib = std::chrono::high_resolution_clock::now();
	runtime = start_gedlib - start_init;
	double initialization_time = runtime.count();
	bool add_following = false;
	if(args.count("path_structure")>0 && args.at("path_structure")=="true") add_following = true;


	bool run_ged = true;
	if(args.count("match_node_map")>0 && args.at("match_node_map")=="true") run_ged = false;

	std::vector<std::pair<std::size_t, std::string>> att_1;
	std::vector<std::pair<std::size_t, std::string>> att_2;

	for(i_par = limits.first; i_par<limits.second; i_par++){
		g1 = env_coded.get_graph(i_par, true, false, true);

		if(!run_ged){
			att_1.clear();
			for(std::size_t h = 0; h<g1.num_nodes; h++){
				att_1.emplace_back(std::make_pair(h, g1.node_labels.at(h).at(args.at("match_node_map_by"))));
			}			
			std::sort(att_1.begin(), att_1.end(), pair_sort);			
		}
		
		// get the k graphs to calculate the distance to 
		if(!complete && i_par != empty_id){
			subset.clear();
			subset = random_sample(population, graph_sample_size, i_par);
			// If we are working with the msts that have a "path" relationship, add the following graph
			if(add_following && i_par+1 != empty_id && std::find(subset.begin(), subset.end(), i_par+1)==subset.end()) subset.emplace_back(i_par+1);

		}
		else{
			subset = population; // Does not matter
		}
		for(j_par = limits.first; j_par<limits.second-1; j_par++){
			if (stdout >4) std::cout << "Computing " << i_par << "->" << j_par << std::endl;
			if(j_par!=empty_id){
				if(i_par==j_par){
					upper_bounds[i_par + (limits.second)*j_par] = max_arc_cost;
				}
				else{
					if(i_par != empty_id && !complete && std::find(subset.begin(), subset.end(), j_par) == subset.end()){
						upper_bounds[i_par + (limits.second)*j_par] = max_arc_cost;
					}
					else{	
						g2 = env_coded.get_graph(j_par, true, false, true);	
						if(run_ged==false){
							att_2.clear();
							for(std::size_t h = 0; h<g2.num_nodes; h++){
								att_2.emplace_back(std::make_pair(h, g2.node_labels.at(h).at(args.at("match_node_map_by"))));								
							}					
							std::sort(att_2.begin(), att_2.end(), pair_sort);							
							// Find trivial node_maps

							ged::NodeMap trivial_node_map = get_trivial_node_map<UserNodeID, UserNodeLabel, UserEdgeLabel>(att_1, att_2);
							env_coded.set_calculation_values(i_par, j_par, trivial_node_map, 0, 0, runtime);

						}
						else{				
							if(edit_modified){							
								// Need to compare the two graphs to define the edit cost constants.
								V1 = env_coded.get_num_nodes(i_par);
								V2 = env_coded.get_num_nodes(j_par);
								this_iter = (V1>=V2)? 1 : 2;
								if(this_iter != last_iter)  {	
									modified_costs(comp_costs,V1, V2, 
									c_nd, c_ni, c_ns, c_ed, c_ei, c_es, c_es_id, omega,b_ni, b_na, b_ei,b_ea);						
									env_coded.set_edit_costs(ged::Options::EditCosts::COMPRESSION, comp_costs);

								}
								last_iter = this_iter;
						
							}						
											
							env_coded.run_method(i_par,j_par);	
						}

						//already_calculated.at(i_par).at(j_par) = true;	
						upper_bounds[i_par + (limits.second)*j_par] = static_cast<ged_type>(compute_induced_compression_cost(env_coded.get_node_map(i_par,j_par), g1, g2, b_ni, b_na, b_ei,b_ea)+ to_add);						
						progress.increment();
						if (stdout >0) std::cout << "\rComputing GED: " << progress << std::flush;
					}
					
				}				
			}			
		}
	}
	if (stdout >0) std::cout<<std::endl;

	if(args.count("write_ged_matrix")>0 && args.at("write_ged_matrix")=="true"){
		aux_string = output_other + "/"+ file_preffix +"_GEDmatrix_k_" + std::to_string(graph_sample_size) + ".csv";
		write_matrix(aux_string, upper_bounds, limits.second);		
	}

	auto start_arb = std::chrono::high_resolution_clock::now();
	runtime = start_arb - start_gedlib;
	double gedlib_time = runtime.count();

	// 4. Spanning arborescence of minimum weight
	std::vector<std::size_t> arborescence;
	std::size_t cost_arb = 0;
	std::size_t root = empty_id;
	spanning_arborescence_of_minimum_weight(arborescence, cost_arb, upper_bounds,root, max_arc_cost);

	if(args.count("write_arb")>0 && args.at("write_arb")=="true"){
		aux_string = output_other +"/"+ file_preffix +"_arb_k_" + std::to_string(graph_sample_size) + ".csv";
		write_to_file<std::size_t>(aux_string, arborescence);
	}
	std::vector<std::size_t> arborescence_ref(arborescence);

	auto start_ref = std::chrono::high_resolution_clock::now();
	runtime = start_ref - start_arb;
	double arb_time = runtime.count();

	// 5. Refinement
	std::size_t refinement_size = 0;
	if(args.count("refinement_size")>0){
		refinement_size = std::stoi(args.at("refinement_size"));
	}	
	if(stdout>0) std::cout<<"Refinement size: "<<refinement_size<<std::endl;


	std::size_t cost_arb_ref = 0;
	double ref_time=0;
	double ref_arb_time=0;

	double lb=0, ub=0;
	ged::Seconds last_runtime = runtime;
	ged::NodeMap last_node_map = ged::NodeMap(1,1);
	ged_type aux_value=0;

	if(refinement_size>0 && run_ged){

		// 5.1 Set method
		std::string ged_method_refinement_options = "";
		if(args.count("ged_method_refinement_options")>0){
			ged_method_refinement_options = args.at("ged_method_refinement_options");
		}
		else{
			// Default number of threads
			ged_method_refinement_options = "16";	
		}

		// Currently working with BRANCH_UNIFORM. All other methods have not been thoroughly tested
		if(args.count("ged_method_refinement")>0){
			if(stdout) std::cout<<"GED method (Refinement): "<<args.at("ged_method_refinement")<<std::endl;

			if(args.at("ged_method_refinement") == "branch_uniform"){
				env_coded.set_method(ged::Options::GEDMethod::BRANCH_UNIFORM, "--threads " + ged_method_refinement_options);
			}
			else{
				if(args.at("ged_method_refinement") == "branch_fast"){
					env_coded.set_method(ged::Options::GEDMethod::BRANCH_FAST, "--threads " + ged_method_refinement_options);
				}
				else{
					if(args.at("ged_method_refinement") == "ipfp"){
						env_coded.set_method(ged::Options::GEDMethod::IPFP, "--threads " + ged_method_refinement_options);
					}
					else{
						std::cout<<"No valid ged_method in args. Setting it to BRANCH_UNIFORM"<<std::endl;
						env_coded.set_method(ged::Options::GEDMethod::BRANCH_UNIFORM, "--threads " + ged_method_refinement_options);
					}
				}
			}
		}
		else{
			std::cout<<"No field ged_method in args. Setting it to BRANCH_UNIFORM"<<std::endl;
			env_coded.set_method(ged::Options::GEDMethod::BRANCH_UNIFORM, "--threads " + ged_method_refinement_options);
			
		}

		// 5.2 GED calculations

		std::size_t step, node;
		this_iter = 0;
		last_iter =0;
		
		to_add = cost_extra_constant;
		if(relaxed) to_add += b_ni;


		std::size_t cont = 0;

		for(std::size_t n =0; n<arborescence.size(); n++){
			step=0;
			node = n;
			while(step<refinement_size && node!=root){
				cont++;
				step++;
				node = arborescence.at(node);
			}
		}

		ged::ProgressBar progress( cont );
		if (stdout >0) std::cout << "\rComputing GED (refinement): " << progress << std::flush;
		for(std::size_t n =0; n<arborescence.size(); n++){
			step=0;
			node = n;
			j_par = n; // Fix the bottom node
			while(step<refinement_size && node!=root){
				step++;
				i_par = arborescence.at(node);
				if(i_par != j_par && j_par != empty_id){
					if(edit_modified){
						// Need to compare the two graphs to define the edit cost constants.
						V1 = env_coded.get_num_nodes(i_par);
						V2 = env_coded.get_num_nodes(j_par);
						this_iter = (V1>=V2)? 1 : 2;
						if(this_iter != last_iter)  {					
							modified_costs(comp_costs,V1, V2, 
							c_nd, c_ni, c_ns, c_ed, c_ei, c_es, c_es_id, omega,b_ni, b_na, b_ei,b_ea);				
							env_coded.set_edit_costs(ged::Options::EditCosts::COMPRESSION, comp_costs);
						}
						last_iter = this_iter;
						
					}
					
						// Save last values to replace them in the environment if they are better
						lb = env_coded.get_lower_bound(i_par, j_par);
						ub = env_coded.get_upper_bound(i_par, j_par);
						last_runtime = ged::Seconds(env_coded.get_runtime(i_par, j_par));
						last_node_map = env_coded.get_node_map(i_par, j_par);	
					
			
					env_coded.run_method(i_par, j_par);	
					
					g1 = env_coded.get_graph(i_par, true, false, true);
					g2 = env_coded.get_graph(j_par, true, false, true);	
					
					aux_value = static_cast<ged_type>(compute_induced_compression_cost(env_coded.get_node_map(i_par,j_par),g1,g2, b_ni, b_na, b_ei,b_ea)+ to_add);
					if(aux_value < upper_bounds[i_par + (limits.second)*j_par] ){
						upper_bounds[i_par + (limits.second)*j_par] = aux_value;
					}
					else{
						
							// Put initial values and forget about the new calculation
							env_coded.set_calculation_values(i_par, j_par, last_node_map, lb, ub, last_runtime);	
						
					}
				}
				node = arborescence.at(node);
				progress.increment();
				if (stdout >0) std::cout << "\rComputing GED (refinement): " << progress << std::flush;
			}
		}

		if (stdout >0) std::cout << std::endl;

		if(args.count("write_ged_matrix")>0 && args.at("write_ged_matrix")=="true"){
			aux_string = output_other +"/"+ file_preffix +"_GEDmatrix_refined_k_" + std::to_string(graph_sample_size) + ".csv";
			write_matrix(aux_string, upper_bounds, limits.second);
		}

		auto end_ref = std::chrono::high_resolution_clock::now();
		runtime = end_ref - start_ref;
		ref_time = runtime.count();

		// 6. Spanning arborescence of minimum weight
		spanning_arborescence_of_minimum_weight(arborescence_ref, cost_arb_ref, upper_bounds,root, max_arc_cost);

		if(args.count("write_arb")>0 && args.at("write_arb")=="true"){
			aux_string = output_other + "/"+ file_preffix +"_arb_refined_k_" + std::to_string(graph_sample_size) + ".csv";
			write_to_file<std::size_t>(aux_string, arborescence_ref);
		}

		auto end_ref_arb = std::chrono::high_resolution_clock::now();
		runtime = end_ref_arb - end_ref;
		ref_arb_time = runtime.count();
	}
	
	delete[] upper_bounds;	

	std::size_t info_file_size=0;
	auto end_ref_arb = std::chrono::high_resolution_clock::now();



	if(args.count("test_mode")>0 && args.at("test_mode")=="true"){
		// skip the writing of the files
	}
	else{
		// 7. Encode
		// 7.1 info_file
		
		create_info_file<UserNodeID,UserNodeLabel,UserEdgeLabel>(binary, output_encoded, file_preffix,
				env_coded, alphabets, attr_sizes,
				attr_types,
				b_ni, b_ei,
				arborescence_ref, root,
				fast_node_translate, fast_edge_translate,
			 	stdout, separator_1, separator_2);

		info_file_size = get_file_size(output_encoded + "/_000_" + file_preffix + ".met");

		// 7.2 encode graphs using arborescence
		encode_arborescence<UserNodeID,UserNodeLabel,UserEdgeLabel>(binary, relaxed, output_encoded, file_preffix,
				env_coded, attr_sizes, b_ni, b_ei, arborescence_ref, root,
				stdout, separator_1);

	}
	
	auto end_all = std::chrono::high_resolution_clock::now();
	runtime = end_all - end_ref_arb;
	double final_time = runtime.count();



	// 8. If results are asked, fill them
	if(args.count("write_results")>0 && args.at("write_results")=="true"){
		headers = std::vector<std::string>();
		values = std::vector<std::string>();

		headers.emplace_back("file_preffix");
		values.emplace_back(file_preffix);

		headers.emplace_back("edit_costs");
		if(edit_modified){
			values.emplace_back("modified");
		}
		else{
			values.emplace_back("traditional");
		}

		headers.emplace_back("relaxed_compression");
		if(relaxed){
			values.emplace_back("true");
		}
		else{
			values.emplace_back("false");
		}

		headers.emplace_back("binary_mode");
		if(binary){
			values.emplace_back("true");
		}
		else{
			values.emplace_back("false");
		}

		headers.emplace_back("ged_method");
		values.emplace_back(args.at("ged_method"));

		headers.emplace_back("b_ni");
		values.emplace_back(std::to_string(b_ni));

		headers.emplace_back("b_na");
		values.emplace_back(std::to_string(b_na));

		headers.emplace_back("b_ei");
		values.emplace_back(std::to_string(b_ei));

		headers.emplace_back("b_ea");
		values.emplace_back(std::to_string(b_ea));

		headers.emplace_back("graph_sample_size");
		values.emplace_back(std::to_string(graph_sample_size));

		headers.emplace_back("info_file_size");
		values.emplace_back(std::to_string(info_file_size));


		headers.emplace_back("cost_arborescence");
		values.emplace_back(std::to_string(cost_arb));


		headers.emplace_back("k_complexity_a");
		values.emplace_back(std::to_string(base_compr_cost_triangular_matrix<UserNodeID, UserNodeLabel, UserEdgeLabel>(env_coded, b_ni, b_na, b_ei, b_ea, 0, true)));

		headers.emplace_back("k_complexity_b");
		values.emplace_back(std::to_string(base_compr_cost_edge_pairs<UserNodeID, UserNodeLabel, UserEdgeLabel>(env_coded, b_ni, b_na, b_ei, b_ea, 0, true)));

		headers.emplace_back("load_time");
		values.emplace_back(std::to_string(load_time));

		headers.emplace_back("initialization_time");
		values.emplace_back(std::to_string(initialization_time));

		headers.emplace_back("gedlib_runtime_initial");
		values.emplace_back(std::to_string(gedlib_time));

		headers.emplace_back("spanning_arb_runtime");
		values.emplace_back(std::to_string(arb_time));

		headers.emplace_back("refinement_size");
		values.emplace_back(std::to_string(refinement_size));

		headers.emplace_back("cost_arborescence_refined");
		values.emplace_back(std::to_string(cost_arb_ref));

		headers.emplace_back("gedlib_runtime_refinement");
		values.emplace_back(std::to_string(ref_time));

		headers.emplace_back("refine_arb_runtime");
		values.emplace_back(std::to_string(ref_arb_time));

		headers.emplace_back("final_time");
		values.emplace_back(std::to_string(final_time));

		
		std::map<std::size_t, std::vector<std::size_t>> depth_degrees;
		get_arborescence_info(headers, values, depth_degrees, arborescence_ref, root);	

	}

}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void 
GED_ABC::
decode_single_graph(bool binary, bool relaxed, std::string path, std::string graph_name, std::string graph_class,
	ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env,
	std::size_t parent_num, std::size_t child,
	std::size_t root,
	std::map<std::size_t, std::size_t> &pos_to_id,
	std::size_t b_ni, std::size_t b_ei,
	std::map<std::string, std::vector<std::string>> &ordered_attributes,
	std::map<std::string, std::map<std::string, std::size_t>> &attr_sizes,
	int stdout, char separator_1){

	unsigned short int type_read=0;

	ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> g1, g2;

	ged::NodeMap node_map_aux = ged::NodeMap(1,1);
	
	std::vector<std::size_t> v_d,v_s,v_rest, v_i, v_is, e_ed, e_is, e_s;
	std::vector<std::pair<std::size_t,std::size_t>> e_ed_v, e_is_v, e_s_v;
	
	std::map<std::pair<std::size_t, std::size_t>, bool> marker;
	std::map<std::pair<std::size_t, std::size_t>, UserEdgeLabel> marker_labels;
	std::size_t num_nodes, num_subs, size_read, graph_id, from, to, num_edges, aux_read, first, second, cont, limit;
	std::size_t num_subs_is, num_ins; // relaxed
	std::vector<UserNodeID> aux_node_ids;
	std::vector<UserNodeLabel> aux_node_labels, new_node_labels;
	
	std::vector<UserEdgeLabel> aux_edge_labels;
	std::vector<std::pair<std::size_t, std::size_t>> aux_edges;
	std::map<std::size_t, std::size_t> map_indices;
	std::map<std::string, std::string> label;

	typename std::list<std::pair<std::pair<std::size_t, std::size_t>,UserEdgeLabel>>::iterator iter;
	std::pair<std::pair<std::size_t, std::size_t>,UserEdgeLabel> edge;

	std::ifstream graph_in;
	std::string line;

	if(binary){
		graph_in.open(path.c_str(), ios::binary);	
	}
	else{
		graph_in.open(path.c_str());
	}
	

	if(graph_in.is_open()){
		if(parent_num == root){
			ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> empty_env;
			ged::GEDGraph::GraphID empty_id = empty_env.add_graph("empty","");
			g1 = empty_env.get_graph(empty_id, false, false, true);
			g2 = empty_env.get_graph(empty_id, false, false, true);
		}
		else{					
			g1 = env.get_graph(pos_to_id.at(parent_num), false, false, true);
			g2 = env.get_graph(pos_to_id.at(parent_num), false, false, true);
		}
		if(stdout>3){
			std::cout<<" ____________ BEGIN SOURCE GRAPH ______________"<<std::endl;
			describe_graph(g1);
			std::cout<<" ____________ END SOURCE GRAPH ______________"<<std::endl;					
		}
		
		// V
		// Line 1: # of vertices
		if(binary){
			aux_read = get_size_t_from_bytes(graph_in, b_ni);	
		}
		else{
			std::getline(graph_in,line,separator_1);
			aux_read = std::stoi(line);
		}
		if(stdout>3) std::cout<<aux_read<<std::endl;
		num_nodes = aux_read;
		
		//Line 2: # subs	
		if(binary){
			aux_read = get_size_t_from_bytes(graph_in, b_ni);	
		}
		else{
			std::getline(graph_in,line,separator_1);
			aux_read = std::stoi(line);
		}
		if(stdout>3) std::cout<<aux_read<<std::endl;			
		num_subs = aux_read;


		if(relaxed){
			//Line 3: # is subs (only if relaxed)	
			if(binary){
				aux_read = get_size_t_from_bytes(graph_in, b_ni);	
			}
			else{
				std::getline(graph_in,line,separator_1);
				aux_read = std::stoi(line);
			}
			if(stdout>3) std::cout<<aux_read<<std::endl;			
			num_subs_is = aux_read;	
		}
		else{
			num_subs_is = 0;
		}
		

		// Determine the situation							
		if(relaxed){
			if(g1.num_nodes - num_subs - num_subs_is < num_subs_is){
				type_read = 0;	// read v_d
			}
			else{
				type_read = 1;	// read v_is
			}		
		}
		else{
			if(g1.num_nodes >= num_nodes){
				if(g1.num_nodes - num_nodes < num_nodes - num_subs){
					type_read = 0;	// read v_d
				}
				else{
					type_read = 1;	// read v_is
				}
			}
			else{
				type_read = 2; // insertions
			}
		}
		
		
		if(stdout>2) std::cout<<"Type read: "<<type_read<<std::endl;

		cont=0;
		// Read either v_d or v_is or insertions
		switch(type_read){
			case 0:
				// Read v_d

			 	limit = (relaxed) ? g1.num_nodes - num_subs - num_subs_is : g1.num_nodes - num_nodes;
				for (std::size_t i = 0; i < limit; i++)
				{				
					if(binary){
						aux_read = get_size_t_from_bytes(graph_in, b_ni);	
					}
					else{
						std::getline(graph_in,line,separator_1);
						aux_read = std::stoi(line);
					}
					if(stdout>3) std::cout<<aux_read<<std::endl;
					v_d.emplace_back(aux_read);
					v_rest.emplace_back(aux_read);
				}

				break;
			case 1:
				// Read v_is
				limit = (relaxed) ? num_subs_is : num_nodes - num_subs;
				for (std::size_t i = 0; i < limit; i++)
				{
					if(binary){
						aux_read = get_size_t_from_bytes(graph_in, b_ni);	
					}
					else{
						std::getline(graph_in,line,separator_1);
						aux_read = std::stoi(line);
					}
					v_is.emplace_back(aux_read);
					v_rest.emplace_back(aux_read);
				} 

				break;
			case 2:
				// Read Node insertions					
				if(stdout>2) std::cout<<"Insertions "<<num_nodes - g1.num_nodes<<std::endl;
				if(stdout>2) std::cout<<g1.num_nodes<<" to "<< num_nodes<<std::endl;

				cont = g1.num_nodes;
				for (std::size_t i = 0; i < num_nodes - g1.num_nodes; i++)
				{							
											
					v_i.emplace_back(cont);
					cont++;
					
					label.clear();
					for(const auto attr: ordered_attributes.at("node_attr")){
						if(binary){
							aux_read = get_size_t_from_bytes(graph_in, attr_sizes.at("node_attr").at(attr));	
						}
						else{
							std::getline(graph_in,line,separator_1);
							aux_read = std::stoi(line);
						}
						label.emplace(std::make_pair(attr, std::to_string(aux_read)));
					}							
					new_node_labels.emplace_back(label);							
				}

				break;
		}

		// Read Node substitutions
		if(stdout>2) std::cout<<"Substitutions "<<num_subs<<std::endl;
		
		for (std::size_t i = 0; i < num_subs; i++){

			if(binary){
				aux_read = get_size_t_from_bytes(graph_in, b_ni);	
			}
			else{
				std::getline(graph_in,line,separator_1);
				aux_read = std::stoi(line);
			}; // index in g1 tilda
			if(stdout>3) std::cout<<aux_read<<std::endl;
			
			v_s.emplace_back(aux_read);
			v_rest.emplace_back(aux_read);
			map_indices.emplace(aux_read, aux_node_labels.size());					

			label.clear();
			for(const auto attr: ordered_attributes.at("node_attr")){
				if(binary){
					aux_read = get_size_t_from_bytes(graph_in, attr_sizes.at("node_attr").at(attr));	
				}
				else{
					std::getline(graph_in,line,separator_1);
					aux_read = std::stoi(line);
				}
				
				label.emplace(std::make_pair(attr, std::to_string(aux_read)));
			}
			aux_node_labels.emplace_back(label);
		}

		if(relaxed){
			// Node insertions
			num_ins = num_nodes - num_subs_is - num_subs;
			if(stdout>2) std::cout<<"Insertions "<<num_ins<<std::endl;
			if(stdout>2) std::cout<<g1.num_nodes<<" to "<< num_nodes<<std::endl;
			cont = num_subs_is + num_subs;
			for (std::size_t i = 0; i < num_ins; i++)
			{
				v_i.emplace_back(cont);
				cont++;	
				label.clear();
				for(const auto attr: ordered_attributes.at("node_attr")){
					if(binary){
						aux_read = get_size_t_from_bytes(graph_in, attr_sizes.at("node_attr").at(attr));	
					}
					else{
						std::getline(graph_in,line,separator_1);
						aux_read = std::stoi(line);
					}
					
					label.emplace(std::make_pair(attr, std::to_string(aux_read)));
				}							
				new_node_labels.emplace_back(label);
				
			}
		}

		// End V

		if(stdout>2) std::cout<< "Deduce v_d os v_is"<<std::endl;
		// Now deduce v_d to create the auxiliary node map
		if(type_read!=2){ // V > V'
			if(type_read==1){
				for(std::size_t n = 0; n<g1.num_nodes; n++){
					if(std::find(v_rest.begin(), v_rest.end(), n) == v_rest.end()){
						v_d.emplace_back(n);
						if(stdout>3) std::cout<< n << " to v_d"<<std::endl;
					}
				}
			}
		}

		if(type_read!=1){
			for(std::size_t n = 0; n<g1.num_nodes; n++){
				if(std::find(v_rest.begin(), v_rest.end(), n) == v_rest.end()){
					v_is.emplace_back(n);
					if(stdout>3) std::cout<< n << " to v_is"<<std::endl;
				}
			}
		}

		if(stdout>2) std::cout<< "Get auxiliary node map"<<std::endl;
		// NodeMap between g1 tilda and g2 tilda
		node_map_aux = get_aux_node_map(g1.num_nodes, num_nodes, v_d, v_i);
		if(stdout>2) std::cout<< node_map_aux<<std::endl;


		if(stdout>2) std::cout<< "Index correction"<<std::endl;
		// Correct the indices (g1 tilda to g2 tilda) when needed
		g2.node_labels.clear();
		g2.original_node_ids.clear();
		g2.num_nodes=num_nodes;

		label.clear();
		g2.node_labels = std::vector<UserNodeLabel>(num_nodes, label);
		g2.original_node_ids = std::vector<UserNodeID>(num_nodes, "");
		if(stdout>3) std::cout<< "g2.node_labels.size(): " << g2.node_labels.size() << std::endl;
		if(stdout>3) std::cout<< "aux_node_labels.size(): " << aux_node_labels.size() << std::endl;
		if(stdout>3) std::cout<< "new_node_labels.size(): " << new_node_labels.size() << std::endl;

		if(stdout>2) std::cout<< " ==== V_S ==="<<std::endl;
		for(auto const &v : v_s){
			// Send to g2 tilda
			if(stdout>3) std::cout<< v << " -> "<< node_map_aux.image(v) << " - > " << map_indices.at(v) <<std::endl;
			g2.node_labels.at(node_map_aux.image(v)) = aux_node_labels.at(map_indices.at(v));
		}

		if(stdout>2) std::cout<< " ==== V_IS ==="<<std::endl;
		for(auto const &v : v_is){
			// Send to g2 tilda
			if(stdout>3) std::cout<< v << " -> "<< node_map_aux.image(v) << " - > " << v <<std::endl;
			g2.node_labels.at(node_map_aux.image(v)) = g1.node_labels.at(v);
		}				

		if(stdout>2) std::cout<< " ==== V_I ==="<<std::endl;
		cont=0;
		for(auto const &v : v_i){
			// Already g2 tilda
			if(stdout>3) std::cout<< v <<std::endl;
			g2.node_labels.at(v) = new_node_labels.at(cont);
			cont++;
		}

		// Fix problem with node ids not corresponding to position in the node list
		for(std::size_t n=0; n< num_nodes; n++){
			g2.original_node_ids.at(n) = std::to_string(n);
		}
		
		// E
		if(stdout>2) std::cout<<"---------edges-------------------"<<std::endl;

		// init. Get current edges and labels (Need to be improved)
		if(stdout>2) std::cout<<"-------- Existing edges"<<std::endl;				
		for(iter = g1.edge_list.begin(); iter != g1.edge_list.end(); iter++){
			edge = (*iter);
			
			if(stdout>3) std::cout<<"Edge (in g1 tilde): "<<edge.first.first << " - "<< edge.first.second<<std::endl;
			label.clear();
			for(auto const & ll : edge.second){
				if(stdout>3) std::cout<<"\t"<<ll.first<< " = "<<ll.second<<std::endl;
				label.emplace(ll.first, ll.second);
			}
		
			if(stdout>3) std::cout<<"Edge (in g2 tilde): "<<node_map_aux.image(edge.first.first) << " - "<< node_map_aux.image(edge.first.second)<<std::endl;
			// marker in g1 tilde
			from = edge.first.first;
			to = edge.first.second;
			marker.emplace(std::make_pair(std::make_pair(from, to), true));
			marker_labels.emplace(std::make_pair(std::make_pair(from, to), label));
		}
		if(stdout>2) std::cout<<"-------- end Existing edges"<<std::endl;

		if(stdout>3){
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

		// Line 1: # edges
		if(binary){
			aux_read = get_size_t_from_bytes(graph_in, b_ei);	
		}
		else{
			std::getline(graph_in,line,separator_1);
			aux_read = std::stoi(line);
		}
		if(stdout>3) std::cout<<aux_read<<std::endl;
		num_edges = aux_read;			

		// Line 2: # subs
		if(binary){
			aux_read = get_size_t_from_bytes(graph_in, b_ei);	
		}
		else{
			std::getline(graph_in,line,separator_1);
			aux_read = std::stoi(line);
		}
		if(stdout>3) std::cout<<aux_read<<std::endl;
		num_subs = aux_read;

		//Line 3: type
		if(binary){
			aux_read = get_size_t_from_bytes(graph_in, 1);	
		}
		else{
			std::getline(graph_in,line,separator_1);
			aux_read = std::stoi(line);
		}
		if(stdout>3) std::cout<<aux_read<<std::endl;
		type_read = aux_read;

		//Line 4: size of the set to read (either e_ed or e_is)
		if(binary){
			aux_read = get_size_t_from_bytes(graph_in, b_ei);	
		}
		else{
			std::getline(graph_in,line,separator_1);
			aux_read = std::stoi(line);
		}
		if(stdout>3) std::cout<<aux_read<<std::endl;
		size_read = aux_read;

		// Read either e_ed or e_is or insertions
		if(stdout>2) std::cout<<"Read e_ed or e_is"<<std::endl;
		for (std::size_t i = 0; i < size_read; i++){
			if(binary){
				aux_read = get_size_t_from_bytes(graph_in, b_ni);	
			}
			else{
				std::getline(graph_in,line,separator_1);
				aux_read = std::stoi(line);
			}					
			first = aux_read;

			if(binary){
				aux_read = get_size_t_from_bytes(graph_in, b_ni);	
			}
			else{
				std::getline(graph_in,line,separator_1);
				aux_read = std::stoi(line);
			}					
			second = aux_read;

			// Nodes already in g2 tilde					
			if(stdout>3) std::cout<<first<<" -> "<<second<<std::endl;
			from = node_map_aux.pre_image(first);
			to = node_map_aux.pre_image(second);
			
			if(type_read == 0){
				e_ed_v.emplace_back(std::make_pair(first,second));						
			}
			else{	
				e_is_v.emplace_back(std::make_pair(first,second));
				e_is.emplace_back(0);  // For the limits in the other loops

				// add information in the temp vectors						
				aux_edges.emplace_back(std::make_pair(first,second));

				if(marker_labels.count(std::make_pair(from,to))>0){
					aux_edge_labels.emplace_back(marker_labels.at(std::make_pair(from,to)));
				}
				else{
					if(marker_labels.count(std::make_pair(to,from))>0){
						aux_edge_labels.emplace_back(marker_labels.at(std::make_pair(to,from)));
					}
					else{						
						throw compression_exception( "Edge not found: " + std::to_string(from) + " -> " + std::to_string(to) );
					}
				}
			}
			
			// Must change for directed graphs
			if(marker.count(std::make_pair(from,to))>0){
				marker.at(std::make_pair(from,to))=false;
			}
			else{
				if(marker.count(std::make_pair(to,from))>0){
					marker.at(std::make_pair(to, from))=false;
				}
				else{							
					throw compression_exception( "Edge not found: " + std::to_string(from) + " -> " + std::to_string(to) );
				}
			}
		}

		// Read Edge substitutions
		if(stdout>2) std::cout<<"Substitutions "<<num_subs<<std::endl;

		for (std::size_t i = 0; i < num_subs; i++){		
			if(binary){
				aux_read = get_size_t_from_bytes(graph_in, b_ni);	
			}
			else{
				std::getline(graph_in,line,separator_1);
				aux_read = std::stoi(line);
			}					
			first = aux_read;

			if(binary){
				aux_read = get_size_t_from_bytes(graph_in, b_ni);	
			}
			else{
				std::getline(graph_in,line,separator_1);
				aux_read = std::stoi(line);
			}					
			second = aux_read;

			from = first;
			to = second;
			
			aux_edges.emplace_back(std::make_pair(from, to));

			from = node_map_aux.pre_image(first);
			to = node_map_aux.pre_image(second);

			// Must change for directed graphs
			if(marker.count(std::make_pair(from,to))>0){
				marker.at(std::make_pair(from,to))=false;
			}
			else{
				if(marker.count(std::make_pair(to,from))>0){
					marker.at(std::make_pair(to, from))=false;
				}
				else{							
					throw compression_exception( "Edge not found: " + std::to_string(from) + " -> " + std::to_string(to) );
				}
			}
			
			label.clear();	
			for(const auto attr: ordered_attributes.at("edge_attr")){
				if(binary){
					aux_read = get_size_t_from_bytes(graph_in, attr_sizes.at("edge_attr").at(attr));	
				}
				else{
					std::getline(graph_in,line,separator_1);
					aux_read = std::stoi(line);
				}				
				label.emplace(std::make_pair(attr, std::to_string(aux_read)));
			}
			aux_edge_labels.emplace_back(label);

		}

		// Delete edges by node deletion
		if(stdout>2) std::cout<<"Deletions by node deletion "<<std::endl;
		cont=0;
		for(iter = g1.edge_list.begin(); iter != g1.edge_list.end(); iter++){
			edge = (*iter);
			for(auto const & v: v_d){
				if(edge.first.first == v || edge.first.second == v){
					if(stdout>3) std::cout<<"Edge deleted in g1 tilde: "<<edge.first.first << " - "<< edge.first.second<<std::endl;
					e_ed.emplace_back(cont);
					
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
							throw compression_exception( "Edge not found: " + std::to_string(from) + " -> " + std::to_string(to) );
						}
					}
					break;
				}
			}
			cont++;
		}


		// Deduce the set that was not read
		if(stdout>2) std::cout<<"Deduce other sets "<<std::endl;
		cont=0;
		if(type_read == 0){
			if(stdout>3) std::cout<<"deduce IS "<<marker.size()<<std::endl;
			for(const auto & m : marker){
				
				if(m.second){
					if(stdout>3) std::cout<< "Edge in IS in g1 tilde: "<<m.first.first << " - " << m.first.second << std::endl;
					
					from = node_map_aux.image(m.first.first);
					to = node_map_aux.image(m.first.second);

					if(stdout>3) std::cout<<" - after node_map_aux: " << from << " - " << to << std::endl;	

					marker.at(std::make_pair(m.first.first, m.first.second))=false;
					
					e_is_v.emplace_back(std::make_pair(from,to));
					aux_edges.emplace_back(std::make_pair(from,to));

					aux_edge_labels.emplace_back(marker_labels.at(m.first));
					e_is.emplace_back(cont); // For the limits in the other loops
				}
			}

			if(stdout>2) std::cout<<"END ______ deduce IS "<<std::endl;



		}
		else{
			if(stdout>2) std::cout<<"deduce ED "<<std::endl;				
			for(const auto & m : marker){
				
				if(m.second){
					if(stdout>3) std::cout<<"Edge in ED: "<<m.first.first << " - "<< m.first.second<<std::endl;
					from = node_map_aux.image(m.first.first);
					to = node_map_aux.image(m.first.second);
					marker.at(std::make_pair(m.first.first,m.first.second))=false;
					e_ed_v.emplace_back(std::make_pair(from,to));
					
				}
			}
			if(stdout>2) std::cout<<"END ______ deduce ED "<<std::endl;			
		}

		// Read Edge insertions
		if(num_edges > num_subs + e_is.size()){ 
			if(stdout>2) std::cout<<"Edge Insertions: "<< num_edges - num_subs - e_is.size() <<std::endl;
			
			for (std::size_t i = 0; i < num_edges - num_subs - e_is.size(); i++)
			{		
				if(binary){
					aux_read = get_size_t_from_bytes(graph_in, b_ni);	
				}
				else{
					std::getline(graph_in,line,separator_1);
					aux_read = std::stoi(line);
				}						
				first = aux_read;

				if(binary){
					aux_read = get_size_t_from_bytes(graph_in, b_ni);	
				}
				else{
					std::getline(graph_in,line,separator_1);
					aux_read = std::stoi(line);
				}						
				second = aux_read;

				if(stdout>3) std::cout<<first<<" -> "<<second<<std::endl;						
				aux_edges.emplace_back(std::make_pair(first,second));
				
				label.clear();
				for(const auto attr: ordered_attributes.at("edge_attr")){
					if(binary){
						aux_read = get_size_t_from_bytes(graph_in, attr_sizes.at("edge_attr").at(attr));	
					}
					else{
						std::getline(graph_in,line,separator_1);
						aux_read = std::stoi(line);
					}	
					label.emplace(std::make_pair(attr, std::to_string(aux_read)));
				}
				aux_edge_labels.emplace_back(label);
			}

		}

		g2.edge_list.clear();
		if (stdout>2) std::cout<<"Build edges: "<<aux_edges.size()<<", "<<aux_edge_labels.size()<<std::endl;
		for(std::size_t i =0; i<aux_edges.size(); i++){
			g2.edge_list.push_front(std::make_pair(aux_edges.at(i), aux_edge_labels.at(i)));
		}

		g2.num_edges = g2.edge_list.size();
		
		// End E
		if(stdout>3){
			std::cout<<" ____________ BEGIN RESULTING GRAPH ______________"<<std::endl;
			describe_graph(g2);
			std::cout<<" ____________ END RESULTING GRAPH ______________"<<std::endl;					
		}

		graph_id = env.load_exchange_graph(g2, ged::undefined(), ged::Options::ExchangeGraphType::EDGE_LIST, 
			graph_name, graph_class);
		if (stdout>2) std::cout<<"Done: "<<env.num_graphs()<<std::endl;
		pos_to_id.emplace(child, graph_id);
		graph_in.close();
		
	}
	else{				
		throw compression_exception( "Unable to open graph file " + path );
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void 
GED_ABC::
decode_single_graph_relaxed(bool binary, std::string path, std::string graph_name, std::string graph_class,
	ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> &env,
	std::size_t parent_num, std::size_t child,
	std::size_t root,
	std::map<std::size_t, std::size_t> &pos_to_id,
	std::size_t b_ni, std::size_t b_ei,
	std::map<std::string, std::vector<std::string>> &ordered_attributes,
	std::map<std::string, std::map<std::string, std::size_t>> &attr_sizes,
	int stdout){

	unsigned short int type_read=0;

	ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> g1, g2;

	ged::NodeMap node_map_aux = ged::NodeMap(1,1);
	
	std::vector<std::size_t> v_d,v_s,v_rest, v_i, v_is, e_ed, e_is, e_s;
	std::vector<std::pair<std::size_t,std::size_t>> e_ed_v, e_is_v, e_s_v;
	
	std::map<std::pair<std::size_t, std::size_t>, bool> marker;
	std::map<std::pair<std::size_t, std::size_t>, UserEdgeLabel> marker_labels;
	std::size_t num_nodes, num_subs, size_read, graph_id, from, to, num_edges, aux_read, first, second, cont;
	std::size_t num_subs_is, num_ins; // relaxed
	std::vector<UserNodeID> aux_node_ids;
	std::vector<UserNodeLabel> aux_node_labels, new_node_labels;
	
	std::vector<UserEdgeLabel> aux_edge_labels;
	std::vector<std::pair<std::size_t, std::size_t>> aux_edges;
	std::map<std::size_t, std::size_t> map_indices;
	std::map<std::string, std::string> label;

	typename std::list<std::pair<std::pair<std::size_t, std::size_t>,UserEdgeLabel>>::iterator iter;
	std::pair<std::pair<std::size_t, std::size_t>,UserEdgeLabel> edge;

	std::ifstream graph_in;
	graph_in.open(path.c_str());
	if(graph_in.is_open()){
		if(parent_num == root){
			ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> empty_env;
			ged::GEDGraph::GraphID empty_id = empty_env.add_graph("empty","");
			g1 = empty_env.get_graph(empty_id, false, false, true);
			g2 = empty_env.get_graph(empty_id, false, false, true);
		}
		else{					
			g1 = env.get_graph(pos_to_id.at(parent_num), false, false, true);
			g2 = env.get_graph(pos_to_id.at(parent_num), false, false, true);
		}
		if(stdout>3){
			std::cout<<" ____________ BEGIN SOURCE GRAPH ______________"<<std::endl;
			describe_graph(g1);
			std::cout<<" ____________ END SOURCE GRAPH ______________"<<std::endl;					
		}
		
		// V
		// Line 1: # of vertices
		aux_read = get_size_t_from_bytes(graph_in, b_ni);
		if(stdout>3) std::cout<<aux_read<<std::endl;
		num_nodes = aux_read;
		
		//Line 2: # subs	
		aux_read = get_size_t_from_bytes(graph_in, b_ni);
		if(stdout>3) std::cout<<aux_read<<std::endl;			
		num_subs = aux_read;

		//Line 3: # subs_is
		aux_read = get_size_t_from_bytes(graph_in, b_ni);
		if(stdout>3) std::cout<<aux_read<<std::endl;
		num_subs_is = aux_read;

		// Determine the situation							
		if(g1.num_nodes - num_subs - num_subs_is < num_subs_is){
			type_read = 0;	// read v_d
		}
		else{
			type_read = 1;	// read v_is
		}
	
		if(stdout>3) std::cout<<"Type read: "<<type_read <<std::endl;

		cont=0;
		// Read either v_d or v_is or insertions
		switch(type_read){
			case 0:
				// Read v_d
				for (std::size_t i = 0; i < g1.num_nodes - num_subs - num_subs_is; i++)
				{
					aux_read = get_size_t_from_bytes(graph_in, b_ni);
					if(stdout>3) std::cout<<aux_read<<std::endl;
					v_d.emplace_back(aux_read);
					v_rest.emplace_back(aux_read);
				}

				break;
			case 1:
				// Read v_is
				for (std::size_t i = 0; i < num_subs_is; i++)
				{
					aux_read = get_size_t_from_bytes(graph_in, b_ni);
					if(stdout>3) std::cout<<aux_read<<std::endl;
					v_is.emplace_back(aux_read);
					v_rest.emplace_back(aux_read);
				} 

				break;					
		}


		// Read Node substitutions
		if(stdout>2) std::cout<<"Substitutions "<<num_subs<<std::endl;
		
		for (std::size_t i = 0; i < num_subs; i++)
		{
			aux_read = get_size_t_from_bytes(graph_in, b_ni); // index in g1 tilda
			if(stdout>3) std::cout<<aux_read<<std::endl;
			
			v_s.emplace_back(aux_read);
			v_rest.emplace_back(aux_read);
			map_indices.emplace(aux_read, aux_node_labels.size());					

			label.clear();
			for(const auto attr: ordered_attributes.at("node_attr")){
				aux_read = get_size_t_from_bytes(graph_in, attr_sizes.at("node_attr").at(attr));
				label.emplace(std::make_pair(attr, std::to_string(aux_read)));
			}
			aux_node_labels.emplace_back(label);					

		}

		// Node insertions
		num_ins = num_nodes - num_subs_is - num_subs;
		if(stdout>2) std::cout<<"Insertions "<<num_ins<<std::endl;
		if(stdout>2) std::cout<<g1.num_nodes<<" to "<< num_nodes<<std::endl;
		cont = num_subs_is + num_subs;
		for (std::size_t i = 0; i < num_ins; i++)
		{
		
			v_i.emplace_back(cont);
			cont++;
			
			label.clear();
			for(const auto attr: ordered_attributes.at("node_attr")){
				aux_read = get_size_t_from_bytes(graph_in, attr_sizes.at("node_attr").at(attr));
				label.emplace(std::make_pair(attr, std::to_string(aux_read)));
			}							
			new_node_labels.emplace_back(label);
			
		}

		// End V

		if(stdout>2) std::cout<< "Deduce v_d os v_is"<<std::endl;
		// Now deduce v_d to create the auxiliary node map				
		if(type_read!=0){
			for(std::size_t n = 0; n<g1.num_nodes; n++){
				if(std::find(v_rest.begin(), v_rest.end(), n) == v_rest.end()){
					v_d.emplace_back(n);
					if(stdout>3) std::cout<< n << " to v_d"<<std::endl;
				}
			}
		}
		else{
			for(std::size_t n = 0; n<g1.num_nodes; n++){
				if(std::find(v_rest.begin(), v_rest.end(), n) == v_rest.end()){
					v_is.emplace_back(n);
					if(stdout>3) std::cout<< n << " to v_is"<<std::endl;
				}
			}
		}

		if(stdout>2) std::cout<< "Get auxiliary node map - Decoding"<<std::endl;
		// NodeMap between g1 tilda and g2 tilda
		node_map_aux = get_aux_node_map(g1.num_nodes, num_nodes, v_d, v_i);
		if(stdout>2) std::cout<< node_map_aux<<std::endl;


		if(stdout>2) std::cout<< "Index correction"<<std::endl;
		// Correct the indices (g1 tilda to g2 tilda) when needed
		g2.node_labels.clear();
		g2.original_node_ids.clear();
		g2.num_nodes=num_nodes;

		label.clear();
		g2.node_labels = std::vector<UserNodeLabel>(num_nodes, label);
		g2.original_node_ids = std::vector<UserNodeID>(num_nodes, "");


		if(stdout>2) std::cout<< " ==== V_S ==="<<std::endl;
		for(auto const &v : v_s){
			if(stdout>3) std::cout<< v << " -> "<< node_map_aux.image(v) <<std::endl;
			// Send to g2 tilda
			g2.node_labels.at(node_map_aux.image(v)) = aux_node_labels.at(map_indices.at(v));					
		}

		if(stdout>2) std::cout<< " ==== V_IS ==="<<std::endl;
		for(auto const &v : v_is){
			if(stdout>3) std::cout<< v << " -> "<< node_map_aux.image(v) <<std::endl;
			// Send to g2 tilda
			g2.node_labels.at(node_map_aux.image(v)) = g1.node_labels.at(v);
		}				

		if(stdout>2) std::cout<< " ==== V_I ==="<<std::endl;
		cont=0;
		for(auto const &v : v_i){
			if(stdout>3) std::cout<< v <<std::endl;
			// Already g2 tilda
			g2.node_labels.at(v) = new_node_labels.at(cont);
			cont++;					
		}

		// Fix problem with node ids not corresponding to position in the node list
		for(std::size_t n=0; n < num_nodes; n++){
			g2.original_node_ids.at(n) = std::to_string(n);
		}

		// -------------------------------------------------------------------
		// E
		if(stdout>2) std::cout<<"---------edges-------------------"<<std::endl;

		// init. Get current edges and labels
		if(stdout>3) std::cout<<"-------- Existing edges"<<std::endl;	

		for(iter = g1.edge_list.begin(); iter != g1.edge_list.end(); iter++){
			edge = (*iter);
			 
			if(stdout>3) std::cout<<"Edge (in g1 tilde): "<<edge.first.first << " - "<< edge.first.second<<std::endl;
			label.clear();
			for(auto const & ll : edge.second){
				if(stdout>3) std::cout<<"\t"<<ll.first<< " = "<<ll.second<<std::endl;
				label.emplace(ll.first, ll.second);
			}
		
			if(stdout>3) std::cout<<"Edge (in g2 tilde): "<<node_map_aux.image(edge.first.first) << " - "<< node_map_aux.image(edge.first.second)<<std::endl;
			// marker in g1 tilde
			from = edge.first.first;
			to = edge.first.second;
			marker.emplace(std::make_pair(std::make_pair(from, to), true));
			marker_labels.emplace(std::make_pair(std::make_pair(from, to), label));
		}
		if(stdout>3) std::cout<<"-------- end Existing edges"<<std::endl;

		if(stdout>3){
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

		// Line 1: # edges
		aux_read = get_size_t_from_bytes(graph_in, b_ei);
		if(stdout>3) std::cout<<aux_read<<std::endl;
		num_edges = aux_read;			

		// Line 2: # subs
		aux_read = get_size_t_from_bytes(graph_in, b_ei);
		if(stdout>3) std::cout<<aux_read<<std::endl;
		num_subs = aux_read;

		//Line 3: type
		aux_read = get_size_t_from_bytes(graph_in, 1);
		if(stdout>3) std::cout<<aux_read<<std::endl;
		type_read = aux_read;

		//Line 4: size of the set to read (either e_ed or e_is)
		aux_read = get_size_t_from_bytes(graph_in, b_ei);
		if(stdout>3) std::cout<<aux_read<<std::endl;
		size_read = aux_read;
		
		// Read either e_ed or e_is or insertions
		if(stdout>2) std::cout<<"Read e_ed or e_is"<<std::endl;
		for (std::size_t i = 0; i < size_read; i++)
		{
			aux_read = get_size_t_from_bytes(graph_in, b_ni);					
			first = aux_read;

			aux_read = get_size_t_from_bytes(graph_in, b_ni);					
			second = aux_read;

			from = node_map_aux.pre_image(first);
			to = node_map_aux.pre_image(second);

			if(type_read == 0){									
				e_ed_v.emplace_back(std::make_pair(first,second));						
			}
			else{
				e_is_v.emplace_back(std::make_pair(first,second));
				e_is.emplace_back(0);  // For the limits in the other loops
				
				// add information in the temp vectors						
				aux_edges.emplace_back(std::make_pair(first,second));
				if(marker_labels.count(std::make_pair(from,to))>0){
					aux_edge_labels.emplace_back(marker_labels.at(std::make_pair(from,to)));
				}
				else{
					if(marker_labels.count(std::make_pair(to,from))>0){						
						aux_edge_labels.emplace_back(marker_labels.at(std::make_pair(to,from)));
					}
					else{								
						throw compression_exception( "Edge not found: " + std::to_string(from) + " -> " + std::to_string(to) );
					}
				}
				

			}

			// Must change for directed graphs
			if(marker.count(std::make_pair(from,to))>0){
				marker.at(std::make_pair(from,to))=false;
			}
			else{
				if(marker.count(std::make_pair(to,from))>0){
					marker.at(std::make_pair(to, from))=false;
				}
				else{							
					throw compression_exception( "Edge not found: " + std::to_string(from) + " -> " + std::to_string(to) );
				}
			}

		}

		// Read Edge substitutions
		if(stdout>2) std::cout<<"Substitutions "<<num_subs<<std::endl;

		for (std::size_t i = 0; i < num_subs; i++)
		{		
			aux_read = get_size_t_from_bytes(graph_in, b_ni);					
			first = aux_read;

			aux_read = get_size_t_from_bytes(graph_in, b_ni);					
			second = aux_read;

			from = first;
			to = second;
			
			aux_edges.emplace_back(std::make_pair(from, to));

			from = node_map_aux.pre_image(first);
			to = node_map_aux.pre_image(second);

			// Must change for directed graphs
			if(marker.count(std::make_pair(from,to))>0){
				marker.at(std::make_pair(from,to))=false;
			}
			else{
				if(marker.count(std::make_pair(to,from))>0){				
					marker.at(std::make_pair(to, from))=false;
				}
				else{							
					throw compression_exception( "Edge not found: " + std::to_string(from) + " -> " + std::to_string(to) );
				}
			}
			

			label.clear();	
			for(const auto attr: ordered_attributes.at("edge_attr")){
				aux_read = get_size_t_from_bytes(graph_in, attr_sizes.at("edge_attr").at(attr));
				label.emplace(std::make_pair(attr, std::to_string(aux_read)));
			}
			aux_edge_labels.emplace_back(label);

		}

		// Delete edges by node deletion
		if(stdout>2) std::cout<<"Deletions by node deletion "<<std::endl;
		cont=0;
		for(iter = g1.edge_list.begin(); iter != g1.edge_list.end(); iter++){
			edge = (*iter);
			for(auto const & v: v_d){
				if(edge.first.first == v || edge.first.second == v){
					if(stdout>3) std::cout<<"Edge deleted in g1 tilde: "<<edge.first.first << " - "<< edge.first.second<<std::endl;
					e_ed.emplace_back(cont);							
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
							throw compression_exception( "Edge not found: " + std::to_string(from) + " -> " + std::to_string(to) );
						}
					}
					break;
				}
			}
			cont++;
		}


		// Deduce the set that was not read
		if(stdout>2) std::cout<<"Deduce other sets "<<std::endl;
		cont=0;
		if(type_read == 0){
			if(stdout>2) std::cout<<"deduce IS "<<marker.size()<<std::endl;
			// Now deduce IS					

			for(const auto & m : marker){
				
				if(m.second){
					if(stdout>3) std::cout<< "Edge in IS in g1 tilde: "<<m.first.first << " - " << m.first.second << std::endl;
					
					from = node_map_aux.image(m.first.first);
					to = node_map_aux.image(m.first.second);

					if(stdout>3) std::cout<<" - after node_map_aux: " << from << " - " << to << std::endl;	

					marker.at(std::make_pair(m.first.first, m.first.second))=false;
					
					e_is_v.emplace_back(std::make_pair(from,to));
					aux_edges.emplace_back(std::make_pair(from,to));
					
					aux_edge_labels.emplace_back(marker_labels.at(m.first));

					e_is.emplace_back(cont); // For the limits in the other loops
				}
			}

			if(stdout>2) std::cout<<"END ______ deduce IS "<<std::endl;



		}
		else{
			if(stdout>2) std::cout<<"deduce ED "<<std::endl;
			// Now deduce ED
			for(const auto & m : marker){
				
				if(m.second){
					if(stdout>3) std::cout<<"Edge in ED: "<<m.first.first << " - "<< m.first.second<<std::endl;
					from = node_map_aux.image(m.first.first);
					to = node_map_aux.image(m.first.second);
					marker.at(std::make_pair(m.first.first,m.first.second))=false;
					e_ed_v.emplace_back(std::make_pair(from,to));
					
				}
			}

			if(stdout>2) std::cout<<"END ______ deduce ED "<<std::endl;			
		}

		// Read Edge insertions
		if(num_edges > num_subs + e_is.size()){ 
			if(stdout>2) std::cout<<"Edge Insertions: "<< num_edges - num_subs - e_is.size() <<std::endl;
			
			for (std::size_t i = 0; i < num_edges - num_subs - e_is.size(); i++)
			{		
				aux_read = get_size_t_from_bytes(graph_in, b_ni);					
				first = aux_read;

				aux_read = get_size_t_from_bytes(graph_in, b_ni);					
				second = aux_read;
				
				aux_edges.emplace_back(std::make_pair(first,second));
				
				label.clear();
				for(const auto attr: ordered_attributes.at("edge_attr")){
					aux_read = get_size_t_from_bytes(graph_in, attr_sizes.at("edge_attr").at(attr));
					label.emplace(std::make_pair(attr, std::to_string(aux_read)));
				}
				aux_edge_labels.emplace_back(label);
			}
		}

		g2.edge_list.clear();
		if (stdout>2) std::cout<<"Build edges: "<<aux_edges.size()<<", "<<aux_edge_labels.size()<<std::endl;
		for(std::size_t i =0; i<aux_edges.size(); i++){
			g2.edge_list.push_front(std::make_pair(aux_edges.at(i), aux_edge_labels.at(i)));
		}

		g2.num_edges = g2.edge_list.size();
		
		// End E
		if(stdout>3){
			std::cout<<" ____________ BEGIN RESULTING GRAPH ______________"<<std::endl;
			describe_graph(g2);
			std::cout<<" ____________ END RESULTING GRAPH ______________"<<std::endl;					
		}

		graph_id = env.load_exchange_graph(g2, ged::undefined(), ged::Options::ExchangeGraphType::EDGE_LIST, 
			graph_name, graph_class);
		if (stdout>2) std::cout<<"Done: "<<env.num_graphs()<<std::endl;			
		pos_to_id.emplace(child, graph_id);
		graph_in.close();	
		
	}
	else{				
		throw compression_exception( "Unable to open graph file " + path );
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void 
GED_ABC::
decode_collection(bool binary, bool relaxed, ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> & env,
	std::string path,
	std::string file_preffix,
	std::map<std::string, std::string> &args,
	int stdout,
	char separator_1,
	char separator_2
	){	


	if(stdout>1) std::cout<<"decode_collection: Start DECODING"<<std::endl;

	// Read info_file to get the decoding structure
	std::map<std::string, std::map<std::string, std::map<std::string,std::string>>> decoded_attributes;
	decoded_attributes.emplace("graphs", std::map<std::string, std::map<std::string,std::string>> ());
	decoded_attributes.emplace("node_attr", std::map<std::string, std::map<std::string,std::string>> ());
	decoded_attributes.emplace("edge_attr", std::map<std::string, std::map<std::string,std::string>> ());
	
	std::ifstream in;
	std::string line, attr;
	std::map<std::string,std::string> value_map;
	
	std::size_t num_graphs, num_node_attr, num_edge_attr, aux_long;
	std::vector<std::size_t> arborescence;

	std::map<std::string, std::vector<std::string>> ordered_attributes;
	std::map<std::string, std::map<std::string, std::size_t>> attr_sizes;

	std::vector<std::string> node_attr_name, edge_attr_name;
	std::vector<std::size_t> node_attr_size, edge_attr_size;	
	
	std::size_t b_ni, b_ei, cont=0;

	std::vector<std::string> graph_names, graph_classes;

	std::string ext = ".cmp";
	graph_names = get_sorted_file_names(path, ext);

	std::string graph_file;

	std::string aux_string;
	int aux_int;
	float aux_float;
	double aux_double;
	char aux_char;
	char type;
	std::vector<std::string> aux_map_labels = {"node_attr", "edge_attr"};
	
	if(stdout>1) std::cout<<"decode_collection: info_file"<<std::endl;
	std::string info_file = path + "/_000_" + file_preffix + ".met";
	in.open(info_file.c_str());

	if(in.is_open()){

		num_graphs = graph_names.size();

		// graph classes
		if(binary) in.get(type); // get the type

		for(std::size_t i=0; i<num_graphs; i++){
			if(binary){
				switch(type){ 
		        case 'i': // int		        	
		            in.read(reinterpret_cast<char*>(&aux_int), sizeof(int));
		            aux_string = std::to_string(aux_int);
		            break;
		        case 'f': // float		        	
		            in.read(reinterpret_cast<char*>(&aux_float), sizeof(float));
		            aux_string = std::to_string(aux_float);
		            break;
		        case 'd': // double		        	
		            in.read(reinterpret_cast<char*>(&aux_double), sizeof(double));
		            aux_string = std::to_string(aux_double);
		            break;
		        case 'c': // char   
		            in.get(aux_char);
		            aux_string = aux_char;
		            break;	
		        case 's': // string
		        	aux_string = read_word_binary(in, separator_1);
		    	}	
			}
			else{
				std::getline(in,aux_string,separator_1);								
			}
			if(stdout>3) std::cout<<aux_string<<std::endl;
			graph_classes.emplace_back(aux_string);				
		}


		// arborescence
		arborescence.clear();
		if(binary){
			if(stdout>3) std::cout<<"arb begin"<<std::endl;
			// get size to read arborescence			
			// Fixed to 1
			aux_long = get_size_t_from_bytes(in, 1);
			if(stdout>3) std::cout<<aux_long<<" bytes to read arb"<<std::endl;
			
		}
		for(std::size_t i=0; i<num_graphs; i++){ // no empty graph in env
			if(binary){
				cont = get_size_t_from_bytes(in, aux_long);
			}
			else{
				std::getline(in,line,separator_1);
				cont = std::stoi(line);	
			}
			if(stdout>3) std::cout<<i<<" - "<<cont<<std::endl;
			arborescence.emplace_back(cont);
		}

		// # of attributes
		if(stdout>3) std::cout<<"attrs"<<std::endl;
		if(binary){
			num_node_attr = get_size_t_from_bytes(in, 1);
			if(stdout>3) std::cout<<num_node_attr<<std::endl;
			num_edge_attr = get_size_t_from_bytes(in, 1);
			if(stdout>3) std::cout<<num_edge_attr<<std::endl;
			b_ni = get_size_t_from_bytes(in, 1);
			if(stdout>3) std::cout<<b_ni<<std::endl;
			b_ei = get_size_t_from_bytes(in, 1);
			if(stdout>3) std::cout<<b_ei<<std::endl;
		}
		else{
			//nodes
			std::getline(in,line,separator_1);
			if(stdout>3) std::cout<<line<<std::endl;
			num_node_attr = std::stoi(line);
			//edges
			std::getline(in,line,separator_1);
			if(stdout>3) std::cout<<line<<std::endl;
			num_edge_attr = std::stoi(line);

			//nodes -> b_ni
			std::getline(in,line,separator_1);
			if(stdout>3) std::cout<<line<<std::endl;
			b_ni = std::stoi(line);
			//edges -> b_ei
			std::getline(in,line,separator_1);
			if(stdout>3) std::cout<<line<<std::endl;
			b_ei = std::stoi(line);
		}
			
		

		// size of attr
		if(stdout>3) std::cout<<"attr_sizes"<<std::endl;
		if(stdout>3) std::cout<<"nodes"<<std::endl;
		for(std::size_t i=0; i<num_node_attr; i++){
			if(binary){
				cont = get_size_t_from_bytes(in,1);
				node_attr_size.emplace_back(cont);
				if(stdout>3) std::cout<<cont<<std::endl;
			}
			else{
				std::getline(in,line,separator_1);
				if(stdout>3) std::cout<<line<<std::endl;
				node_attr_size.emplace_back(std::stoi(line));
			}
		}
		if(stdout>3) std::cout<<"edges"<<std::endl;
		for(std::size_t i=0; i<num_edge_attr; i++){
			if(binary){
				cont = get_size_t_from_bytes(in,1);
				edge_attr_size.emplace_back(cont);
				if(stdout>3) std::cout<<cont<<std::endl;
			}
			else{
				std::getline(in,line,separator_1);
				if(stdout>3) std::cout<<line<<std::endl;
				edge_attr_size.emplace_back(std::stoi(line));
			}
		}
		

		// Read alphabets
		if(binary){
			//nodes
			if(stdout>3) std::cout<<"nodes"<<std::endl;
			value_map.clear();
			for(std::size_t i=0; i<num_node_attr; i++){				
				//read the name
				line = read_word_binary(in, separator_1);
				if(stdout>3) std::cout<<line<<std::endl;
				node_attr_name.emplace_back(line);
				// get the type
				in.get(type);
				// NO more // read until separator_2 is found
				// Get # of values to read (coded in 4 bytes for now)
				cont = get_size_t_from_bytes(in, 4);				
				for(std::size_t j = 0; j< cont; j++){
					switch(type){ 
			        case 'i': // int		        	
			            in.read(reinterpret_cast<char*>(&aux_int), sizeof(int));
			            aux_string = std::to_string(aux_int);
			            break;
			        case 'f': // float		        	
			            in.read(reinterpret_cast<char*>(&aux_float), sizeof(float));
			            aux_string = std::to_string(aux_float);
			            break;
			        case 'd': // double		        	
			            in.read(reinterpret_cast<char*>(&aux_double), sizeof(double));
			            aux_string = std::to_string(aux_double);
			            break;
			        case 'c': // char   
			            in.get(aux_char);
			            aux_string = aux_char;
			            break;	
			        case 's': // string
			        	aux_string = read_word_binary(in, separator_1);
			    	}
			    	if(stdout>4) std::cout<<aux_string<<std::endl;
					value_map.emplace(std::to_string(j), aux_string);
				}

				value_map.emplace(std::to_string(cont), "dummy");
				if(stdout>3) std::cout<<"insert "<<node_attr_name.at(i)<<" , "<<value_map.size()<<std::endl;
				decoded_attributes.at("node_attr").emplace(std::make_pair(node_attr_name.at(i), value_map));
				value_map.clear();
			}

			//edges
			if(stdout>3) std::cout<<"edges"<<std::endl;
			value_map.clear();
			for(std::size_t i=0; i<num_edge_attr; i++){		
				//read the name
				line = read_word_binary(in, separator_1);
				if(stdout>3) std::cout<<line<<std::endl;
				edge_attr_name.emplace_back(line);
				// get the type
				in.get(type);
				// NO more // read until separator_2 is found
				// Get # of values to read (coded in 4 bytes for now)
				cont = get_size_t_from_bytes(in, 4);				
				for(std::size_t j = 0; j< cont; j++){
					switch(type){ 
			        case 'i': // int		        	
			            in.read(reinterpret_cast<char*>(&aux_int), sizeof(int));
			            aux_string = std::to_string(aux_int);
			            break;
			        case 'f': // float		        	
			            in.read(reinterpret_cast<char*>(&aux_float), sizeof(float));
			            aux_string = std::to_string(aux_float);
			            break;
			        case 'd': // double		        	
			            in.read(reinterpret_cast<char*>(&aux_double), sizeof(double));
			            aux_string = std::to_string(aux_double);
			            break;
			        case 'c': // char   
			            in.get(aux_char);
			            aux_string = aux_char;
			            break;	
			        case 's': // string
			        	aux_string = read_word_binary(in, separator_1);
			    	}
			    	if(stdout>4) std::cout<<aux_string<<std::endl;
					value_map.emplace(std::to_string(j), aux_string);
				}

				value_map.emplace(std::to_string(cont), "dummy");
				if(stdout>3) std::cout<<"insert "<<edge_attr_name.at(i)<<" , "<<value_map.size()<<std::endl;
				decoded_attributes.at("edge_attr").emplace(std::make_pair(edge_attr_name.at(i), value_map));
				value_map.clear();
			}

		}
		else{
			//nodes
			value_map.clear();
			for(std::size_t i=0; i<num_node_attr; i++){
				// attribute name
				std::getline(in,line,separator_1);
				if(stdout>3) std::cout<<line<<std::endl;
				node_attr_name.emplace_back(line);
				// number of values
				std::getline(in,line,separator_1);
				if(stdout>3) std::cout<<line<<std::endl;
				cont = std::stoi(line); 
				// Values
				for(std::size_t j=0; j < cont; j++){
					std::getline(in,line,separator_1);
					value_map.emplace(std::to_string(j), line);
				}
				value_map.emplace(std::to_string(cont), "dummy");
				decoded_attributes.at("node_attr").emplace(std::make_pair(node_attr_name.at(i), value_map));
				value_map.clear();
			}
			//edges
			value_map.clear();
			for(std::size_t i=0; i<num_edge_attr; i++){
				// attribute name
				std::getline(in,line,separator_1);
				if(stdout>3) std::cout<<line<<std::endl;
				edge_attr_name.emplace_back(line);
				// number of values
				std::getline(in,line,separator_1);
				if(stdout>3) std::cout<<line<<std::endl;
				cont = std::stoi(line); 
				// Values
				for(std::size_t j=0; j < cont; j++){
					std::getline(in,line,separator_1);
					value_map.emplace(std::to_string(j), line);
				}
				value_map.emplace(std::to_string(cont), "dummy");
				decoded_attributes.at("edge_attr").emplace(std::make_pair(edge_attr_name.at(i), value_map));
				value_map.clear();
			}	
		}
			

		ordered_attributes = get_ordered_attributes<std::map<std::string,std::string>>(decoded_attributes);
		attr_sizes.emplace(std::make_pair("node_attr",std::map<std::string, std::size_t>() ));
		attr_sizes.emplace(std::make_pair("edge_attr",std::map<std::string, std::size_t>() ));

		for(std::size_t i=0 ; i<node_attr_name.size(); i++){
			if(stdout>3) std::cout<<node_attr_name.at(i)<<" - "<<node_attr_size.at(i)<<std::endl;
			attr_sizes.at("node_attr").emplace(std::make_pair(node_attr_name.at(i), node_attr_size.at(i)));
		}
		for(std::size_t i=0 ; i<edge_attr_name.size(); i++){
			if(stdout>3) std::cout<<edge_attr_name.at(i)<<" - "<<edge_attr_size.at(i)<<std::endl;
			attr_sizes.at("edge_attr").emplace(std::make_pair(edge_attr_name.at(i), edge_attr_size.at(i)));
		}

		in.close();
	}
	else{
		std::cout<<"Unable to open info_file. Stopping decode execution"<<std::endl;
		throw compression_exception("Unable to open info_file. Stopping decode execution");
		
	}

	// Decode graphs

	if(stdout>1) std::cout<<"decode_collection: children structure"<<std::endl;
	// Get children structure
	std::map<std::size_t, std::vector<std::size_t>> children;
	std::size_t num_nodes = arborescence.size()+1; 
	for(std::size_t i=0; i<num_nodes; i++){
		children.emplace(std::make_pair(i, std::vector<std::size_t> ()));
	}	

	for(std::size_t i=0; i<arborescence.size(); i++){
		children.at(arborescence.at(i)).emplace_back(i);			
	}

	std::size_t root = arborescence.size();
	std::size_t parent_num;


	std::list<std::size_t> to_do;
	to_do.emplace_front(root);

	// Graphs are not decompressed in the same order they were read initially 
	// This maps the position in the environment to the actual graph id stored in the arborescence
	std::map<std::size_t, std::size_t> pos_to_id;

	std::ifstream graph_in;
	std::pair<ged::GEDGraph::GraphID, ged::GEDGraph::GraphID> limits;
	

	typename std::list<std::pair<std::pair<std::size_t, std::size_t>,UserEdgeLabel>>::iterator iter;
	std::pair<std::pair<std::size_t, std::size_t>,UserEdgeLabel> edge;


	std::string decomp_graph_name;
	while(!to_do.empty()){
		if(stdout>3) std::cout<<"START WHILE"<<std::endl;
		parent_num = to_do.front();
		to_do.pop_front();
		if(stdout>3) std::cout<<"star for"<<std::endl;
		for(auto const child : children.at(parent_num)){			
			to_do.emplace_front(child);
			if(stdout>2) std::cout<<"Parent: "<<parent_num<<", Child: "<<child<<std::endl;
			graph_file = path + "/" + graph_names.at(child);
			if(stdout>2) std::cout<<"File: "<<graph_file<<std::endl;
			if(stdout>2) std::cout<<"relaxed?: "<<relaxed<<std::endl;
			// Decompress graph
			decomp_graph_name = graph_names.at(child).substr(0,graph_names.at(child).find(".cmp"));

			decode_single_graph<UserNodeID,UserNodeLabel,UserEdgeLabel>(
				binary, relaxed,
					graph_file, decomp_graph_name, graph_classes.at(child),
					env, parent_num, child, root,
					pos_to_id, b_ni, b_ei,
					ordered_attributes,	attr_sizes, stdout, separator_1);
			

		}
	if(stdout>3) std::cout<<"END WHILE"<<std::endl;
	}
	
	// Decode environement
	if(stdout>2) 
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

	if(stdout>1) std::cout<<"Num graphs: "<<env.num_graphs()<<std::endl;
	if(stdout>1) std::cout<<"translate_env"<<std::endl;

	translate_env<UserNodeID,UserNodeLabel,UserEdgeLabel>(decoded_attributes, env, env_decoded, true, true, stdout);
	
	env = env_decoded;
	if(stdout>1) std::cout<<"END DECODE"<<std::endl;

}


}
#endif /* COMPRESS_SRC_ABC_IPP_ */
