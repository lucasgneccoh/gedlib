#define GXL_GEDLIB_SHARED
#include "../src/arborescence_based_compression.hpp"
#include <cstdlib>

/*
	Util functions
*/

class test_exception : public std::exception{
	public:
	    test_exception(const std::string& msg) : m_msg(msg)
	    {}

	   ~test_exception()
	   {
	        std::cout << "test_exception::~test_exception" << std::endl;
	   }

	   virtual const char* what() const throw () 
	   {
	        std::cout << "test_exception - what:" << std::endl;
	        return m_msg.c_str();
	   }

	   const std::string m_msg;
};

template<class T>
double mean(std::vector<T> x){
	T sum = 0;
	for(T v : x){
		sum += v;
	}
	return static_cast<double>(sum) / x.size();
}

template<class T>
bool is_sorted(std::vector<T> x){
	for(std::size_t i =0; i< x.size()-1; i++){
		if(x.at(i)>x.at(i+1)) return false;
	}
	return true;
}

template<class T>
T get_percentile(std::vector<T> x, double perc){
	if(! is_sorted<T>(x)) std::sort(x.begin(), x.end());
	std::size_t pos = static_cast<std::size_t>(perc * x.size());
	return(x.at(min(pos,x.size())));
}

template<class T>
void basic_stats_from_vector(std::vector<T> vec, std::string name, 
		std::vector<std::string> &headers, std::vector<std::string> &values){

		headers.emplace_back("min_" + name);
		values.emplace_back(std::to_string(*std::min_element(vec.begin(), vec.end())));

		headers.emplace_back("max_" + name);
		values.emplace_back(std::to_string(*std::max_element(vec.begin(), vec.end())));

		headers.emplace_back("mean_" + name);
		values.emplace_back(std::to_string(mean<int>(vec)));

		for(int i = 1; i <= 9; i++){
			headers.emplace_back("p" + std::to_string(i*10) + "_" + name);
			values.emplace_back(std::to_string(get_percentile<int>(vec, static_cast<double>(i)/10)));	
		}
		
}

/*
	Returns the sum of the sizes of all the files inside a folder
*/
std::size_t get_folder_size(std::string path, bool verbose = false){
	std::vector<std::string> files;
	ged::GED_ABC abc;
	std::size_t total_size;
	std::size_t f_size;
	
	files = abc.get_sorted_file_names(path, "");
	total_size = 0;
	for(auto file : files){
		if(file=="." || file =="..") continue;
		f_size = abc.get_file_size(path+ "/" +file);
		if(verbose) std::cout<<"\t"<<file<<": "<<f_size<<std::endl;
		total_size += f_size;
	}
	
	return total_size;

}


/*
	Very important to define here the type of each attribute of each collection to be used!
	This will change the behaviour of the binary encoding only
*/
std::map<std::string, std::map<std::string, char>>
get_dataset_attr_types(std::string dataset){
	std::map<std::string, std::map<std::string, char>> ans;
	ans.emplace(std::make_pair("node_attr", std::map<std::string, char>()));
	ans.emplace(std::make_pair("edge_attr", std::map<std::string, char>()));
	ans.emplace(std::make_pair("graph_attr", std::map<std::string, char>()));

	if(dataset=="acyclic"){
		ans.at("node_attr").emplace(std::make_pair("chem", 'i'));
		ans.at("edge_attr").emplace(std::make_pair("valence", 'i'));
		ans.at("graph_attr").emplace(std::make_pair("class", 'c'));
		return ans;
	}
	if(dataset=="AIDS"){
		ans.at("node_attr").emplace(std::make_pair("charge", 'i'));
		ans.at("node_attr").emplace(std::make_pair("chem", 'i'));
		ans.at("node_attr").emplace(std::make_pair("symbol", 's'));
		ans.at("node_attr").emplace(std::make_pair("x", 'f'));
		ans.at("node_attr").emplace(std::make_pair("y", 'f'));
		ans.at("edge_attr").emplace(std::make_pair("valence", 'i'));
		ans.at("graph_attr").emplace(std::make_pair("class", 'c'));
		return ans;
	}
	if(dataset=="Letter"){
		ans.at("node_attr").emplace(std::make_pair("x", 'f'));
		ans.at("node_attr").emplace(std::make_pair("y", 'f'));
		ans.at("graph_attr").emplace(std::make_pair("class", 'c'));
		return ans;
	}
	if(dataset=="mao"){
		ans.at("node_attr").emplace(std::make_pair("chem", 'i'));
		ans.at("edge_attr").emplace(std::make_pair("valence", 'i'));
		ans.at("graph_attr").emplace(std::make_pair("class", 'i'));
		return ans;
	}
	if(dataset=="Mutagenicity"){
		ans.at("node_attr").emplace(std::make_pair("chem", 's'));
		ans.at("edge_attr").emplace(std::make_pair("valence", 'i'));
		ans.at("graph_attr").emplace(std::make_pair("class", 's'));
		return ans;
	}
	if(dataset=="pah"){
		ans.at("node_attr").emplace(std::make_pair("chem", 'i'));
		ans.at("edge_attr").emplace(std::make_pair("valence", 'i'));
		ans.at("graph_attr").emplace(std::make_pair("class", 'i'));
		return ans;
	}

	if(dataset=="Protein"){
		ans.at("node_attr").emplace(std::make_pair("aaLength", 'i'));
		ans.at("node_attr").emplace(std::make_pair("sequence", 's'));
		ans.at("node_attr").emplace(std::make_pair("type", 'i'));
		ans.at("edge_attr").emplace(std::make_pair("distance0", 'd'));
		ans.at("edge_attr").emplace(std::make_pair("distance1", 'd'));
		ans.at("edge_attr").emplace(std::make_pair("frequency", 'i'));
		ans.at("edge_attr").emplace(std::make_pair("type0", 'i'));
		ans.at("edge_attr").emplace(std::make_pair("type1", 'i'));
		ans.at("graph_attr").emplace(std::make_pair("class", 'i'));
		return ans;
	}
	if(dataset=="pcba"){
		ans.at("node_attr").emplace(std::make_pair("element", 's'));
		ans.at("node_attr").emplace(std::make_pair("charge", 'i'));
		ans.at("node_attr").emplace(std::make_pair("aromatic", 's'));
		ans.at("node_attr").emplace(std::make_pair("hcount", 'i'));
		ans.at("edge_attr").emplace(std::make_pair("order", 'f'));
		ans.at("graph_attr").emplace(std::make_pair("class", 'c'));
		return ans;
	}

	if(dataset=="msts_float_w" || dataset=="msts_float_w_un"){
		ans.at("node_attr").emplace(std::make_pair("stock", 's'));
		ans.at("edge_attr").emplace(std::make_pair("w", 'f'));
		ans.at("graph_attr").emplace(std::make_pair("class", 'c'));
		return ans;
	}
	if(dataset=="msts_int_w" || dataset=="msts_int_w_un"){
		ans.at("node_attr").emplace(std::make_pair("stock", 's'));
		ans.at("edge_attr").emplace(std::make_pair("w", 'i'));
		ans.at("graph_attr").emplace(std::make_pair("class", 'c'));
		return ans;
	}
	if(dataset=="msts_no_w" || dataset=="msts_no_w_un"){
		ans.at("node_attr").emplace(std::make_pair("stock", 's')); // Very important to have nodes LABELED
		ans.at("graph_attr").emplace(std::make_pair("class", 'c'));
		return ans;
	}

	return ans;
		
}

/*
	Does the compression, decompression and writing of decompressed graphs
*/
void treat_dataset(std::map<std::string, std::string> &args){

	std::vector<std::string> headers;
	std::vector<std::string> values;
	std::vector<std::string> gxl_file_names;
	std::vector<std::string> graph_classes;

	int stdout = std::stoi(args.at("stdout"));

	std::string output_root = args.at("output_root");
	std::string file_preffix = args.at("file_preffix");


	bool binary = false;
	bool relaxed = true;
	bool decomp_only = false;
	bool tar_size_only = false;
	// WARNING: The choice of separators can impact the execution of the program. 
	// I suggest \n and ;
	char separator_1 = ';';
	char separator_2 = '@';
	if(args.count("binary_encoding")>0 && args.at("binary_encoding")=="true") binary = true;
	if(args.count("relaxed_compression")>0 && args.at("relaxed_compression")=="false") relaxed = false;
	if(args.count("decomp_only")>0 && args.at("decomp_only")=="true") decomp_only = true;
	if(args.count("test_mode")>0 && args.at("test_mode")=="tar_size") tar_size_only = true;
	if(args.count("separator_1")>0) separator_1 = args.at("separator_1")[0];
	if(args.count("separator_2")>0) separator_2 = args.at("separator_2")[0];

	std::map<std::string, std::map<std::string, char>> attr_types = get_dataset_attr_types(file_preffix);

	// Specify dataset and paths
	std::string folder_suffix = (binary)? "_bin":"_text";
	output_root = output_root + "/" + file_preffix;

	std::string dir = output_root + "/encoded" + folder_suffix;

	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env_decoded;
	ged::GED_ABC abc;
	std::size_t num_trials = std::stoi(args.at("num_trials"));

	double comp_time, decomp_time, write_time;
	
	auto start = std::chrono::high_resolution_clock::now();
	auto end = std::chrono::high_resolution_clock::now();
	ged::Seconds runtime;



	ged::Options::GXLNodeEdgeType node_type = attr_types.at("node_attr").size()>0 ? ged::Options::GXLNodeEdgeType::LABELED : ged::Options::GXLNodeEdgeType::UNLABELED;
	ged::Options::GXLNodeEdgeType edge_type = attr_types.at("edge_attr").size()>0 ? ged::Options::GXLNodeEdgeType::LABELED : ged::Options::GXLNodeEdgeType::UNLABELED;
	std::unordered_set<std::string>  irrelevant_node_attributes;
	std::unordered_set<std::string>  irrelevant_edge_attributes;

	std::cout<<"Working on : "<<file_preffix<<std::endl;

	int sys_ans = 0;
	std::size_t tar_size;

	for(std::size_t num_t = 0; num_t < num_trials; num_t++){
		std::cout<<"Trial #"<<num_t+1<< " of "<< num_trials<<std::endl;
		headers.clear();
		values.clear();
		if(!decomp_only && !tar_size_only){
		if(stdout>0) std::cout<<"**********    COMPRESS   ************"<<std::endl;
		try{
			
			start = std::chrono::high_resolution_clock::now();

			abc.set_omega(0.2);

			abc.compress_collection<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>(
				binary, relaxed,
				args.at("graph_dir"), 
				args.at("collection_file"), node_type, edge_type, 
				irrelevant_node_attributes, irrelevant_edge_attributes,
				attr_types,
				output_root + "/encoded" + folder_suffix,
				output_root,
				file_preffix, args, stdout, headers, values,
				separator_1, separator_2);

			
			end = std::chrono::high_resolution_clock::now();

			runtime = end - start;
			comp_time = runtime.count();
			std::cout<<"Comp :"<<comp_time<<" , "<<runtime.count()<<std::endl;

		}
		catch(const std::exception& e){ throw; }

		}	

		if(!tar_size_only){
		if(stdout>0) std::cout<<"************    DECOMPRESS   **************"<<std::endl;		
		try{		
			
			start = std::chrono::high_resolution_clock::now();

			env_decoded = ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> ();
			abc.decode_collection<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>(
				binary, relaxed,
				env_decoded,
				output_root + "/encoded" + folder_suffix, file_preffix, args, stdout,
				separator_1, separator_2);
			
			end = std::chrono::high_resolution_clock::now();

			runtime = end - start;
			decomp_time = runtime.count();
			std::cout<<"decomp :"<<decomp_time<<" , "<<runtime.count()<<std::endl;
		}
		catch(const std::exception& e){ throw; }

		if(!decomp_only){
		if(stdout>0) std::cout<<"**********    WRITE GRAPHS    ************"<<std::endl;	
		
		try{
			
			start = std::chrono::high_resolution_clock::now();

			std::pair<ged::GEDGraph::GraphID, ged::GEDGraph::GraphID> graph_ids;
			graph_ids = env_decoded.graph_ids();
			std::vector<std::string> gxl_file_names, graph_classes;
			
			for(std::size_t i=graph_ids.first; i<graph_ids.second; i++){
				gxl_file_names.emplace_back(env_decoded.get_graph_name(i));
				graph_classes.emplace_back(env_decoded.get_graph_class(i));
				env_decoded.save_as_gxl_graph(i,  output_root + "/decoded" + folder_suffix + "/" + env_decoded.get_graph_name(i));	
			}
			env_decoded.save_graph_collection(output_root + "/decoded" + folder_suffix + "/" + file_preffix + ".xml",  gxl_file_names,  graph_classes);
						
			end = std::chrono::high_resolution_clock::now();

			runtime = end - start;
			write_time = runtime.count();
			std::cout<<"Write :"<<write_time<<" , "<<runtime.count()<<std::endl;
		}
		catch(const std::exception& e){
			throw;
		}

		}

		if(!decomp_only){
			headers.emplace_back("comp_time");
			values.emplace_back(std::to_string(comp_time));	
		}
		
		headers.emplace_back("decomp_time");
		values.emplace_back(std::to_string(decomp_time));

		if(!decomp_only){
			headers.emplace_back("write_time");
			values.emplace_back(std::to_string(write_time));
		}

		if(stdout>0) std::cout<<"**********    COMPRESS WITH TAR    ************"<<std::endl;

		
		std::string folder_to_tar = "tar -cjf " + dir + ".tar.bz --directory=" + output_root + " " +  "encoded" + folder_suffix;

		std::string remove = "rm " + dir + ".tar.bz";
		sys_ans = std::system(remove.c_str());
		if (stdout>0 && sys_ans==0) std::cout<<"Old tar file removed"<<std::endl;
		sys_ans = std::system(folder_to_tar.c_str());
		
		if (sys_ans==0){
			if(stdout>0) std::cout<<"Collections compressed"<<std::endl;
			tar_size = abc.get_file_size(dir + ".tar.bz");
			headers.emplace_back("tar_compressed_size");
			values.emplace_back(std::to_string(tar_size));
		}
		else{
			std::cout<<"Possible error while using tar system call. Returned value: "<<sys_ans<<std::endl;
			headers.emplace_back("tar_compressed_size");
			values.emplace_back(std::to_string(-1));
		}

		}

		tar_size = abc.get_file_size(dir + ".tar.bz");
		headers.emplace_back("tar_compressed_size");
		values.emplace_back(std::to_string(tar_size));
		headers.emplace_back("tar_compressed_sys_ans");
		values.emplace_back(std::to_string(sys_ans));


		if(stdout>0) std::cout<<"**********    WRITE TEST RESULTS    ************"<<std::endl;	


		if(decomp_only || tar_size_only){
			headers.emplace_back("dataset");
			values.emplace_back(file_preffix);
			headers.emplace_back("sample_size");
			values.emplace_back(args.at("graph_sample_size"));

		}

		if(args.count("write_results")>0 && args.at("write_results")=="true"){

			std::ofstream output_file;
			std::string out_file_str = args.at("output_root") + "/" + args.at("output_results_file");
			output_file.open(out_file_str.c_str(), ios::out | ios::app);
			
			if(output_file.is_open()){
				if(args.at("first_iteration")=="true"){
					abc.write_to_file<std::string>(output_file, headers);
					args.at("first_iteration")="false";
				}
				abc.write_to_file<std::string>(output_file, values);		
				output_file.close();	
			}
			else{
				std::cout<<"Error when opening output file"<<std::endl;
			}
		}
	}

}

/*
	Get statistics of a dataset (mean, min, max, percentiles)
		num_nodes per graph
		num_edges per graph
		node_degre
		For each attribute, the number of unique values and the distribution of the frequencies of each attribute
*/
void get_statistics(std::vector<std::string> collection_files, std::vector<std::string> graph_dirs,
	std::vector<std::string> file_preffixes, std::string output_root
	){

	        
	std::string path_datasets;
	std::string path_attributes;
	std::vector<std::string> headers;
	std::vector<std::string> values;
	std::string collection;
	std::string graph_dir;
	std::string file_preffix;


	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
	ged::ExchangeGraph<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> g;
	ged::Seconds runtime;

	ged::GED_ABC abc;

	

	// Get statistics
	std::vector<int> num_nodes;
	std::vector<int> num_edges;
	std::vector<int> node_degrees;
	std::vector<int> aux;
	std::map<std::string,std::map<std::string,std::map<std::string,int>>> attr_freq;
	

	
	std::pair<std::pair<std::size_t, std::size_t>, ged::GXLLabel> edge;
	typename std::list<std::pair<std::pair<std::size_t, std::size_t>, ged::GXLLabel>>::iterator iter; 

	bool first = true;
	for(std::size_t file = 0; file<collection_files.size(); file++){

		collection = collection_files.at(file);
		graph_dir = graph_dirs.at(file);
		file_preffix = file_preffixes.at(file);

		attr_freq.clear();
		attr_freq.emplace(std::make_pair("node_attr",std::map<std::string,std::map<std::string,int>>()));
		attr_freq.emplace(std::make_pair("edge_attr",std::map<std::string,std::map<std::string,int>>()));

		std::cout<< "Start: " << file_preffix ;
		auto start = std::chrono::high_resolution_clock::now();

		path_datasets = output_root + "/" + "stats_datasets.csv";
		path_attributes = output_root + "/" + "stats_attributes.csv";

		headers.clear();
		values.clear();
		env = ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>();
		std::vector<ged::GEDGraph::GraphID> graph_ids(env.load_gxl_graphs(
			graph_dir, collection,
				ged::Options::GXLNodeEdgeType::LABELED, ged::Options::GXLNodeEdgeType::LABELED));
		
		for(std::size_t i : graph_ids){
			num_nodes.emplace_back(env.get_num_nodes(i));
			num_edges.emplace_back(env.get_num_edges(i));
			g = env.get_graph(i, false, true, true);
			// Get node degree information from adjacency lists
			for(const auto &l : g.adj_lists){
				node_degrees.emplace_back(l.size());
			}
			// Get attribute information
			// Nodes

			for(const auto &node_label : g.node_labels){
				for(const auto & attr: node_label){
					if(attr_freq.at("node_attr").count(attr.first)>0){
						if(attr_freq.at("node_attr").at(attr.first).count(attr.second)>0){
							attr_freq.at("node_attr").at(attr.first).at(attr.second)++;
						}
						else{
							attr_freq.at("node_attr").at(attr.first).emplace(std::make_pair(attr.second,1));	
						}

					}
					else{
						attr_freq.at("node_attr").emplace(std::make_pair(attr.first,std::map<std::string,int>()));
					}
				}

			}
			// Edges
			for(iter = g.edge_list.begin(); iter != g.edge_list.end(); iter ++ ){
				edge = (*iter);
				for(const auto & attr: edge.second){
					if(attr_freq.at("edge_attr").count(attr.first)>0){
						if(attr_freq.at("edge_attr").at(attr.first).count(attr.second)>0){
							attr_freq.at("edge_attr").at(attr.first).at(attr.second)++;
						}
						else{
							attr_freq.at("edge_attr").at(attr.first).emplace(std::make_pair(attr.second,1));	
						}

					}
					else{
						attr_freq.at("edge_attr").emplace(std::make_pair(attr.first,std::map<std::string,int>()));
					}
				}
			}
		}

		headers.emplace_back("dataset");
		values.emplace_back(file_preffix);

		headers.emplace_back("num_graphs");
		values.emplace_back(std::to_string(env.num_graphs()));

		// nodes
		basic_stats_from_vector(num_nodes, "num_nodes", headers, values);

		// edges
		basic_stats_from_vector(num_edges, "num_edges", headers, values);

		// node_degrees
		basic_stats_from_vector(node_degrees, "node_degrees", headers, values);

		if(first){
			abc.write_to_file(path_datasets, headers);
			
		}
		abc.write_to_file(path_datasets, values);
		

		
		num_nodes.clear();
		num_edges.clear();
		node_degrees.clear();

		
		for(const auto &attr: attr_freq.at("node_attr")){
			// Put info into vector
			aux.clear();
			headers.clear();
			values.clear();
			for(const auto &vals : attr.second){
				aux.emplace_back(vals.second);
			}
			headers.emplace_back("dataset");
			values.emplace_back(file_preffix);

			headers.emplace_back("type_attr");
			values.emplace_back("node_attr");

			headers.emplace_back("attr_name");
			values.emplace_back(attr.first);

			headers.emplace_back("unique_values");
			values.emplace_back(std::to_string(attr.second.size()));

			basic_stats_from_vector(aux, "freq", headers, values);

			if(first){
				abc.write_to_file(path_attributes, headers);
				first = false;
			}
			abc.write_to_file(path_attributes, values);

		}

		for(const auto &attr: attr_freq.at("edge_attr")){
			// Put info into vector
			aux.clear();
			headers.clear();
			values.clear();
			for(const auto &vals : attr.second){
				aux.emplace_back(vals.second);
			}
			headers.emplace_back("dataset");
			values.emplace_back(file_preffix);

			headers.emplace_back("type_attr");
			values.emplace_back("edge_attr");

			headers.emplace_back("attr_name");
			values.emplace_back(attr.first);

			headers.emplace_back("unique_values");
			values.emplace_back(std::to_string(attr.second.size()));

			basic_stats_from_vector(aux, "freq", headers, values);

			if(first){
				abc.write_to_file(path_attributes, headers);
				first = false;
			}
			abc.write_to_file(path_attributes, values);

		}

		auto end = std::chrono::high_resolution_clock::now();
		runtime = end - start;
		std::cout<< " ... Done (" << runtime.count() << " s)"<<std::endl;
	}

}

/*
	Given the xml files listing the collection and the directory with the graphs, calculate the total size of the 
	collection (gxl files + xml file) and the size of the compressed .tar.bz file
*/
void get_gxl_sizes(std::vector<std::string> collection_files, std::vector<std::string> graph_dirs,
	std::vector<std::string> file_preffixes, std::string output_root, std::string tar_path
	){

	        
	std::string path_datasets;
	std::vector<std::string> headers;
	std::vector<std::string> values;
	std::string collection;
	std::string graph_dir;
	std::string file_preffix;

	ged::GED_ABC abc;

	path_datasets = output_root + "/" + "original_gxl_sizes.csv";

	// Get statistics
	std::size_t total_size;
	std::list<std::pair<std::string, std::string>> graphs;
	typename std::list<std::pair<std::string, std::string>>::iterator iter;
	std::string g_name;
	bool first = true;
	for(std::size_t file = 0; file<collection_files.size(); file++){

		headers.clear();
		values.clear();

		collection = collection_files.at(file);
		graph_dir = graph_dirs.at(file);
		file_preffix = file_preffixes.at(file);
		total_size=0;
		graphs.clear();

		abc.read_xml_graph_collection(collection, graphs);
		// add xml collection file
		total_size += abc.get_file_size(collection);

		for(iter = graphs.begin(); iter != graphs.end(); iter++){
			g_name = (*iter).first;
			total_size += abc.get_file_size(graph_dir + "/" + g_name);
		}

		headers.emplace_back("dataset");
		values.emplace_back(file_preffix);

		headers.emplace_back("total_size");
		values.emplace_back(std::to_string(total_size));

		headers.emplace_back("tar_size");
		values.emplace_back(std::to_string(abc.get_file_size(tar_path + "/" + file_preffix + ".tar.bz")));


		if(first){
			abc.write_to_file(path_datasets, headers);
			first = false;			
		}
		abc.write_to_file(path_datasets, values);	
		
	}

}


/*
	Given the paths:
		- gxl_orig_path: Folder with the datasets in original format (xml file + gxl files)
		- separate_files_path_attr: Folder with datasets in "separate files" format using dictionaries for attributes
		- separate_files_path_no_attr: Folder with datasets in "separate files" format with attributes inside the label files (no dictionary)
		- encoded_path: Folder with the datasets compressed by the abc method (binary and text forms). There

	This function will go through the folders getting the size of each representation and its .tar.bz form, to create a table comparing each size.

*/
void create_table_sizes(std::string gxl_orig_path, std::string separate_files_path_attr, std::string separate_files_path_no_attr, std::string encoded_path,
	std::vector<std::string> datasets, std::string output_file){

	ged::GED_ABC abc;
	std::vector<std::string> headers;
	std::vector<std::string> values;
	bool first = true;

	for(auto ds : datasets){
		headers.clear();
		values.clear();

		headers.emplace_back("dataset");
		values.emplace_back(ds);

		// gxl size
		headers.emplace_back("gxl_orig");
		values.emplace_back(std::to_string(get_folder_size(gxl_orig_path + "/" + ds)));

		headers.emplace_back("gxl_orig_tar_bz");
		values.emplace_back(std::to_string(abc.get_file_size(gxl_orig_path + "/compressed/" + ds + ".tar.bz")));


		// separate files
		headers.emplace_back("sep_files_attr");
		values.emplace_back(std::to_string(get_folder_size(separate_files_path_attr + "/" + ds)));

		headers.emplace_back("sep_files_attr_tar_bz");
		values.emplace_back(std::to_string(abc.get_file_size(separate_files_path_attr + "/compressed/" + ds + ".tar.bz")));

		headers.emplace_back("sep_files_no_attr");
		values.emplace_back(std::to_string(get_folder_size(separate_files_path_no_attr + "/" + ds)));

		headers.emplace_back("sep_files_no_attr_tar_bz");
		values.emplace_back(std::to_string(abc.get_file_size(separate_files_path_no_attr + "/compressed/" + ds + ".tar.bz")));

		// encoded
		headers.emplace_back("encoded_bin");
		values.emplace_back(std::to_string(get_folder_size(encoded_path + "/" + ds + "/encoded_bin")));

		headers.emplace_back("encoded_text");
		values.emplace_back(std::to_string(get_folder_size(encoded_path + "/" + ds + "/encoded_text")));

		headers.emplace_back("encoded_bin_tar_bz");
		values.emplace_back(std::to_string(abc.get_file_size(encoded_path + "/" + ds + "/encoded_bin.tar.bz")));

		headers.emplace_back("encoded_text_tar_bz");
		values.emplace_back(std::to_string(abc.get_file_size(encoded_path + "/" + ds + "/encoded_text.tar.bz")));

		if(first){
			abc.write_to_file(output_file, headers);
			first = false;			
		}
		abc.write_to_file(output_file, values);	
	}
	
}

/*
	Converts a collection in xml + gxl format to the "separate files" format. You can choose whether to use dictionaries for attributes or not
*/
void write_in_separate_files(std::string collection_file, std::string graph_dir, 
	std::string dataset, std::string output_root,
	char sep_line = '\n', char sep_col = ';', bool attr_file = true){

	ged::GED_ABC abc;
	std::list<std::pair<std::string, std::string>> graphs;
	typename std::list<std::pair<std::string, std::string>>::iterator iter;
	// Get names and classes
	abc.read_xml_graph_collection(collection_file, graphs);
	

	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;

	std::unordered_set<std::string>  irrelevant_node_attributes;
	std::unordered_set<std::string>  irrelevant_edge_attributes;
	
	for(iter = graphs.begin(); iter != graphs.end(); iter ++){
		env.load_gxl_graph(graph_dir, (*iter).first, 
			ged::Options::GXLNodeEdgeType::LABELED, ged::Options::GXLNodeEdgeType::LABELED, 
			irrelevant_node_attributes, irrelevant_edge_attributes,
			ged::undefined(), (*iter).second);		
	}


	std::map<std::string, std::map<std::string, std::vector<std::string>>> distribution;
	std::map<std::string, std::map<std::string, std::set<std::string>>> alphabets;
	std::map<std::string, std::map<std::string, std::size_t>> attr_sizes;
	std::size_t b_ni, b_na, b_ei, b_ea;
	bool fast_node_translate=false, fast_edge_translate=false;

	// Start writing data
	std::string path = output_root + "/" + dataset;
	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env_coded;
	abc.get_graphs_structure<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>(env, distribution, alphabets, attr_sizes, b_ni, b_na, b_ei, b_ea, fast_node_translate, fast_edge_translate);	
	distribution.clear();

	std::map<std::string, std::vector<std::string>> ordered_attributes;
	ordered_attributes = abc.get_ordered_attributes<std::set<std::string>>(alphabets);

	if(attr_file){
		// Attributes in metadata files		
		std::map<std::string, std::map<std::string, std::map<std::string,std::string>>> encoded_attributes;
		encoded_attributes = abc.get_attribute_encoding(alphabets);

		// Translate
		abc.translate_env<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>(encoded_attributes, env, env_coded, fast_node_translate, fast_edge_translate, 0);

		// deallocate first env to save space (is this really deallocating?) The initial env is then "lost"
		env = ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>();

		std::ofstream node_attr, edge_attr;
		node_attr.open(path + "/" + dataset + ".node_attr");
		edge_attr.open(path + "/" + dataset + ".edge_attr");
		if(node_attr.is_open() && edge_attr.is_open()){

			// Nodes
			for(const auto & attr : ordered_attributes.at("node_attr")){
				// Attr info: name, num_values
				node_attr << attr << sep_col << alphabets.at("node_attr").at(attr).size() << sep_line;
				for(const auto & val : alphabets.at("node_attr").at(attr)){
					node_attr << val << sep_line;
				}
			}

			// Edges
			for(const auto & attr : ordered_attributes.at("edge_attr")){
				// Attr info: name, num_values
				edge_attr << attr << sep_col << alphabets.at("edge_attr").at(attr).size() << sep_line;
				for(const auto & val : alphabets.at("edge_attr").at(attr)){
					edge_attr << val << sep_line;
				}
			}

			node_attr.close();
			edge_attr.close();
		}
		else{
			std::cout<<"Error opening output attribute files"<<std::endl;
			return;
		}

	}
	else{
		env_coded = env;
	}
		

	
	std::ofstream graph_labels, graph_idx, node_labels, edges, edge_labels;
	graph_labels.open(path + "/" + dataset + ".graph_labels");
	graph_idx.open(path + "/" + dataset + ".graph_idx");
	node_labels.open(path + "/" + dataset + ".node_labels");
	edges.open(path + "/" + dataset + ".edges");
	edge_labels.open(path + "/" + dataset + ".edge_labels");
	
	std::size_t node_count = 1;
	std::size_t graph_cont = 0;
	std::size_t node;
	ged::ExchangeGraph<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> g;
	bool first;
	std::string value_no_attr = "?";
	typename std::list<std::pair<std::pair<std::size_t, std::size_t>,ged::GXLLabel>>::iterator edge_iter;

	if(graph_labels.is_open() && graph_idx.is_open() &&  node_labels.is_open() &&  edges.is_open() &&  edge_labels.is_open()){
		for(iter = graphs.begin(); iter != graphs.end(); iter ++){
			graph_labels << (*iter).second << sep_line; // Graph class

			g = env_coded.get_graph(graph_cont, false, false, true);

			graph_cont ++;

			// Node info
			for(node = 0; node < g.num_nodes ; node ++){
				graph_idx << graph_cont << sep_line;
				first = true;
				for(const auto & attr : ordered_attributes.at("node_attr")){
					if(!first) node_labels << sep_col;
					first = false; 
					if (g.node_labels.at(node).count(attr)>0) {
						node_labels << g.node_labels.at(node).at(attr);
					}
					else{
						node_labels << value_no_attr;	
					}
				}
				node_labels << sep_line;
			}

			// Edge info
			for(edge_iter = g.edge_list.begin(); edge_iter != g.edge_list.end() ; edge_iter ++){
				edges << (*edge_iter).first.first + node_count << sep_col << (*edge_iter).first.second + node_count << sep_line;
				first = true;
				for(const auto & attr : ordered_attributes.at("edge_attr")){
					if(!first) edge_labels << sep_col;
					first = false;
					if ((*edge_iter).second.count(attr)>0 && static_cast<std::size_t>(std::stoi((*edge_iter).second.at(attr))) < alphabets.at("edge_attr").at(attr).size()) {
						edge_labels << (*edge_iter).second.at(attr);
					}
					else{
						edge_labels << value_no_attr;	
					}
				}
				edge_labels << sep_line;
			}

			node_count += g.num_nodes;

		}	
	}
	else{
		std::cout<<"Error opening output files"<<std::endl;
		return;
	}

}

/*
	Print values on the screen:
		Kolmogorov complexities for two representations (upper triangular matrix and edge pairs)
		ABC estimated costs if compressed from empty graph
*/
void print_info(std::map<std::string, std::string> &args){


	std::string file_preffix = args.at("file_preffix");
	
	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
	ged::GED_ABC abc;
	
	std::map<std::string, std::map<std::string, std::vector<std::string>>> distribution;
	std::map<std::string, std::map<std::string, std::set<std::string>>> alphabets;
	std::map<std::string, std::map<std::string, std::size_t>> attr_sizes;
	std::size_t b_ni, b_na, b_ei, b_ea;
	bool fast_node_translate, fast_edge_translate;

	std::vector<ged::GEDGraph::GraphID> graph_ids(env.load_gxl_graphs(
		args.at("graph_dir"), args.at("collection_file"),
			ged::Options::GXLNodeEdgeType::LABELED, ged::Options::GXLNodeEdgeType::LABELED));


	abc.get_graphs_structure<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> (
		env, distribution, alphabets, attr_sizes,b_ni, b_na, b_ei, b_ea, 
		fast_node_translate,fast_edge_translate);


	std::cout<<file_preffix<<":   "<<env.num_graphs()<<" graphs"<<std::endl;
	std::cout<<"\tK Complexity a: "<< abc.base_compr_cost_triangular_matrix(env, b_ni, b_na, b_ei, b_ea,0) <<std::endl;
	std::cout<<"\tK Complexity b: "<< abc.base_compr_cost_edge_pairs(env, b_ni, b_na, b_ei, b_ea,0) <<std::endl;
	std::cout<<"\tFrom empty graph (abc) (relaxed): "<< abc.base_compr_cost_abc(env, b_ni, b_na, b_ei, b_ea,1) <<std::endl;
	std::cout<<"\tFrom empty graph (abc) (non-relaxed): "<< abc.base_compr_cost_abc(env, b_ni, b_na, b_ei, b_ea,0) <<std::endl;
	std::cout<<"----------------------------------------------------------------------"<<std::endl;	

}

std::string get_graph_dir(std::string ds){
	std::string base = "../../data/datasets/";
	if(ds=="Letter" ){
		return base + ds + "/MED";  
	}
	if(ds=="AIDS" || ds=="Mutagenicity" ||ds=="Protein"){
		return base + ds + "/data";
	}
	return base + ds;
	
	
}

std::string get_collection_file(std::string ds){
	return "../../data/collections/" + ds +".xml";
}


std::vector<std::string> split_string(std::string s, std::string d){
	std::vector<std::string> res;
	size_t pos = 0;
	std::string token;
	while ((pos = s.find(d)) != std::string::npos) {
	    res.emplace_back(s.substr(0, pos));
	    s.erase(0, pos + d.length());
	}
	res.emplace_back(s.substr(0, pos));
	return res;
}


int main(int argc, char* argv[]){



	std::map<std::string, std::string> args;
	short int param = 1; 

	// Handle the execution to do other tasks
	// "dataset_stats" to only compute the descriptive statistics of nodes, edges and attributes
	// "gxl_sizes" to compute the size of the collections in the original gxl format and its .tar.bz compressed file
	// "write_separate_files" to trasnform the collections into separate files format. Attention with the directory names that are fixed
	// "table_sizes" to create the table comparing sizes of original files, abc-compressed files, tar.bz files, etc
	// WARNING: again, paths and directory names are fixed for the moment
	// "print_only" to iterate over datasets and show information on the terminal. Does not run the compression algorithm
	// "tar_size" only looks for the encoded file compressed using tar.bz and gets its size
	// This was made to overcome the issue of the tar system call not working and thus getting the tar file size of the last and not the current iteration
	// This also supposes that the test is being executed using a bash file that first call the original compression, then from the bash the tar command is called, and then this part again to get the file size only

	// Any other value results in running the compression and decompression algorithms for all datasets and graph_sample_sizes
	args.emplace(std::make_pair("test_mode",argv[param++]));

	// stdout: Number indicating the amount of text shown on terminal	
	args.emplace(std::make_pair("stdout",argv[param++]));

	std::vector<std::string> datasets_names = split_string(argv[param++], ":");

	// Number of times each dataset and k_sample will be executed 
	args.emplace(std::make_pair("num_trials", argv[param++]));

	// Set to "true" to always add the edge between the graph and the following one (with the stock msts for example) 
	args.emplace(std::make_pair("path_structure",argv[param++]));


	args.emplace(std::make_pair("output_root",argv[param++]));

	args.emplace(std::make_pair("output_results_file",argv[param++]));

	args.emplace(std::make_pair("ged_method",argv[param++]));

	args.emplace(std::make_pair("ged_method_options",argv[param++]));

	args.emplace(std::make_pair("ged_method_refinement",argv[param++]));

	args.emplace(std::make_pair("ged_method_refinement_options",argv[param++]));

	args.emplace(std::make_pair("refinement_size",argv[param++]));
	

	args.emplace(std::make_pair("write_ged_matrix",argv[param++]));
	
	args.emplace(std::make_pair("write_arb",argv[param++]));

	args.emplace(std::make_pair("write_results",argv[param++]));

	// edit_cost_type: String telling if the modified compression costs should be used. 
	// "mod" meand transformed costs, any other value means the more "conservative" costs
	args.emplace(std::make_pair("edit_cost_type",argv[param++])); 

	// relaxed_compression: "true" if the compression format should allow simultaneous node insertions and deletions.
	// File sizes will bigger be in general 
	args.emplace(std::make_pair("relaxed_compression",argv[param++])); 
	
	// Parameter to control the density of the collection graph. 
	// Out degree of each node. Select parameter "graph_sample_type" to chose if
	// this parameter should be reas as a % of the total number of graphs or 
	// simply as the total degree
	// Example: 50 for 50% if "graph_sample_type"=="%", if not, then 50 means 
	// every graph has 50 edges starting at it.

	args.emplace(std::make_pair("graph_sample_size","100"));
	std::vector<std::string> graph_sample_sizes = split_string(argv[param++], ":");

	// "%" if the 
	args.emplace(std::make_pair("graph_sample_type", argv[param++]));

	// "true" for compressed files in binary. Any other value leads to compressed files in text mode
	args.emplace(std::make_pair("binary_encoding",argv[param++]));

	// "true" if the executions will only try the decompression part
	args.emplace(std::make_pair("decomp_only", argv[param++]));


	// "true" if nodes are uniquely identified and they should be matched according to attributes
	args.emplace(std::make_pair("match_node_map", argv[param++]));
	args.emplace(std::make_pair("match_node_map_by", argv[param++]));
	




    // Just to handle output writing
    args.emplace(std::make_pair("first_iteration", "true"));
    // Creating entries that will be used
    args.emplace(std::make_pair("collection_file", "fill"));
    args.emplace(std::make_pair("graph_dir", "fill"));
    args.emplace(std::make_pair("file_preffix", "fill"));
    

    std::vector<std::string> collection_files;
	std::vector<std::string> graph_dirs;
	for(const auto d : datasets_names) {
        collection_files.emplace_back(get_collection_file(d));
        graph_dirs.emplace_back(get_graph_dir(d));
    }

	// Only get the statistics. No compression
	if(args.at("test_mode") == "dataset_stats"){		
		get_statistics(collection_files, graph_dirs, datasets_names, args.at("output_root"));
		return 0;
	}

	std::string tar_path = "../data/orig_datasets_to_tar/compressed";
	// get the real size of original datasets
	if(args.at("test_mode") == "gxl_sizes"){
		get_gxl_sizes(collection_files, graph_dirs, datasets_names, args.at("output_root"), tar_path);
		return 0;
	}


	// Write in separate files format
	if(args.at("test_mode") == "write_separate_files"){
		for(const auto d : datasets_names){
			std::cout << d << "  1 ..." <<std::flush;
	        write_in_separate_files(get_collection_file(d), get_graph_dir(d), 
				d, args.at("output_root") + "/separate_files",
				'\n', ';', true);

	        std::cout <<"  1 ..." <<std::flush;

	        write_in_separate_files(get_collection_file(d), get_graph_dir(d), 
				d, args.at("output_root") + "/separate_files",
				'\n', ';', false);
	        std::cout<<" done" <<std::endl;
	    }		
		return 0;
	}

	
	std::string gxl_orig_path = "../data/orig_datasets_to_tar";
	std::string separate_files_path_attr = "../data/output/separate_files";
	std::string separate_files_path_no_attr = "../data/output/separate_files_2";
	std::string encoded_path = "../data/output";
	std::string output_file = "../data/output/table_sizes_plus.csv";
	std::vector<std::string> folders = datasets_names;
	if(args.at("test_mode") == "table_sizes"){
		create_table_sizes(gxl_orig_path, separate_files_path_attr, separate_files_path_no_attr, encoded_path, folders,  output_file);
		return 0;
	}
	
	// Test compression
    for(std::size_t ds = 0; ds < datasets_names.size(); ds ++) {
    	
        args.at("collection_file") = collection_files.at(ds);
        args.at("graph_dir") = graph_dirs.at(ds);
        args.at("file_preffix") = datasets_names.at(ds);
        

        if(args.at("test_mode")=="print_only"){
        	// Just get the info in the screen
        	print_info(args);
        }
        else{
	        for(const auto k_sample:graph_sample_sizes){
		        std::cout<<"*** START: "<<args.at("file_preffix")<<", k_sample: "<<k_sample<<"  "<<args.at("graph_sample_type")<<" ***"<<std::endl;
				try{
					args.at("graph_sample_size") = k_sample;				
					treat_dataset(args);
					args.at("first_iteration") = "false";	
					
				}
				catch(const std::exception& e){
					std::cout<<"main: Error during execution: "<<e.what()<<std::endl;
					return 1;
				}
				std::cout<<"*** END: "<<args.at("file_preffix")<<" ***"<<std::endl;
			}	
		}
		
    }
	return 0;

}