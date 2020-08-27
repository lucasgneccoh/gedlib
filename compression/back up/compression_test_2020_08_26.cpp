#include "../src/arborescence_based_compression.hpp"

void treat_dataset(std::map<std::string, std::string> &args){

	std::vector<std::string> headers;
	std::vector<std::string> values;
	std::vector<std::string> gxl_file_names;
	std::vector<std::string> graph_classes;

	int stdout = std::stoi(args.at("stdout"));

	std::string output_root = args.at("output_root");
	std::string file_preffix = args.at("file_preffix");
	std::string folder_for_encoded = "encoded";
	// Specify dataset
	output_root = output_root + "/" + file_preffix;

	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
	ged::GED_ABC abc;

	
	if(stdout>0) std::cout<<"**********    LOAD DATA   ************"<<std::endl;
	
	std::vector<ged::GEDGraph::GraphID> graph_ids(env.load_gxl_graphs(
		args.at("graph_dir"), args.at("collection_file"),
			ged::Options::GXLNodeEdgeType::LABELED, ged::Options::GXLNodeEdgeType::LABELED));


	if(stdout>0) std::cout<<"**********    COMPRESS   ************"<<std::endl;
	try{
		
		abc.set_omega(0.2);
		abc.compress_collection<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>(env,
				output_root, file_preffix, folder_for_encoded, args, stdout, headers, values);

		abc.compress_collection_from_empty<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>(env,
				output_root, file_preffix, folder_for_encoded + "_from_empty", args, stdout);
	}
	catch(const std::exception& e){ throw; }

	if(stdout>0) std::cout<<"************    DECOMPRESS   **************"<<std::endl;
	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env_decoded;
	try{		
		abc.decode_collection<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>(env_decoded,
				output_root, file_preffix, args, stdout);
	}
	catch(const std::exception& e){ throw; }


	if(stdout>0) std::cout<<"**********    WRITE GRAPHS    ************"<<std::endl;	
	
	try{
		std::pair<ged::GEDGraph::GraphID, ged::GEDGraph::GraphID> graph_ids;
		graph_ids = env_decoded.graph_ids();
		std::vector<std::string> gxl_file_names, graph_classes;
		
		for(std::size_t i=graph_ids.first; i<graph_ids.second; i++){
			gxl_file_names.emplace_back(env_decoded.get_graph_name(i));
			graph_classes.emplace_back(env_decoded.get_graph_class(i));
			env_decoded.save_as_gxl_graph(i,  output_root + "/decoded/" + env_decoded.get_graph_name(i));	
		}
		env_decoded.save_graph_collection(output_root + "/decoded/" + file_preffix + ".xml",  gxl_file_names,  graph_classes);
		
	}
	catch(const std::exception& e){
		throw;
	}


	if(stdout>0) std::cout<<"**********    WRITE TEST RESULTS    ************"<<std::endl;	

	if(args.count("write_results")>0 && args.at("write_results")=="true"){
		std::ofstream output_file;
		output_file.open(args.at("output_results_file").c_str(), ios::out | ios::app);
		
		if(output_file.is_open()){
			if(args.at("first_iteration")=="true"){
				abc.write_to_file<std::string>(output_file, headers);			
			}
			abc.write_to_file<std::string>(output_file, values);		
			output_file.close();	
		}
		else{
			std::cout<<"Error when opening output file"<<std::endl;
		}
	}
	

}

int main(int argc, char* argv[]){


	std::map<std::string, std::string> args;
	
	// LLenar args y lanzar

	//std::cout<<"main: START, num of args: "<<argc<<" - "<<argv[argc-1]<<std::endl;
	args.emplace(std::make_pair("stdout",argv[1]));

	args.emplace(std::make_pair("collection_file",argv[2]));
	args.emplace(std::make_pair("graph_dir",argv[3]));
	args.emplace(std::make_pair("output_root",argv[4]));
	args.emplace(std::make_pair("file_preffix",argv[5]));

	args.emplace(std::make_pair("output_results_file",argv[6]));

	args.emplace(std::make_pair("ged_method",argv[7]));
	args.emplace(std::make_pair("ged_method_options",argv[8]));

	args.emplace(std::make_pair("graph_sample_size",argv[9]));
	
	args.emplace(std::make_pair("ged_method_refinement",argv[10]));
	args.emplace(std::make_pair("ged_method_refinement_options",argv[11]));
	args.emplace(std::make_pair("refinement_size",argv[12]));

	args.emplace(std::make_pair("write_ged_matrix",argv[13]));
	args.emplace(std::make_pair("write_arb",argv[14]));
	args.emplace(std::make_pair("write_results",argv[15]));

	// Not used for now
	args.emplace(std::make_pair("train_set",argv[16]));
	args.emplace(std::make_pair("train_path",argv[17]));
	args.emplace(std::make_pair("ring_method",argv[18]));



	args.emplace(std::make_pair("edit_cost_type",argv[19]));
	args.emplace(std::make_pair("relaxed_compression",argv[20]));
	
	args.emplace(std::make_pair("k_sample_file",argv[21]));

	

	//std::cout<<"main: files"<<std::endl;
	std::ifstream in_file_collections(args.at("collection_file").c_str());
	std::ifstream in_file_graphs(args.at("graph_dir").c_str());
	std::ifstream in_file_output_root(args.at("output_root").c_str());
	std::ifstream in_file_dataset(args.at("file_preffix").c_str());
	
	std::ifstream in_file_k_sample(args.at("k_sample_file").c_str());
	
	std::string input_collection_file;
    std::string input_graph_dir;
    std::string input_output_root;
    std::string input_dataset;

	if (!in_file_k_sample) {
        std::cout << "main: Unable to open k_sample_file"<<endl;
        return(1); // terminate with error
    }
    std::vector<std::string> graph_sample_sizes;
    while (in_file_k_sample >> input_collection_file) {
    	graph_sample_sizes.emplace_back(input_collection_file);
    }

    if (!in_file_collections || !in_file_graphs || !in_file_output_root || !in_file_dataset) {
        std::cout << "main: Unable to open files"<<endl;
        return(1); // terminate with error
    }

    

    args.emplace(std::make_pair("first_iteration", "true"));
    //std::cout<<"main: start while"<<std::endl;
    while (in_file_collections >> input_collection_file) {
        in_file_graphs >> input_graph_dir;
        in_file_output_root >> input_output_root;
        in_file_dataset >> input_dataset;


        args.at("collection_file") = input_collection_file;
        args.at("graph_dir") = input_graph_dir;
        args.at("output_root") = input_output_root;
        args.at("file_preffix") = input_dataset;

        for(const auto k_sample:graph_sample_sizes){
	        std::cout<<"*** START: "<<input_dataset<<", k_sample: "<<k_sample<<" % ***"<<std::endl;
			try{
				args.at("graph_sample_size") = k_sample;
				treat_dataset(args); 	
			}
			catch(const std::exception& e){
				std::cout<<"main: Error during execution: "<<e.what()<<std::endl;
				return 1;
			}
			
			std::cout<<"*** END: "<<input_dataset<<" ***"<<std::endl;
			args.at("first_iteration") = "false";
		}
		
    }
   
	return 0;
	
}