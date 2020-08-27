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
	std::size_t num_trials = std::stoi(args.at("num_trials"));

	std::cout<<"Working on : "<<file_preffix<<std::endl;

	if(stdout>0) std::cout<<"**********    LOAD DATA   ************"<<std::endl;
		
	std::vector<ged::GEDGraph::GraphID> graph_ids(env.load_gxl_graphs(
		args.at("graph_dir"), args.at("collection_file"),
			ged::Options::GXLNodeEdgeType::LABELED, ged::Options::GXLNodeEdgeType::LABELED));


	for(std::size_t i = 0; i < num_trials; i++){
		std::cout<<"Trial #"<<i<< " of "<< num_trials<<std::endl;


		if(stdout>0) std::cout<<"**********    COMPRESS   ************"<<std::endl;
		try{
			
			abc.set_omega(0.2);
			abc.compress_collection<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>(env,
					output_root, file_preffix, folder_for_encoded, args, stdout, headers, values);

			if(args.count("test_mode")>0 && args.at("test_mode")=="true"){
				// skip the writing of the files
			}
			else{
				abc.compress_collection_from_empty<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>(env,
						output_root, file_preffix, folder_for_encoded + "_from_empty", args, stdout);
			}
		}
		catch(const std::exception& e){ throw; }

		if(args.count("test_mode")>0 && args.at("test_mode")=="true"){
			// skip the writing of the files
		}
		else{

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
		}
	

		if(stdout>0) std::cout<<"**********    WRITE TEST RESULTS    ************"<<std::endl;	

		if(args.count("write_results")>0 && args.at("write_results")=="true"){

			std::ofstream output_file;
			output_file.open(args.at("output_results_file").c_str(), ios::out | ios::app);
			
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

int main(int argc, char* argv[]){


	std::map<std::string, std::string> args;
	
	// Define values for the experiment
	// For an easier and more flexible testing, the tests are executed from a bash file including the
	// variables in the following order

	// stdout: Number indicating the amount of text shown on terminal	
	args.emplace(std::make_pair("stdout",argv[1]));

	// collection_file: File listing the paths to the collection files (.xml)
	// There will be one execution (test) for each collection path in the file
	// Ex of collection path: /home/usr/gedlib/data/collections/mao.xml
	args.emplace(std::make_pair("collection_file",argv[2]));

	// graph_dir: File listing the directories containing the graphs (.gxl). 
	// Each line corresponds to the line of same number in the collection_file
	// Ex of graph_dir: /home/usr/gedlib/data/datasets/mao
	args.emplace(std::make_pair("graph_dir",argv[3]));
	
	// output_root: File listing the directoy where to put all the resulting compressed graphs, metadata files and csv
	// files with information of the testing (if demanded).
	// Ex of output_root: /home/usr/gedlib/compression/data/output
	args.emplace(std::make_pair("output_root",argv[4]));

	// file_preffix: File listing the strings representing the dataset and test done. 
	// It will be used for naming the metadata file
	// and for choosing the output folder inside output_root  
	// Ex of file_preffix: For the mao dataset, the preffix would be mao. This would imply that the results
	// will be organized as follows:
	// 		output_root/mao/encoded -> for the encoded graphs
	// 		output_root/mao/encoded_from_empty -> for the encoded graphs from the empty graph to do a comparisson
	// 		output_root/mao/decoded -> for the decompressed graphs. This is to test if the decompression yields isomorphic graphs
	//      output_root/mao/mao.info_file -> file with metadata
	//		In this folder you will find the cost matrices and arborescences if the outputs were demanded
	args.emplace(std::make_pair("file_preffix",argv[5]));

	// output_results_file: path to a file that will store the results of the test, in csv format.
	// Ex: /home/usr/gedlib/compression/data/output/results_compression.csv
	args.emplace(std::make_pair("output_results_file",argv[6]));

	// ged_method: String telling the method to be used for GED calculation
	// For the moment the only tested method is "branch_uniform". Other options are "branch_fast" and "ipfp"
	args.emplace(std::make_pair("ged_method",argv[7]));

	// ged_method_options: This variable is supposed to carry the options to initialize the method with.
	// For the moment it contains only the number of threads to work with
	args.emplace(std::make_pair("ged_method_options",argv[8]));
	

	// ged_method_refinement: String telling the method to be used for GED calculation in the refinement stage
	// For the moment the only tested method is "branch_uniform". Other options are "branch_fast" and "ipfp"
	args.emplace(std::make_pair("ged_method_refinement",argv[9]));

	// ged_method_refinement_options: This variable is supposed to carry the options to initialize the method with.
	// For the moment it contains only the number of threads to work with
	args.emplace(std::make_pair("ged_method_refinement_options",argv[10]));

	// refinement_size: Number of steps up of each node to recalculate in the arborescence
	args.emplace(std::make_pair("refinement_size",argv[11]));

	// write_ged_matrix: true or false depending on wether the cost matrices should be written as output in csv format
	args.emplace(std::make_pair("write_ged_matrix",argv[12]));

	// write_arb: true or false depending on wether the arborescences should be written as output in csv format
	args.emplace(std::make_pair("write_arb",argv[13]));

	// write_results: true or false depending on wether the results of the test should be written as output in csv format
	args.emplace(std::make_pair("write_results",argv[14]));

	// edit_cost_type: String telling if the modified compression costs should be used. 
	// "mod" meand transformed costs, any other value means the more "conservative" costs
	args.emplace(std::make_pair("edit_cost_type",argv[15])); 

	// relaxed_compression: "true" if the compression format should allow simultaneous node insertions and deletions.
	// File sizes will bigger be in general 
	args.emplace(std::make_pair("relaxed_compression",argv[16])); 
	
	// k_sample_file: file containing the different sample sizes to test. They are expressed in percentage but as an integer
	// Example: 50 for 50%
	args.emplace(std::make_pair("k_sample_file",argv[17]));

	args.emplace(std::make_pair("num_trials", argv[18]));

	// Create encoding files??
	args.emplace(std::make_pair("test_mode",argv[19]));



	
	// To get only the values, without creating the actual files
	



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
        std::cout << "main: Unable to open some of the files"<<endl;
        return(1); // terminate with error
    }

    
    args.emplace(std::make_pair("graph_sample_size","100"));
    args.emplace(std::make_pair("first_iteration", "true"));
    


	std::vector<std::string> edit_cost_type = {"trad", "mod"};
	
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

				if(args.at("test_mode")=="true"){
					for(const auto edit : edit_cost_type){
						args.at("edit_cost_type") = edit;
						if(edit == "mod"){
							args.at("relaxed_compression") = "false";
						} 
						else{
							args.at("relaxed_compression") = "true";
						}
						treat_dataset(args); 
						args.at("first_iteration") = "false";
					}	
				}
				else{
					treat_dataset(args);
					args.at("first_iteration") = "false";
				}

				

				
			}
			catch(const std::exception& e){
				std::cout<<"main: Error during execution: "<<e.what()<<std::endl;
				return 1;
			}
			
			std::cout<<"*** END: "<<input_dataset<<" ***"<<std::endl;
			
		}
		
    }
   
	return 0;
	
}