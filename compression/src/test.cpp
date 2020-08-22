//Test
#include <iostream>
#include <fstream>


int main(int argc, char* argv[]){
	std::string filename = "binary_file.bin";
	std::ofstream out_file;
	out_file.open(filename.c_str(),ios::binary);
	std::vector<int> vec;
	for(int i =0; i<5; i++) vec.emplace_back(i);
	
	if(out_file.is_open()){
		for(int i =0; i<5; i++) out_file << vec.at(i);	
	}
	out_file.close();

}