//Test
#include<iostream>

using namespace std;

int main(int argc, char* argv[]){
	string str1 = "#x:empty";
	std::size_t pos=0;
	pos = str1.find(":");
	string str2;
	string str3;

	str2 = str1.substr(1, pos-1);
	str3 = str1.substr(pos+1);

	cout<<"str2: "<<str2<<endl;
	cout<<"str3: "<<str3<<endl;

}