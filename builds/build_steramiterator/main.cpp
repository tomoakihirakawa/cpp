#include "fundamental.hpp"

#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <algorithm>

int main(){
  
  std::ifstream fin("./main.cpp");

  std::string a;

  //HERE
  
  std::string content( (std::istreambuf_iterator<char>(fin) ),
		       (std::istreambuf_iterator<char>()    ) );
  //HERE

  
  std::cout << StringSplit(content, {"//HERE","//HERE"})[1] << std::endl;

  
  // while(!fin.eof()){
  //   std::getline(fin,a);
  //   std::cout << a << std::endl;;
  // }

  // std::ofstream fout;
  // fout.open("./");
  
}
