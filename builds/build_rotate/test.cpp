
#include <iostream>
#include <vector>
#include <algorithm>

int main(){

  std::vector<double> vec{1.,203.,45.};

  std::rotate(vec.begin(),vec.begin()+1,vec.end());
  
  for(const auto& d:vec)
    std::cout << d << std::endl; 
  
};
