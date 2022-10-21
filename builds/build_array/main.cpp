#include <array>
#include <iostream>

int main() {

	std::array<double, 5> arr{1,2,3,4,5};
  for(auto i=0; i<arr.size();i++)
  {
    std::cout << arr[i] << std::endl;
  }
};
