#include <iostream>
#include <numeric>
#include <unordered_set>
#include <vector>
#include "basic.hpp"

class MyClass {
  private:
   std::vector<int> data;
   bool on_device;

  public:
   MyClass(const std::vector<int> &d) : data(d), on_device(false) {}

   void update() {
#pragma acc enter data copyin(this->data)
#pragma acc update device(this->data)
      on_device = true;
   }

   bool find(int target, bool acc = false) {
      bool found = false;

      if (acc) {
#pragma acc parallel loop reduction(|| \
                                    : found)
         for (auto i = 0; i < this->data.size(); ++i) {
            if (this->data[i] == target)
               found = true;
         }
         std::cout << "openacc" << found << std::endl;
      } else {
         for (auto i = 0; i < data.size(); ++i) {
            if (data[i] == target)
               found = true;
         }
         std::cout << "not openacc" << found << std::endl;
      }
      return found;
   }
};

int main() {
   int N = 100000000;
   std::vector<int> data(N);
   std::iota(data.begin(), data.end(), 1);
   MyClass my_class(data);
   int target = N - 1;

   my_class.update();

   Timer timer;
   std::cout << "timer : " << timer() << std::endl;
   if (my_class.find(target, true)) {
      std::cout << "The target value " << target << " was found in the unordered set." << std::endl;
   } else {
      std::cout << "The target value " << target << " was not found in the unordered set." << std::endl;
   }
   std::cout << "timer : " << timer() << std::endl;
   if (my_class.find(target, false)) {
      std::cout << "The target value " << target << " was found in the unordered set." << std::endl;
   } else {
      std::cout << "The target value " << target << " was not found in the unordered set." << std::endl;
   }
   std::cout << "timer : " << timer() << std::endl;

   return 0;
}