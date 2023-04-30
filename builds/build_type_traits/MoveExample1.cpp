#include <chrono>
#include <iostream>
#include <vector>

int main() {
   std::vector<double> v1(20000 * 20000, 42);  // A large vector

   auto start = std::chrono::high_resolution_clock::now();
   auto end = std::chrono::high_resolution_clock::now();

   start = std::chrono::high_resolution_clock::now();
   std::vector<double> v2 = v1;  // Copying the vector
   end = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double> copy_time = end - start;

   start = std::chrono::high_resolution_clock::now();
   std::vector<double> v3 = std::move(v1);  // Moving the vector
   end = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double> move_time = end - start;

   std::cout << "Copy time: " << copy_time.count() << " seconds" << std::endl;
   std::cout << "Move time: " << move_time.count() << " seconds" << std::endl;

   return 0;
}
