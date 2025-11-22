// stdafx.h
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <filesystem>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <random>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

int main() {
   std::vector<int> vec = {1, 2, 3, 4, 5};
   std::cout << "Vector elements: ";
   for (const auto& v : vec) {
      std::cout << v << " ";
   }
   std::cout << std::endl;

   std::map<std::string, int> myMap = {{"apple", 1}, {"banana", 2}};
   std::cout << "Map elements: ";
   for (const auto& pair : myMap) {
      std::cout << pair.first << ":" << pair.second << " ";
   }
   std::cout << std::endl;

   return 0;
}