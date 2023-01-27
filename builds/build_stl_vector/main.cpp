#include <iostream>
#include <vector>

int main() {
   /*
   vectorのpushでサイズが変わる．
   占有するメモリ容量は変わらない
   */
   std::vector<int> v;
   for (auto i = 0; i < 10; i++) {
      v.emplace_back(i);
      std::cout << v.size() << std::endl;
   }
   for (auto i = 0; i < 10; i++) {
      v.pop_back();
      std::cout << v.size() << std::endl;
   }
   /* -------------------------------------------------------------------------- */
   for (auto i = 0; i < 10; i++) {
      v.insert(v.begin(), i);
      std::cout << v.size() << std::endl;
   }
   for (auto i = 0; i < 10; i++) {
      v.pop_back();
      std::cout << v.size() << std::endl;
   }
};