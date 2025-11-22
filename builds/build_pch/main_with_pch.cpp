
/*

```shell
g++-14 -std=c++17 -x c++-header -o stdafx.h.gch stdafx.h -Ofast
```

```shell
time g++-14 -std=c++17 main_without_pch.cpp -Ofast
```

g++-14 -std=c++17 main_without_pch.cpp -Ofast  0.56s user 0.10s system 109% cpu 0.600 total

```shell
time g++-14 -std=c++17 main_with_pch.cpp -Ofast
```
g++-14 -std=c++17 main_with_pch.cpp -Ofast  0.24s user 0.09s system 115% cpu 0.284 total



*/

#include "stdafx.h"  // プリコンパイル済みヘッダーをインクルード

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