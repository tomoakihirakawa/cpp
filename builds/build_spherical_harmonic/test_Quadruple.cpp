#include "basic.hpp"
#include "basic_IO.hpp"

/*DOC_EXTRACT 0_0_0_translation_of_a_multipole_expansion

# ４倍精度の実装

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Debug ../ -DSOURCE_FILE=test_Quadruple.cpp
make
./test_Quadruple
```

*/

int main() {

   std::array<std::array<std::float128_t, 3>, 3> A = {{{1.0q, 2.0q, 3.0q},
                                                       {4.0q, 5.0q, 6.0q},
                                                       {7.0q, 8.0q, 9.0q}}};
   std::array<std::array<std::float128_t, 3>, 3> B = {{{M_PIq, 8.0q, 7.0q},
                                                       {6.0q, 5.0q, 4.0q},
                                                       {3.0q, 2.0q, 1.0q}}};
   auto C = A * B;
   std::cout << "A * B = " << C << std::endl;
}