#include "fundamental.hpp"

int main()
{

  int i0 = 0, i1 = 1, i2 = 2, i3 = 3;
  int *i0p = &i0;
  int *i1p = &i1;
  int *i2p = &i2;
  int *i3p = &i3;

  std::vector<int *> vi0{i0p, i1p, i2p, i1p};
  std::vector<std::vector<int *>> vi1{{i3p, i3p, i2p, i2p},{i3p, i3p, i2p, i2p},{i1p, i3p, i2p, i0p}};
  std::vector<int *> vi2;

  std::cout << Chain(vi0, vi1, vi2) << std::endl;
  std::cout << vi0 << std::endl;
  std::cout << vi1 << std::endl;
  std::cout << vi2 << std::endl;

};
