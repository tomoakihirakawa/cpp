#include "fundamental.hpp"

// template <typename T>
// std::vector<T> Subsets(const std::vector<T> &set, const int n)
// {
//   std::vector<T> ret;

//   for(auto i=0; i<set.size();i++)
//   {

//   };
// }
// ;

int main()
{

  auto add = [](const std::vector<int> &set, const std::vector<std::vector<int>> &tmp, int s) {
    std::vector<std::vector<int>> ret;
    for (auto A : tmp)
    {
      for (auto j = s; j < set.size(); j++)
      {
        auto a = A;
        a.insert(a.end(), set[j]);
        ret.insert(ret.end(), a);
      };
    }
    return ret;
  };

  std::vector<int> set = {1, 2, 3};
  std::vector<std::vector<int>> ret;

  ret = add(set, {{}}, 0);
  Print(ret);

  ret = add(set, add(set, {{}}, 0), 1);
  Print(ret);  
};
