// #include <iostream>
// #include <vector>
// #include <variant>
// #include <string>
// #include <map>

// using vec_name_comp = std::vector<std::variant<std::string, int>>;
// using map_var_vec = std::map<vec_name_comp, std::vector<double>>;

// using vec_name_comp_vd = std::vector<std::variant< std::string, int, std::vector<double> >>;

// int main(){
//   int i=0;
//   double d=1.1;
//   std::vector<std::variant<int, double>> vec{i,d};

//   std::cout << std::get<int>(vec[0]) << std::endl;
//   std::cout << std::get<double>(vec[1]) << std::endl;

//   vec_name_comp v{"test",1};
//   std::vector<double> vd={1.,2.};

//   map_var_vec m={{v,vd}};

//   vec_name_comp_vd M={{"test",1,vd}};
// };

/////////////////////////////////////////////////

// #include <iostream>
// #include <vector>
// #include <variant>
// #include <string>

// using SorI = std::variant<std::string, int>;
// using V_SorI = std::vector<SorI>;
// using VV_SorI = std::vector<V_SorI>;

// int main(){

//   VV_SorI vv{{"key1",2},{"key2","2"}};

//   for(const auto& v:vv){
//     auto str = std::get<std::string>(v[0]);
//     auto si = v[1];
//     auto i = si.index();
//     if(i==0){
//       std::cout << str << ", " <<  std::get<std::string>(si) << std::endl;
//     }else if(i==1){
//       std::cout << str << ", " << std::get<int>(si) << std::endl;
//     }
//   }

// };

/* ------------------------------------------------------ */
// 2021/07/11

// #include <iostream>
// #include <vector>
// #include <variant>
// #include <string>

// using Var_SID = std::variant<std::string, int, double>;
// using V_Var_SID = std::vector<Var_SID>;

// int main()
// {
//   V_Var_SID Vec_Var{"name", 3.2, 1.2, "number", 5., 1, 2, 3};

//   for (const auto &v : Vec_Var)
//   {
//     if (std::holds_alternative<int>(v))
//       std::cout << "int : " << std::get<int>(v) << std::endl;
//     if (std::holds_alternative<double>(v))
//       std::cout << "double : " << std::get<double>(v) << std::endl;
//     if (std::holds_alternative<std::string>(v))
//       std::cout << "string : " << std::get<std::string>(v) << std::endl;
//   }
// };

/* ------------------------------------------------------ */

#include <iostream>
#include <vector>
#include <variant>
#include <string>
#include <sstream>
#include <iomanip>

using Var_SID = std::variant<std::string, int, double>;
using V_Var_SID = std::vector<Var_SID>;
using VV_Var_SID = std::vector<std::vector<Var_SID>>;

std::string Grid(const V_Var_SID &V, const int i = 10)
{
  std::stringstream ss;
  for (const auto &var : V)
  {
    ss << std::setw(i) << std::right;
    if (std::holds_alternative<int>(var))
      ss << std::get<int>(var);
    else if (std::holds_alternative<double>(var))
      ss << std::get<double>(var);
    else if (std::holds_alternative<std::string>(var))
      ss << std::get<std::string>(var);
  }
  return ss.str();
};
std::string Grid(const VV_Var_SID &VV, const int i = 10)
{
  std::stringstream ss;
  for (const auto &V : VV)
  {
    ss << Grid(V);
    ss << std::endl;
  }
  return ss.str();
};

int main()
{
  V_Var_SID Vec_Var{"name", 3.2, 1.2, "number", 5., 1, 2, 3};
  std::cout << Grid(Vec_Var) << std::endl;
};
