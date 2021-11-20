////////////////////////////////////////////////////////////////

#include "fundamental.h"
#include "GNUPLOT.h"
#include <sstream>
int main(){


  // LoadObj obj;
  // obj.ConstructFromString("# The units used in this file are centimeters.\nv 2.229345 -0.992723 -0.862826\nv 2.292449 -0.871852 -0.882400\nv 2.410367 -0.777999 -0.841105\nv 2.407309 -0.974980 -0.805091\nv 2.539200 -0.727778 -0.750475");
  // std::cout << obj.JSON() << std::endl;
    
  ////////////////////////////////////////////////////////////////  
  GNUPLOT plot;

  std::string str = "/Users/tomoaki/Dropbox/cpp/source_test/obj/papercup.obj";
  // std::cout << StringSplit("  f 1/1/1 2/2/2 3/3/3 4/4/4  ", {"    ","   ","  "," ","/"}) << std::endl;

  // std::cout << str + "adfadd" << std::endl;

  
  // std::cout << ( str.find_first_not_of(" ")==std::string::npos/*no string except spaces*/) << std::endl;

  // std::vector<std::vector<std::string>> read_line;
  // Load(str, read_line, {"  ","/"});


  // for(const auto& s: read_line)
  //   std::cout << s << std::endl;
  // std::cout << (*(read_line.rbegin()+1)).size() << std::endl;
  // std::cout << (*(read_line.rbegin()+2)).size() << std::endl;
  // std::cout << (*(read_line.rbegin()+3)).size() << std::endl;        

  LoadObj test(str);

  //std::cout << test.f_v << std::endl;  
  
  std::vector<std::vector<std::vector<double>>> ret;

  for(auto m=0; m<test.f_v.size(); m++)
    {
      if(test.f_v[m].size()==3){
  	int max = 3;	
  	std::vector<std::vector<double>> tmp;
  	for(auto i=0; i<max; i++)
  	  {
  	    double x = (double)i/(max-1);
  	    int maxj=max-i;
  	    for(auto j=0; j<maxj; j++)
  	      {
  		double y = (double)j/(max-1);
  		tmp.push_back(test.surface(m,x, y));
  	      }
  	  }
  	//      ret.push_back(tmp);
  	plot.SaveSplotData({tmp},{{"w","l"},{"lc","rgb \"blue\""}});
      }

      if(test.f_v[m].size()==4){
  	int max = 4;
  	std::vector<std::vector<double>> tmp;
  	for(auto i=0; i<max; i++)
  	  {
  	    double x = (double)i/(max-1);
  	    int maxj=max;
  	    for(auto j=0; j<maxj; j++)
  	      {
  		double y = (double)j/(max-1);
  		tmp.push_back(test.surface(m,x, y));
  	      }
  	  }
  	plot.SaveSplotData({tmp},{{"w","l"},{"lc","rgb \"red\""}});      
      }
      
    }

  plot.Splot();
  // std::cout << test.v[test.f[0][0]] << std::endl;
  // plot.SaveSplotData({{test.v[test.f[0][0]]}},{{"w","p"}});
  // plot.Splot();
  
  std::cin.ignore(); 
  
  return 0;
}

// ////////////////////////////////////////////////////////////////

// #include "fundamental.h"
// #include "GNUPLOT.h"
// int main(){
  
//   GNUPLOT plot;

//   std::string str = "/Users/tomoaki/Dropbox/cpp/source_test/obj/cow.obj";
//   LoadObj test(str);
//   plot.SaveSplotData({test.v},{{"w","p"}});
//   plot.Splot();
  
//   std::cin.ignore();
  
//   return 0;
// }

//////////////////////////////////////////////////////////////////////////////////////////

// #include "fundamental.h"
// int main(){

//   std::vector<std::vector<int>> v{{1,2,3,4,5},{10,20,30,40,50}};  
//   std::cout << v << std::endl;
  
//   std::vector<std::vector<double>> u{{1,2,3,4,5},{10,20,30,40,50}};  
//   std::cout << u << std::endl;

//   std::vector<std::vector<std::vector<double>>> w{{{1,2,3,4,5},{10,20,30,40,50}},{{1,2,3,4,5},{10,20,30,40,50}}};  
//   std::cout << w << std::endl;
  
//   std::vector<std::vector<std::string>> str{{"1","2","3","4","5"},{"10","20","30","40","50"}};  
//   std::cout << str << std::endl;
  
//   return 0;
// }

// ////////////////////////////////////////////////////////////////

// #include "fundamental.h"
// #include "GNUPLOT.h"
// int main(){
  
//   GNUPLOT plot;

//   std::string str = "/Users/tomoaki/Dropbox/cpp/source_test/obj/cow.obj";
//   LoadObj test(str);
//   plot.SaveSplotData({test.v},{{"w","p"}});
//   plot.Splot();
  
//   std::cin.ignore();
  
//   return 0;
// }

//////////////////////////////////////////////////////////////////////////////////////////

// #include "fundamental.h"
// int main(){

//   std::string str = "/Users/tomoaki/Dropbox/cpp/source_test/obj/pumpkin.obj";
//   LoadObj test(str);
//   std::cout << test.v << std::endl;
  
//   return 0;
// }

//////////////////////////////////////////////////////////////////////////////////////////

// #include "fundamental.h"
// int main(){
//   std::string str = "/Users/tomoaki/Dropbox/cpp/source_test/obj/bunny.obj";
//   std::cout << StringSplit(str,{""}) << std::endl;  
//   return 0;
// }

///////////////////////////////////////////////////////////////////////////////////////////

// #include "fundamental.h"
// using namespace std;

// int main(){

//   std::string str = "/Users/tomoaki/Dropbox/cpp/source_test/obj/bunny.obj";
//   cout << StringSplit(str,{"/","."}) << endl;
  
//   return 0;

// }

///////////////////////////////////////////////////////////////////////////////////////////

// #include "fundamental.h"
// #include "GNUPLOT.h"

// using namespace std;

// // class triangle{
// // public:
// //   std::vector<double> xyz;
// //   double EPS=1E-15;
// //   int K;
// //   std::vector<double> q;//knot
// //   double itvl=1.;//interval  

// //   std::vector<double> C;
// //   std::vector<Bbasis*> VecBbasis;
// //   std::map<Bbasis*,Bbasis*> LinedBbasis;
// //   std::map<Bbasis*,Bbasis*> BifurcatedBbasis;  
// //   std::map<Bbasis*,std::vector<double>> Mapq;
// // };

// int main(){

//   std::string str = "/Users/tomoaki/Dropbox/cpp/source_test/obj/bunny.obj";
//   cout << StringSplit(str,{"/","."}) << endl;
  
//   return 0;

// }
