class test1_3{
public:
  double length(const std::vector<std::vector<double>>& vv){  
    return Norm(vv[0]-vv[1]);
  };
  double length(const networkLine* l){  
    return length(l->getLocations());
  };
  // SECURE DIVIDING CEHCK
  test1_3(){
    Print("**************** test1_3 ****************", Red);
    NetworkObj obj("./obj/tank.obj");
    NetworkW water({3,3},{1.12,1.15,1/2.},.0);  
    GNUPLOT plot;
    //============================
    plot.Set({{"style","arrow 1 nohead lc \"blue\" lw 1"}});
    std::vector<std::vector<std::vector<double>>> VVV;
    for(auto i=0; i<17; i++){
      Print(obj.Points.size(),Red);
      plot.Set({{"title","\""+std::to_string(i)+"\""}});      	  
      auto L = obj.getLines();      
      pushVectorData(VVV,L);	    

      for(const auto& l:obj3D::takeInsideOfBounds(obj.getLines(),{{-0.55,0.55},{-0.55,0.55},{-100,100}})){
	auto longline = longerLine(l);
	if(length(longline)>(1.-i/15.)){
	  obj.divide(longline);
	}
      }
      
      plot.SaveVectorData(VVV,{{"arrowstyle","1"},{"notitle",""}});
      plot.plot3d();
      //      std::cin.ignore();
      plot.Clear();      
      VVV.clear();
    }
  };
  ~test1_3(){Print(CHECK+" test1_3",Red);};
};
