class test1_4{
public:
  double length(const std::vector<std::vector<double>>& vv){  
    return Norm(vv[0]-vv[1]);
  };
  double length(const networkLine* l){  
    return length(l->getLocations());
  };
  // SECURE DIVIDING CEHCK
  test1_4(){
    Print("**************** test1_4 ****************", Red);
    NetworkObj obj("./obj/tank.obj");
    NetworkW water({3,3},{1.1,1.1,1/2.},.0);
    GNUPLOT plot;
    //============================
    plot.Set({{"style","arrow 1 nohead lc \"blue\" lw 1"}});
    plot.Set({{"style","arrow 2 nohead lc \"magenata\" lw 1"}});    
    plot.Set({{"style","arrow 3 nohead lc \"red\" lw 2"}});
    for(auto i=0; i<14; i++){
      Network corssNetwork(obj,water);          
      Print(obj.Points.size(),Red);
      plot.Set({{"title","\""+std::to_string(i)+"\""}});      	  
      
      {
	std::vector<std::vector<std::vector<double>>> VVV;      
	std::vector<networkLine*> L;
	for(const auto& p:obj.Points[0]->getNeighbors_AvoidCross_recursive())
	  for(const auto& l:p->getLines())
	    L.emplace_back(l);
	pushVectorData(VVV,DeleteDuplicates(L));
	plot.SaveVectorData(VVV,{{"arrowstyle","1"},{"notitle",""}});
	for(const auto& l:takeInsideOfBounds(obj.getLines(),{{-0.55,0.55},{-0.55,0.55},{-100,100}})){
	  auto longline = longerLine(l);
	  if(length(longline)>(1.-i/15.))
	    obj.divide(longline);
	}
      }

      {
	std::vector<std::vector<std::vector<double>>> VVV;      
	std::vector<networkLine*> L;
	for(const auto& p:water.Points[0]->getNeighbors_AvoidCross_recursive())
	  for(const auto& l:p->getLines())
	    L.emplace_back(l);
	pushVectorData(VVV,DeleteDuplicates(L));
	plot.SaveVectorData(VVV,{{"arrowstyle","2"},{"notitle",""}});

	for(const auto& l:takeInsideOfBounds(water.getLines(),{{-0.55,0.55},{-0.55,0.55},{-100,100}})){
	  auto longline = longerLine(l);
	  if(length(longline)>(1.-i/15.))
	    water.divide(longline);
	}	
      }
      
      plot.SaveData(corssNetwork.getLocations(),{{"w","p"},{"lc","\"magenta\""},{"pt","7"},{"title","cross points"}});
      plot.SaveVectorData(getVectorData(corssNetwork.getLines()),{{"arrowstyle","3"},{"title","corssNetwork"}});
      
      plot.plot3d();
      std::cin.ignore();
      plot.Clear();      
    }
  };
  ~test1_4(){Print(CHECK+" test1_4",Red);};
};
