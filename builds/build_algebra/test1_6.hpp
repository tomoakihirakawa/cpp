class test1_6{
  using VVV_double = std::vector<std::vector<std::vector<double>>>;
public:
  template <class T>
  std::vector<T*> takeIfIntersect(std::vector<T*> obj){
    std::vector<T*> ret;
    for(const auto& f:obj){if(f->intersectQ()){ret.emplace_back(f);}}	
    return ret;
  };

  double length(const std::vector<std::vector<double>>& vv){  
    return Norm(vv[0]-vv[1]);
  };
  double length(const networkLine* l){  
    return length(l->getLocations());
  };
  // SECURE DIVIDING CHECK
  test1_6(){
    Print("**************** test1_6 ****************", Red);
    NetworkObj obj("./obj/tank.obj");
    NetworkW water({2,2},{1.1,1.1,1/2.},.2);
    GNUPLOT plot;
    //============================
    for(auto i=0; i<16; i++){
      plot.Set({{"style","arrow 1 nohead lc \"blue\" lw 1"}});
      plot.Set({{"style","arrow 2 nohead lc \"magenta\" lw 1"}});    
      plot.Set({{"style","arrow 3 nohead lc \"red\" lw 2"}});
      plot.Set({{"style","arrow 4 nohead lc \"green\" lw 2 dt 2"}});
      plot.Set({{"style","arrow 5 nohead lc \"orange\" lw 3 dt 2"}});
      plot.Set({{"style","arrow 10 nohead lc \"red\" lw .5"}});            
      plot.Set({{"key",""}});
      
      Print(obj.Points.size(),Red);
      plot.Set({{"title","\""+std::to_string(i)+"\""}});
      
      for(const auto& l:obj3D::takeInsideOfBounds(obj.getLines(),{{-0.55,0.55},{-0.55,0.55},{-100,100}})){
	auto longline = longerLine(l);
	if(length(longline)>(1.-i/15.))
	  obj.divide(longline);
      }
      for(const auto& l:obj3D::takeInsideOfBounds(water.getLines(),{{-0.55,0.55},{-0.55,0.55},{-100,100}})){
	auto longline = longerLine(l);
	if(length(longline)>(1.-i/15.))
	  water.divide(longline);
      }

      obj.displayStates();
      water.displayStates();      
      Network corssNetwork(obj,water);
      
      {
      	for(const auto& f:takeIfIntersect(obj.Faces)){
      	  for(const auto& ps:f->getPointsCutFaces()){
	    
	    VVV_double VVV;		  
	    Print(ps,Green);
	    auto s = ps.size();
	    
	    // for(auto i=0; i<s; i++){
	    //   VVV.push_back({ps[i]->getX(),ps[(i+1)%s]->getX() - ps[i]->getX()});
	    //   Print(f->parameterize(ps[i]->getX()),Red);
	    // }
	    // plot.SaveVectorData(VVV,{{"notitle","getARoute2()"}});	


	    for(const auto& route:f->getPointsCutFaces()){
	      VVV.clear();
	      auto indices = triangulate(route,f->getNormal());
	      for(const auto& index:indices){
		VVV.push_back({route[index[0]]->getX(),route[index[1]]->getX()-route[index[0]]->getX()});
		VVV.push_back({route[index[1]]->getX(),route[index[2]]->getX()-route[index[1]]->getX()});	      
	      }
	      plot.SaveVectorData(VVV,{{"lw",".5"},{"notitle",""}});
	    }
	  }	  
      	}
      }
            
      plot.SaveData(corssNetwork.getLocations(),{{"w","p"},{"lc","\"magenta\""},{"pt","7"},{"title","cross points"}});
      plot.SaveVectorData(getVectorData(corssNetwork.getLines()),{{"title","corssNetwork"}});
      
      plot.plot3d();
      std::cin.ignore();
      plot.Clear();      
    }
  };
  ~test1_6(){Print(CHECK+" test1_6",Red);};
};
