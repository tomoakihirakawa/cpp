class test1_2{
public:
  double length(std::vector<std::vector<double>> vv){  
    return Norm(vv[0]-vv[1]);
  };
  double length(const networkLine* l){  
    return length(l->getLocations());
  };
  
  test1_2(){
    Print("**************** test1_2 **************** ",Red);
    NetworkObj obj("./obj/tank.obj");
    NetworkW water({3,3},{1.12,1.15,1/2.},.0);  
    GNUPLOT plot;
    //============================
    for(auto i=0; i<10; i++)
      for(const auto& l:obj.getLines())
	if(longerLine(l) && length(l)>(1.5-i/10.)){
	  obj.divide(l);
	}
	
    //=========== plot ===========
    plot.Set({{"xrange","[-.8:.8]"},{"yrange","[-.8:.8]"},{"zrange","[0:1.2]"}});
    plot.Set({{"style","arrow 1 nohead lc \"blue\" lw 1"}});
    plot.Set({{"style","arrow 2 nohead lc 2 lw 1"}});
    plot.Set({{"style","arrow 3 nohead lc \"red\" lw 2"}});
    for(auto i=0; i<15; i++){//ここが増えるとより多くcrossをチェックするため，crossの周辺だけをより細かくする
    
      plot.Set({{"key",""},{"title","\"divide "+std::to_string(i)+"\""}});
      Network corssNetwork(obj,water);    
      int counter=0;

      // {
      //   std::vector<std::vector<std::vector<double>>> vecvec; 
      //   for(const auto& p:water.getPoints())
      // 	vecvec.push_back({p->xyz, p->getNormal()});
      //   plot.SaveVectorData(vecvec,{{"notitle",""}});      
      // }    

      // {
      //   std::vector<std::vector<std::vector<double>>> vecvec; 
      //   for(const auto& f:water.Faces)
      // 	vecvec.push_back({f->getMeanLocation(), f->getNormal()});
      //   plot.SaveVectorData(vecvec,{{"notitle",""}});      
      // }

      // for(const auto& p:water.getPoints()){
      // 	double x = p->xyz[0];
      // 	double y = p->xyz[1];
      // 	p->xyz += {0,0,cos(M_PI*(x*y)+i/4.)/30.};
      // }
      
      plot.SaveVectorData(getVectorData(water.getFaces()),{{"arrowstyle","1"},{"title","water"}});      
      //      plot.SaveVectorData(getVectorData(obj.getLines()),{{"arrowstyle","2"},{"title","obj"}});          
      //      plot.SaveVectorData(getVectorData(corssNetwork.getLines()),{{"arrowstyle","3"},{"title","corssNetwork"}});
      //      plot.SaveData(corssNetwork.getLocations(),{{"w","p"},{"lc","\"magenta\""},{"pt","7"},{"title","cross points"}});    


      // plot.SaveData(getData(water.Points[50]->getNeighbors_recursive())
      // 		    ,{{"ps","3"},{"w","p"},{"lc","\"red\""},{"pt","7"},{"notitle",""}});    

      // {
      // 	std::vector<std::vector<std::vector<double>>> vecvec; 		
      // 	for(const auto& p:water.Points[0]->getNeighbors_recursive()){
      // 	  for(const auto& l:p->getLines()){
      // 	    pushVectorData(vecvec,l);
      // 	  }
      // 	}
      // 	plot.SaveVectorData(vecvec,{{"lw","2"},{"notitle",""}});
      // }
      
      {
      	std::vector<std::vector<std::vector<double>>> vecvec; 		
      	for(const auto& p:obj.Points[0]->getNeighbors_recursive()){
      	  for(const auto& l:p->getLines()){
      	    pushVectorData(vecvec,l);
      	  }
      	}
      	plot.SaveVectorData(vecvec,{{"lw","2"},{"notitle",""}});
      }

      //      Print(getData(water.Points[50]->getNeighbors_recursive()),Red);

      
      // setStatus(water.getNeighborsIncludeCross(water.getNearestPoints(water.getMeanLocation())),false);
      // plot.SaveData(getLocations(water.getPoints(true)),
      // 		    {{"w","p"},{"lc","\"black\""},{"pt","7"},{"title","cross points"}});
    
      plot.Plot3D_All();
      //      std::cin.ignore();
      plot.Clear();
    }
  };
  ~test1_2(){Print(CHECK+" test1_2",Red);};
};
