class test1_7{
  using V_double = std::vector<double>;
  using VV_double = std::vector<std::vector<double>>;
  using VVV_double = std::vector<std::vector<std::vector<double>>>;
public:
  
  template <class T>
  bool AnyStatusTrue(const std::vector<T*>& obj){
    for(const auto& o:obj)
      if(o->getStatus())
	return true;
    return false;
  };

  template <class T>
  bool AllStatusTrue(const std::vector<T*>& obj){
    for(const auto& o:obj)
      if(!o->getStatus())
	return false;
    return true;
  };
  
  double length(const std::vector<std::vector<double>>& vv){  
    return Norm(vv[0]-vv[1]);
  };
  double length(const networkLine* l){  
    return length(l->getLocations());
  };
  // SECURE DIVIDING CHECK
  test1_7(){
    Print("**************** test1_7 ****************", Red);
    NetworkObj obj("./obj/tank.obj");
    NetworkW water({1,1},{1.133,1.13,1/2.},.1);
    GNUPLOT plot;   
    //============================
    for(auto i=0; i<16; i++){
      plot.Set({{"style","arrow 1 nohead lc \"blue\" lw .5"}});
      plot.Set({{"style","arrow 2 nohead lc \"magenta\" lw .5 dt 2"}});    
      plot.Set({{"style","arrow 3 nohead lc \"red\" lw .5 dt 2"}});
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

      if(false){//干渉している箇所をプロット
	Network corssNetwork(obj,water);	
	auto save = [this, &plot](Network& obj){
		      for(const auto& f:network::takeIfIntersect(obj.Faces)){
			for(const auto& ps:f->getPointsCutFaces()){
	    
			  Print(ps,Green);
			  auto s = ps.size();
	    
			  for(const auto& routeP:f->getPointsCutFaces()){

			    if(AnyStatusTrue(routeP)){
			      VVV_double VVV;
			      for(const auto& index: triangulate(routeP,f->getNormal()) ){
				VVV.push_back({routeP[index[0]]->getX(),routeP[index[1]]->getX()-routeP[index[0]]->getX()});
				VVV.push_back({routeP[index[1]]->getX(),routeP[index[2]]->getX()-routeP[index[1]]->getX()});
				VVV.push_back({routeP[index[2]]->getX(),routeP[index[0]]->getX()-routeP[index[2]]->getX()});		
			      }
			      plot.SaveVectorData(VVV,{{"arrowstyle","2"},{"notitle",""}});
			    }
			  }
			}	  
		      }
		    };
      
	network::setStatus(obj.Points,true);
	save(obj);

	network::setStatus(water.Points,true);
	save(water);
      }
      
      if(false){//obj.Points[0]から繋がる干渉している箇所をプロット
	Network corssNetwork(obj,water);      
	auto save = [this, &plot](Network& obj){
		      network::setStatus(obj.Points[0]->getNeighbors_AvoidCross_recursive(),true);		      
		      for(const auto& f:network::takeIfIntersect(obj.Faces)){
			for(const auto& ps:f->getPointsCutFaces()){
	    
			  Print(ps,Green);
			  auto s = ps.size();
	    
			  for(const auto& routeP:f->getPointsCutFaces()){
			    if(AnyStatusTrue(routeP)){
			      VVV_double VVV;
			      for(const auto& index: triangulate(routeP,f->getNormal()) ){
				VVV.push_back({routeP[index[0]]->getX(),routeP[index[1]]->getX()-routeP[index[0]]->getX()});
				VVV.push_back({routeP[index[1]]->getX(),routeP[index[2]]->getX()-routeP[index[1]]->getX()});
				VVV.push_back({routeP[index[2]]->getX(),routeP[index[0]]->getX()-routeP[index[2]]->getX()});		
			      }
			      plot.SaveVectorData(VVV,{{"arrowstyle","2"},{"notitle",""}});
			    }
			  }
			}	  
		      }
		    };
      

	save(obj);

	save(water);
      }

      if(false){//制限の流れを確認
	Network corssNetwork(obj,water);      
	auto save = [this, &plot](Network& obj){
		      plot.Set({{"style","arrow 1 nohead lc \"blue\" lw 1 dt 2"}});
		      plot.Set({{"style","arrow 2 nohead lc \"red\" lw 1 dt 1"}});
		      plot.Set({{"style","arrow 3 nohead lc \"magenta\" lw 1 dt 2"}});
		      plot.Set({{"style","arrow 4 nohead lc \"orange\" lw 1 dt 1"}});      		      		      
		      //-----------------
		      VVV_double VVV0;
		      auto faces = obj3D::takeInsideOfBounds(obj.Faces,{{-0.55,0.55},{-0.55,0.55},{0.05,1.5}});
		      pushVectorData(VVV0,faces);
		      plot.SaveVectorData(VVV0,{{"arrowstyle","1"},{"title","takeInsideOfBounds"}});
		      plot.plot3d();// std::cin.ignore();
		      plot.Clear();
		      //------------------
		      plot.SaveVectorData(VVV0,{{"arrowstyle","1"},{"title","takeInsideOfBounds"}});		      
		      VVV_double VVV1;
		      auto face = network::takeNearest(faces,{0,0,0.2});
		      pushVectorData(VVV1,face);
		      plot.SaveVectorData(VVV1,{{"arrowstyle","2"},{"title","takeNearest"}});
		      plot.plot3d();// std::cin.ignore();
		      plot.Clear();
		      //------------------
		      plot.SaveVectorData(VVV0,{{"arrowstyle","1"},{"title","takeInsideOfBounds"}});
		      plot.SaveVectorData(VVV1,{{"arrowstyle","2"},{"title","takeNearest"}});
		      VVV_double VVV2;
		      pushVectorData(VVV2,network::takeIfStatus(face->getNeighbors_AvoidIntersection_recursive(),true));
		      plot.SaveVectorData(VVV2,{{"arrowstyle","3"},{"title","getNeighbors_AvoidIntersection_recursive"}});
		      plot.plot3d();// std::cin.ignore();
		      plot.Clear();
		      //------------------
		      auto startP = network::takeNearest(obj3D::takeInsideOfBounds(obj.Points,{{-0.55,0.55},{-0.55,0.55},{0.1,1.}}),{0,0,0});
		      auto ps = startP->getNeighbors_AvoidCross_recursive();
		      network::setStatus(ps,true);
		      plot.SaveVectorData(VVV0,{{"arrowstyle","1"},{"title","takeInsideOfBounds"}});
		      plot.SaveVectorData(VVV1,{{"arrowstyle","2"},{"title","takeNearest"}});
		      plot.SaveVectorData(VVV2,{{"arrowstyle","3"},{"title","getNeighbors_AvoidIntersection_recursive"}});
		      plot.SaveData(obj3D::extractX((ps)),{{"title","true"},{"pt","7"},{"ps","2"},{"lc","\"pink\""}});
		      plot.SaveData({startP->getX()},{{"title","true"},{"pt","7"},{"ps","2"},{"lc","\"red\""}});		      
		      plot.plot3d();// std::cin.ignore();
		      plot.Clear();
		      //------------------
		      VVV_double VVV;
		      network::setStatus(obj.Points,false);		      
		      network::setStatus(ps,true);
		      for(const auto& f:network::takeIfIntersect(obj.Faces)){
			for(const auto& routeP:f->getPointsCutFaces()){
			  if(AllStatusTrue(routeP)){
			    for(const auto& ind: triangulate(routeP,f->getNormal()) ){
			      VVV.push_back({routeP[ind[0]]->getX(),routeP[ind[1]]->getX()-routeP[ind[0]]->getX()});
			      VVV.push_back({routeP[ind[1]]->getX(),routeP[ind[2]]->getX()-routeP[ind[1]]->getX()});
			      VVV.push_back({routeP[ind[2]]->getX(),routeP[ind[0]]->getX()-routeP[ind[2]]->getX()});		
			    }
			  }
			}
		      }
		      plot.SaveVectorData(VVV0,{{"arrowstyle","1"},{"title","takeInsideOfBounds"}});
		      plot.SaveVectorData(VVV1,{{"arrowstyle","2"},{"title","takeNearest"}});
		      plot.SaveVectorData(VVV2,{{"arrowstyle","3"},{"title","getNeighbors_AvoidIntersection_recursive"}});
		      plot.SaveData(obj3D::extractX(ps),{{"title","true"},{"pt","7"},{"ps","2"},{"lc","\"pink\""}});
		      plot.SaveData({startP->getX()},{{"title","true"},{"pt","7"},{"ps","2"},{"lc","\"red\""}});
		      plot.SaveVectorData(VVV,{{"arrowstyle","4"},{"title",""}});
		      plot.plot3d();// std::cin.ignore();
		      plot.Clear();
		      
		    };
      
	save(obj);
	save(water);
      }

      if(false){//制限の流れを確認
	Network corssNetwork(obj,water);
	network::setStatus(corssNetwork.Points,false);	
	auto save = [this, &plot](Network& obj, Network& obj2){
		      //plot.Set({{"xrange","[0.1:0.55]"},{"yrange","[0.1:-0.55]"},{"zrange","[0.:1.]"}});		      
		      plot.Set({{"style","arrow 1 nohead lc \"blue\" lw 1 dt 2"}});
		      plot.Set({{"style","arrow 2 nohead lc \"red\" lw 1 dt 1"}});
		      plot.Set({{"style","arrow 3 nohead lc \"magenta\" lw 3 dt 2"}});
		      plot.Set({{"style","arrow 4 nohead lc \"orange\" lw 1 dt 1"}});
		      plot.Set({{"style","arrow 5 nohead lc \"web-blue\" lw 1 dt 1"}});      		      		      		      
		      //-----------------
		      VVV_double VVV0;
		      auto faces = obj3D::takeInsideOfBounds(obj.Faces,{{-0.55,0.55},{-0.55,0.55},{0.05,1.5}});
		      pushVectorData(VVV0,faces);
		      plot.SaveVectorData(VVV0,{{"arrowstyle","1"},{"title","takeInsideOfBounds"}});

		      VVV_double YYY0;
		      auto faces1 = obj3D::takeInsideOfBounds(obj2.Faces,{{-0.55,0.55},{-0.55,0.55},{0.05,1.5}});
		      pushVectorData(YYY0,faces1);
		      plot.SaveVectorData(YYY0,{{"arrowstyle","1"},{"title","takeInsideOfBounds"}});

		      plot.plot3d();std::cin.ignore();plot.Clear();
		      //------------------
		      plot.SaveVectorData(VVV0,{{"arrowstyle","1"},{"title","takeInsideOfBounds"}});		      
		      plot.SaveVectorData(YYY0,{{"arrowstyle","1"},{"title","takeInsideOfBounds"}});		      
		      VVV_double VVV1;
		      auto face = network::takeNearest(faces,{0,0,0.2});
		      pushVectorData(VVV1,face);
		      plot.SaveVectorData(VVV1,{{"arrowstyle","2"},{"title","takeNearest"}});

		      VVV_double YYY1;
		      auto face1 = network::takeNearest(faces1,{0,0,0.2});
		      pushVectorData(YYY1,face1);
		      plot.SaveVectorData(YYY1,{{"arrowstyle","2"},{"title","takeNearest"}});

		      plot.plot3d();std::cin.ignore();plot.Clear();
		      //------------------
		      plot.SaveVectorData(VVV0,{{"arrowstyle","1"},{"title","takeInsideOfBounds"}});
		      plot.SaveVectorData(VVV1,{{"arrowstyle","2"},{"title","takeNearest"}});
		      plot.SaveVectorData(YYY0,{{"arrowstyle","1"},{"title","takeInsideOfBounds"}});
		      plot.SaveVectorData(YYY1,{{"arrowstyle","2"},{"title","takeNearest"}});

		      VVV_double VVV2;
		      pushVectorData(VVV2,network::takeIfStatus(face->getNeighbors_AvoidIntersection_recursive(),true));
		      plot.SaveVectorData(VVV2,{{"arrowstyle","3"},{"title","getNeighbors_AvoidIntersection_recursive"}});

		      VVV_double YYY2;
		      pushVectorData(YYY2,network::takeIfStatus(face1->getNeighbors_AvoidIntersection_recursive(),true));
		      plot.SaveVectorData(YYY2,{{"arrowstyle","3"},{"title","getNeighbors_AvoidIntersection_recursive"}});

		      plot.plot3d();std::cin.ignore();plot.Clear();
		      //------------------		      
		      auto startP = network::takeNearest(obj3D::takeInsideOfBounds(obj.Points,{{-0.55,0.55},{-0.55,0.55},{0.1,1.}}),{0,0,0});
		      auto ps = startP->getNeighbors_AvoidCross_recursive();
		      network::setStatus(ps,true);
		      plot.SaveVectorData(VVV0,{{"arrowstyle","1"},{"title","takeInsideOfBounds"}});
		      plot.SaveVectorData(VVV1,{{"arrowstyle","2"},{"title","takeNearest"}});
		      plot.SaveVectorData(VVV2,{{"arrowstyle","3"},{"title","getNeighbors_AvoidIntersection_recursive"}});
		      plot.SaveData(obj3D::extractX((ps)),{{"title","true"},{"pt","7"},{"ps","2"},{"lc","\"pink\""}});
		      plot.SaveData({startP->getX()},{{"title","true"},{"pt","7"},{"ps","2"},{"lc","\"red\""}});		      					      
		      
		      auto startP2 = network::takeNearest(obj3D::takeInsideOfBounds(obj2.Points,{{-0.55,0.55},{-0.55,0.55},{0.1,1.}}),{0,0,0});
		      auto ps2 = startP2->getNeighbors_AvoidCross_recursive();
		      network::setStatus(ps2,true);
		      plot.SaveVectorData(YYY0,{{"arrowstyle","1"},{"title","takeInsideOfBounds"}});
		      plot.SaveVectorData(YYY1,{{"arrowstyle","2"},{"title","takeNearest"}});
		      plot.SaveVectorData(YYY2,{{"arrowstyle","3"},{"title","getNeighbors_AvoidIntersection_recursive"}});
		      plot.SaveData(obj3D::extractX((ps2)),{{"title","true"},{"pt","7"},{"ps","2"},{"lc","\"green\""}});
		      plot.SaveData({startP2->getX()},{{"title","true"},{"pt","7"},{"ps","2"},{"lc","\"sea-green\""}});		      			
		      
		      plot.plot3d();std::cin.ignore();plot.Clear();
		      //------------------
		      VVV_double VVV, VVV_route;
		      // network::setStatus(obj.Points,false); 
		      // network::setStatus(ps,true);

		      startP = network::takeNearest(obj3D::takeInsideOfBounds(obj.Points,{{-0.55,0.55},{-0.55,0.55},{0.1,1.}}),{0,0,0});
		      network::setStatus(startP->getNeighbors_AvoidCross_recursive(),true);

		      for(const auto& f:network::takeIfIntersect(obj.Faces)){
			for(const auto& routeP:f->getPointsCutFacesBehind()){			  
			  plot.SaveData(obj3D::extractX(routeP),{{"w","lp"},{"notitle",""},{"lw","3"},{"lc","1"},{"pt","7"},{"loop",""}});	
			}
		      }
		      
		      // for(const auto& f:network::takeIfIntersect(obj.Faces)){
		      // 	for(const auto& routeP:f->getPointsCutFaces()){			  
		      // 	  if(AllStatusTrue(routeP)){
		      // 	    for(const auto& ind: triangulate(routeP,f->getNormal()) ){
		      // 	      VVV.push_back({routeP[ind[0]]->getX(),routeP[ind[1]]->getX()-routeP[ind[0]]->getX()});
		      // 	      VVV.push_back({routeP[ind[1]]->getX(),routeP[ind[2]]->getX()-routeP[ind[1]]->getX()});
		      // 	      VVV.push_back({routeP[ind[2]]->getX(),routeP[ind[0]]->getX()-routeP[ind[2]]->getX()});		
		      // 	    }
		      // 	  }
		      // 	}
		      // }

		      VVV_double YYY;
		      network::setStatus(obj.Points,false);		      
		      network::setStatus(ps,true);
		      for(const auto& f:network::takeIfIntersect(obj.Faces)){
			for(const auto& routeP:f->getPointsCutFacesBehind()){
			  if(true){
			    for(const auto& ind: triangulate(routeP,f->getNormal()) ){
			      YYY.push_back({routeP[ind[0]]->getX(),routeP[ind[1]]->getX()-routeP[ind[0]]->getX()});
			      YYY.push_back({routeP[ind[1]]->getX(),routeP[ind[2]]->getX()-routeP[ind[1]]->getX()});
			      YYY.push_back({routeP[ind[2]]->getX(),routeP[ind[0]]->getX()-routeP[ind[2]]->getX()});		
			    }
			  }
			}
		      }


		      // int counter=0;
		      // for(const auto& f:obj.Faces){
		      // 	for(const auto& xyz:obj3D::extractX(f->getPointsCutFaces()))
		      // 	  plot.SaveData(xyz,{{"w","lp"},{"notitle",""},{"lc",std::to_string(counter)},{"pt",std::to_string(counter)}});	
		      // 	counter++;								
		      // }
		      
		      plot.SaveVectorData(VVV0,{{"arrowstyle","1"},{"title","takeInsideOfBounds"}});
		      plot.SaveVectorData(VVV1,{{"arrowstyle","2"},{"title","takeNearest"}});
		      plot.SaveVectorData(VVV2,{{"arrowstyle","3"},{"title","getNeighbors_AvoidIntersection_recursive"}});
		      plot.SaveData(obj3D::extractX(ps),{{"title","true"},{"pt","7"},{"ps","2"},{"lc","\"pink\""}});
		      plot.SaveData({startP->getX()},{{"title","true"},{"pt","7"},{"ps","2"},{"lc","\"red\""}});
		      plot.SaveVectorData(VVV,{{"arrowstyle","4"},{"title","VVV"}});

		      plot.SaveVectorData(YYY0,{{"arrowstyle","1"},{"title","takeInsideOfBounds"}});
		      plot.SaveVectorData(YYY1,{{"arrowstyle","2"},{"title","takeNearest"}});
		      plot.SaveVectorData(YYY2,{{"arrowstyle","3"},{"title","getNeighbors_AvoidIntersection_recursive"}});
		      plot.SaveData(obj3D::extractX(ps2),{{"title","true"},{"pt","7"},{"ps","2"},{"lc","\"green\""}});
		      plot.SaveData({startP2->getX()},{{"title","true"},{"pt","7"},{"ps","2"},{"lc","\"sea-green\""}});
		      plot.SaveVectorData(YYY,{{"arrowstyle","5"},{"title","YYY"}});

		      plot.plot3d();std::cin.ignore();plot.Clear();
		      
		    };
	save(obj,water);
      }

      if(false){//制限の流れを確認
	Network corssNetwork(obj,water);
	network::setStatus(corssNetwork.Points,false);	
	auto save = [this, &corssNetwork, &plot](Network& obj, Network& obj2){
		      //plot.Set({{"xrange","[0.1:0.55]"},{"yrange","[0.1:-0.55]"},{"zrange","[0.:1.]"}});		      
		      plot.Set({{"style","arrow 1 nohead lc \"blue\" lw 1 dt 2"}});
		      plot.Set({{"style","arrow 2 nohead lc \"red\" lw 1 dt 1"}});
		      plot.Set({{"style","arrow 3 nohead lc \"magenta\" lw 3 dt 2"}});
		      plot.Set({{"style","arrow 4 nohead lc \"orange\" lw 1 dt 1"}});
		      plot.Set({{"style","arrow 5 nohead lc \"web-blue\" lw 1 dt 1"}});      		      		      		      
		      //-----------------
		      VVV_double VVV0;
		      auto faces = obj3D::takeInsideOfBounds(obj.Faces,{{-0.55,0.55},{-0.55,0.55},{0.05,1.5}});
		      pushVectorData(VVV0,faces);
		      plot.SaveVectorData(VVV0,{{"arrowstyle","1"},{"title","takeInsideOfBounds"}});

		      VVV_double YYY0;
		      auto faces1 = obj3D::takeInsideOfBounds(obj2.Faces,{{-0.6,0.6},{-0.6,0.6},{0.05,1.5}});
		      pushVectorData(YYY0,faces1);
		      plot.SaveVectorData(YYY0,{{"arrowstyle","1"},{"title","takeInsideOfBounds"}});

		      plot.plot3d();/*std::cin.ignore();*/plot.Clear();
		      //------------------
		      plot.SaveVectorData(VVV0,{{"arrowstyle","1"},{"title","takeInsideOfBounds"}});		      
		      plot.SaveVectorData(YYY0,{{"arrowstyle","1"},{"title","takeInsideOfBounds"}});		      
		      VVV_double VVV1;
		      auto face = network::takeNearest(faces,{0,0,0.2});
		      pushVectorData(VVV1,face);
		      plot.SaveVectorData(VVV1,{{"arrowstyle","2"},{"title","takeNearest"}});

		      VVV_double YYY1;
		      auto face1 = network::takeNearest(faces1,{0,0,0.2});
		      pushVectorData(YYY1,face1);
		      plot.SaveVectorData(YYY1,{{"arrowstyle","2"},{"title","takeNearest"}});

		      plot.plot3d();/*std::cin.ignore();*/plot.Clear();
		      //------------------
		      plot.SaveVectorData(VVV0,{{"arrowstyle","1"},{"title","takeInsideOfBounds"}});
		      plot.SaveVectorData(VVV1,{{"arrowstyle","2"},{"title","takeNearest"}});
		      plot.SaveVectorData(YYY0,{{"arrowstyle","1"},{"title","takeInsideOfBounds"}});
		      plot.SaveVectorData(YYY1,{{"arrowstyle","2"},{"title","takeNearest"}});

		      VVV_double VVV2;
		      pushVectorData(VVV2,network::takeIfStatus(face->getNeighbors_AvoidIntersection_recursive(),true));
		      plot.SaveVectorData(VVV2,{{"arrowstyle","3"},{"title","getNeighbors_AvoidIntersection_recursive"}});

		      VVV_double YYY2;
		      pushVectorData(YYY2,network::takeIfStatus(face1->getNeighbors_AvoidIntersection_recursive(),true));
		      plot.SaveVectorData(YYY2,{{"arrowstyle","3"},{"title","getNeighbors_AvoidIntersection_recursive"}});

		      plot.plot3d();/*std::cin.ignore();*/plot.Clear();
		      //------------------		      
		      auto startP = network::takeNearest(obj3D::takeInsideOfBounds(obj.Points,{{-0.55,0.55},{-0.55,0.55},{0.1,1.}}),{0,0,0});
		      auto ps = startP->getNeighbors_AvoidCross_recursive();
		      network::setStatus(ps,true);
		      plot.SaveVectorData(VVV0,{{"arrowstyle","1"},{"title","takeInsideOfBounds"}});
		      plot.SaveVectorData(VVV1,{{"arrowstyle","2"},{"title","takeNearest"}});
		      plot.SaveVectorData(VVV2,{{"arrowstyle","3"},{"title","getNeighbors_AvoidIntersection_recursive"}});
		      plot.SaveData(obj3D::extractX((ps)),{{"title","true"},{"pt","7"},{"ps","2"},{"lc","\"pink\""}});
		      plot.SaveData({startP->getX()},{{"title","true"},{"pt","7"},{"ps","2"},{"lc","\"red\""}});		      					      
		      
		      auto startP2 = network::takeNearest(obj3D::takeInsideOfBounds(obj2.Points,{{-0.55,0.55},{-0.55,0.55},{0.1,1.}}),{0,0,0});
		      auto ps2 = startP2->getNeighbors_AvoidCross_recursive();
		      network::setStatus(ps2,true);
		      plot.SaveVectorData(YYY0,{{"arrowstyle","1"},{"title","takeInsideOfBounds"}});
		      plot.SaveVectorData(YYY1,{{"arrowstyle","2"},{"title","takeNearest"}});
		      plot.SaveVectorData(YYY2,{{"arrowstyle","3"},{"title","getNeighbors_AvoidIntersection_recursive"}});
		      plot.SaveData(obj3D::extractX((ps2)),{{"title","true"},{"pt","7"},{"ps","2"},{"lc","\"green\""}});
		      plot.SaveData({startP2->getX()},{{"title","true"},{"pt","7"},{"ps","2"},{"lc","\"sea-green\""}});		      			
		      
		      plot.plot3d();/*std::cin.ignore();*/plot.Clear();
		      //------------------
		      VVV_double VVV, VVV_route;
		      // network::setStatus(obj.Points,false); 
		      // network::setStatus(ps,true);

		      startP = network::takeNearest(obj3D::takeInsideOfBounds(obj.Points,{{-0.55,0.55},{-0.55,0.55},{0.1,1.}}),{0,0,0});
		      network::setStatus(startP->getNeighbors_AvoidCross_recursive(),true);

		      for(const auto& f:network::takeIfIntersect(obj.Faces)){
			for(const auto& routeP:f->getPointsCutFacesBehind()){			  
			  plot.SaveData(obj3D::extractX(routeP),{{"w","lp"},{"notitle",""},{"lw","3"},{"lc","1"},{"pt","7"},{"loop",""}});	
			}
		      }
		      
		      // for(const auto& f:network::takeIfIntersect(obj.Faces)){
		      // 	for(const auto& routeP:f->getPointsCutFaces()){			  
		      // 	  if(AllStatusTrue(routeP)){
		      // 	    for(const auto& ind: triangulate(routeP,f->getNormal()) ){
		      // 	      VVV.push_back({routeP[ind[0]]->getX(),routeP[ind[1]]->getX()-routeP[ind[0]]->getX()});
		      // 	      VVV.push_back({routeP[ind[1]]->getX(),routeP[ind[2]]->getX()-routeP[ind[1]]->getX()});
		      // 	      VVV.push_back({routeP[ind[2]]->getX(),routeP[ind[0]]->getX()-routeP[ind[2]]->getX()});		
		      // 	    }
		      // 	  }
		      // 	}
		      // }

		      VVV_double YYY;
		      network::setStatus(obj.Points,false);		      
		      network::setStatus(ps,true);
		      for(const auto& f:network::takeIfIntersect(obj.Faces)){
			for(const auto& routeP:f->getPointsCutFacesBehind()){
			  if(true){
			    for(const auto& ind: triangulate(routeP,f->getNormal()) ){
			      YYY.push_back({routeP[ind[0]]->getX(),routeP[ind[1]]->getX()-routeP[ind[0]]->getX()});
			      YYY.push_back({routeP[ind[1]]->getX(),routeP[ind[2]]->getX()-routeP[ind[1]]->getX()});
			      YYY.push_back({routeP[ind[2]]->getX(),routeP[ind[0]]->getX()-routeP[ind[2]]->getX()});		
			    }
			  }
			}
		      }


		      // int counter=0;
		      // for(const auto& f:obj.Faces){
		      // 	for(const auto& xyz:obj3D::extractX(f->getPointsCutFaces()))
		      // 	  plot.SaveData(xyz,{{"w","lp"},{"notitle",""},{"lc",std::to_string(counter)},{"pt",std::to_string(counter)}});	
		      // 	counter++;								
		      // }
		      
		      plot.SaveVectorData(VVV0,{{"arrowstyle","1"},{"title","takeInsideOfBounds"}});
		      plot.SaveVectorData(VVV1,{{"arrowstyle","2"},{"title","takeNearest"}});
		      plot.SaveVectorData(VVV2,{{"arrowstyle","3"},{"title","getNeighbors_AvoidIntersection_recursive"}});
		      plot.SaveData(obj3D::extractX(ps),{{"title","true"},{"pt","7"},{"ps","2"},{"lc","\"pink\""}});
		      plot.SaveData({startP->getX()},{{"title","true"},{"pt","7"},{"ps","2"},{"lc","\"red\""}});
		      plot.SaveVectorData(VVV,{{"arrowstyle","4"},{"title","VVV"}});

		      plot.SaveVectorData(YYY0,{{"arrowstyle","1"},{"title","takeInsideOfBounds"}});
		      plot.SaveVectorData(YYY1,{{"arrowstyle","2"},{"title","takeNearest"}});
		      plot.SaveVectorData(YYY2,{{"arrowstyle","3"},{"title","getNeighbors_AvoidIntersection_recursive"}});
		      plot.SaveData(obj3D::extractX(ps2),{{"title","true"},{"pt","7"},{"ps","2"},{"lc","\"green\""}});
		      plot.SaveData({startP2->getX()},{{"title","true"},{"pt","7"},{"ps","2"},{"lc","\"sea-green\""}});
		      plot.SaveVectorData(YYY,{{"arrowstyle","5"},{"title","YYY"}});

		      plot.plot3d();/*std::cin.ignore();*/plot.Clear();
		      
		    };
	save(obj,water);
      }

      if(false){//制限の流れを確認
	Network corssNetwork(obj,water);
	auto save = [this, &plot](Network& obj, Network& obj2){
		      
		      for(const auto& f:obj.Faces)
			for(const auto& xyz:obj3D::extractX(f->getPointsOnLinesDivided()))
			  plot.SaveData(xyz,{{"w","l"},{"notitle",""}});					

		      for(const auto& f:obj2.Faces)
			for(const auto& xyz:obj3D::extractX(f->getPointsOnLinesDivided()))
			  plot.SaveData(xyz,{{"w","l"},{"notitle",""}});		      
		      
		      plot.plot3d();/*std::cin.ignore();*/plot.Clear();
		      
		    };
	save(obj,water);
      }

      if(true){//制限の流れを確認
	auto save = [this, &plot](Network& obj, Network& obj2){
		      Network corssNetwork(obj,obj2);	
		      //plot.Set({{"xrange","[0.1:0.55]"},{"yrange","[0.1:-0.55]"},{"zrange","[0.:1.]"}});
		      plot.Set({{"zrange","[0.:]"}});		      		      
		      plot.Set({{"style","arrow 1 nohead lc \"blue\" lw 1 dt 2"}});
		      plot.Set({{"style","arrow 2 nohead lc \"red\" lw 1 dt 2"}});
		      plot.Set({{"style","arrow 3 nohead lc \"magenta\" lw 1 dt 1"}});
		      plot.Set({{"style","arrow 4 nohead lc \"web-blue\" lw 1 dt 3"}});
		      plot.Set({{"style","arrow 5 nohead lc \"orange\" lw 1 dt 1"}});      		      		      		      

		      V_netFp faces;
		      V_netPp points;		      
		      VV_double bounds = {{-0.55,0.55},{-0.55,0.55},{0.05,1.5}};

		      // FACES BASED
		      //===================
		      
		      // getNeighbors_AvoidIntersection_recursive()
		      //---------------------------------------------
		      // faces = network::takeNearest(obj3D::takeInsideOfBounds(obj.Faces, bounds), {0,0,0.2})->getNeighbors_AvoidIntersection_recursive();		      
		      // plot.SaveVectorData(getVectorData(faces),{{"arrowstyle","1"},{"title","takeInsideOfBounds"}});
		      // faces = network::takeNearest(obj3D::takeInsideOfBounds(obj2.Faces, bounds), {0,0,0.2})->getNeighbors_AvoidIntersection_recursive();
		      // plot.SaveVectorData(getVectorData(faces),{{"arrowstyle","2"},{"title","takeInsideOfBounds"}});


		      // getNeighbors_UpToIntersection_recursive()
		      //---------------------------------------------
		      // faces = network::takeNearest(obj3D::takeInsideOfBounds(obj.Faces, bounds), {0,0,0.2})->getNeighbors_UpToIntersection_recursive();		      
		      // plot.SaveVectorData(getVectorData(faces),{{"arrowstyle","1"},{"title","takeInsideOfBounds"}});		      
		      // faces = network::takeNearest(obj3D::takeInsideOfBounds(obj2.Faces, bounds), {0,0,0.2})->getNeighbors_UpToIntersection_recursive();
		      // plot.SaveVectorData(getVectorData(faces),{{"arrowstyle","2"},{"title","takeInsideOfBounds"}});
		      
		      // getNeighbors_UpToIntersection_recursive(checkedL,passedL)
		      //---------------------------------------------
		      // V_netLp passedL;
		      // V_netLp checkedL;		      
		      // network::takeNearest(obj3D::takeInsideOfBounds(obj.Faces, bounds),{0,0,0})->getNeighbors_UpToIntersection_recursive(checkedL,passedL);
		      // plot.SaveVectorData(getVectorData(passedL),{{"lc","'red'"},{"title","passedL"}});
		      // plot.SaveVectorData(getVectorData(checkedL),{{"lc","'blue'"},{"title","checkedL"}});

		      
		      // One-Way
		      //-----------------
		      // faces = obj2.Faces[100]->getNeighbors_OneWay_recursive();
		      // plot.SaveVectorData(getVectorData(faces),{{"arrowstyle","1"},{"title","takeInsideOfBounds"}});		      

		      // faces = obj.Faces[0]->getNeighbors_OneWay_recursive();
		      // plot.SaveVectorData(getVectorData(faces),{{"arrowstyle","2"},{"title","takeInsideOfBounds"}});		      


		      // POINTS BASED
		      //===================

		      // getNeighbors_AvoidCross_recursive()
		      //---------------------------------------------
		      // points = network::takeNearest(obj3D::takeInsideOfBounds(obj.Points,bounds),{0,0,0})->getNeighbors_AvoidCross_recursive();
		      // plot.SaveData(getData(points),{{"title","takeInsideOfBounds"}});
		      // points = network::takeNearest(obj3D::takeInsideOfBounds(obj2.Points,bounds),{0,0,0})->getNeighbors_AvoidCross_recursive();
		      // plot.SaveData(getData(points),{{"title","takeInsideOfBounds"}});

		      
		      // getNeighbors_UpToIntersection_recursive(checkedL,passedL)
		      //---------------------------------------------
		      // V_netLp passedL;
		      // V_netLp checkedL;		      
		      // network::takeNearest(obj3D::takeInsideOfBounds(obj.Points, bounds),{0,0,0})->getNeighbors_UpToIntersection_recursive(checkedL,passedL);
		      // plot.Set({{"pm3d",""}});
		      // plot.SaveVectorData(getVectorData(passedL),{{"lc","'red'"},{"title","passedL"}});
		      // plot.SaveVectorData(getVectorData(checkedL),{{"lc","'blue'"},{"title","checkedL"}});
		      // plot.SaveVectorData(getVectorData(extractFaces(checkedL)),{{"lc","'dark-green'"},{"title","extractFaces"}});
		      // plot.SaveVectorData(getVectorData(obj2.Faces),{{"lc","'web-blue'"},{"title","water"}});

		      
		      // One-Way
		      //-----------------
		      // points = obj2.Points[10]->getNeighbors_OneWay_recursive();
		      // plot.SaveData(getData(points),{{"w","lp"},{"title","takeInsideOfBounds"}});		      

		      // points = obj.Points[20]->getNeighbors_OneWay_recursive();
		      // plot.SaveData(getData(points),{{"w","lp"},{"title","takeInsideOfBounds"}});		      

		      // points = corssNetwork.Points[0]->getNeighbors_OneWay_recursive();
		      // plot.SaveData(getData(points),{{"w","lp"},{"title","takeInsideOfBounds"}});		      
		      
		      // INTERSECTION
		      //===================
		      // // true
		      // faces = network::takeIfIntersect(obj.Faces);
		      // plot.SaveVectorData(getVectorData(faces),{{"arrowstyle","3"},{"title","takeInsideOfBounds"}});		      
		      // faces = network::takeIfIntersect(obj2.Faces);
		      // plot.SaveVectorData(getVectorData(faces),{{"arrowstyle","3"},{"title","takeInsideOfBounds"}});
		      // // false
		      // faces = network::takeIfIntersect(obj.Faces,false);
		      // plot.SaveVectorData(getVectorData(faces),{{"arrowstyle","4"},{"title","takeInsideOfBounds"}});		      
		      // faces = network::takeIfIntersect(obj2.Faces,false);
		      // plot.SaveVectorData(getVectorData(faces),{{"arrowstyle","4"},{"title","takeInsideOfBounds"}});

		      
		      // ACCESS TO FACES FROM POINTS
		      // ===================
		      // netP* p; 		      
		      // p = network::takeNearest(obj3D::takeInsideOfBounds(obj.Points,bounds),{0,0,0});
		      // plot.SaveVectorData(getVectorData(p->getFaces()),{{"arrowstyle","3"},{"title","takeInsideOfBounds"}});

		      // p = network::takeNearest(obj3D::takeInsideOfBounds(obj2.Points,bounds),{0,0,0});		      
		      // plot.SaveVectorData(getVectorData(p->getFaces()),{{"arrowstyle","3"},{"title","takeInsideOfBounds"}});

		      // ACCESS TO POINTS FROM FACES
		      // ===================
		      // netF* f; 		      
		      // f = network::takeNearest(obj3D::takeInsideOfBounds(obj.Faces,bounds),{0,0,0});
		      // plot.SaveData(getData(f->getPoints()),{{"loop",""},{"w","lp"},{"title","takeInsideOfBounds"}});

		      // f = network::takeNearest(obj3D::takeInsideOfBounds(obj2.Faces,bounds),{0,0,0});		      
		      // plot.SaveData(getData(f->getPoints()),{{"loop",""},{"w","lp"},{"title","takeInsideOfBounds"}});

		      
		      // BISECTIONING
		      // ===================
		      // netF* f;
		      // for(const auto& l:obj3D::takeInsideOfBounds(obj.getLines(),{{-0.55,0.55},{-0.55,0.55},{-100,100}}))
		      // 	obj.divide(longerLine(l));
		      // plot.SaveVectorData(getVectorData(obj.Faces),{{"arrowstyle","1"},{"title","takeInsideOfBounds"}});


		      // GROUPING INTERSECTION ROUTES
		      // ===================
		      // for(const auto& f:network::takeIfIntersect(obj.Faces)){
		      // 	for(const auto& routeP:f->getPointsCutFacesBehind()){
		      // 	  for(const auto& ind: triangulate(routeP,f->getNormal())){
		      // 	    plot.SaveVectorData({{routeP[ind[0]]->getX(),routeP[ind[1]]->getX()-routeP[ind[0]]->getX()},
		      // 				 {routeP[ind[1]]->getX(),routeP[ind[2]]->getX()-routeP[ind[1]]->getX()},
		      // 				 {routeP[ind[2]]->getX(),routeP[ind[0]]->getX()-routeP[ind[2]]->getX()}},
		      // 	      {{"arrowstyle","4"},{"notitle",""}});		
		      // 	  }
		      // 	}
		      // }

		      // for(const auto& f:network::takeIfIntersect(obj2.Faces)){
		      // 	for(const auto& routeP:f->getPointsCutFacesBehind()){
		      // 	  for(const auto& ind: triangulate(routeP,f->getNormal())){
		      // 	    plot.SaveVectorData({{routeP[ind[0]]->getX(),routeP[ind[1]]->getX()-routeP[ind[0]]->getX()},
		      // 				 {routeP[ind[1]]->getX(),routeP[ind[2]]->getX()-routeP[ind[1]]->getX()},
		      // 				 {routeP[ind[2]]->getX(),routeP[ind[0]]->getX()-routeP[ind[2]]->getX()}},
		      // 	      {{"arrowstyle","5"},{"notitle",""}});		
		      // 	  }
		      // 	}
		      // }
		      
		      // GET CUT FACES TRIANGULATED
		      // ===================
		      // {
			
		      // 	for(const auto& f:network::takeIfIntersect(obj.Faces)){
		      // 	  for(const auto& routeP:f->getPointsCutFacesBehind()){
		      // 	    plot.SaveData(getData(routeP),{{"loop",""},{"ps",".5"},{"w","lp"},{"lc","'magenta'"},{"lw",".5"},{"notitle",""}});			    
		      // 	    for(const auto& ind: triangulate(routeP,f->getNormal())){
		      // 	      int s = ind.size();
		      // 	      for(auto i=0; i<s; i++)
		      // 	    	plot.SaveVectorData({{routeP[ind[i]]->getX(),routeP[ind[(i+1)%s]]->getX()-routeP[ind[i]]->getX()}},{{"lw",".1"},{"lc","'magenta'"},{"notitle",""}});				
		      // 	    }
		      // 	  }
		      // 	}
			
		      // 	for(const auto& f:network::takeIfIntersect(obj2.Faces)){
		      // 	  for(const auto& routeP:f->getPointsCutFacesBehind()){
		      // 	    plot.SaveData(getData(routeP),{{"loop",""},{"ps",".5"},{"w","lp"},{"lc","'blue'"},{"lw",".5"},{"notitle",""}});			    
		      // 	    for(const auto& ind: triangulate(routeP,f->getNormal())){
		      // 	      int s = ind.size();
		      // 	      for(auto i=0; i<s; i++)
		      // 	    	plot.SaveVectorData({{routeP[ind[i]]->getX(),routeP[ind[(i+1)%s]]->getX()-routeP[ind[i]]->getX()}},{{"lw",".1"},{"lc","'blue'"},{"notitle",""}});
		      // 	    }
		      // 	  }
		      // 	}
			
		      // }

		      // GET CUT FACES TRIANGULATED BACKS OF PENETRATING FACES
		      // ===================
		      // {
			
		      // 	for(const auto& f:network::takeIfIntersect(obj.Faces)){
		      // 	  for(const auto& routeP:f->getPointsCutFaces()){
		      // 	    plot.SaveData(getData(routeP),{{"loop",""},{"ps",".5"},{"w","lp"},{"lc","'magenta'"},{"lw",".5"},{"notitle",""}});			    
		      // 	    for(const auto& ind: triangulate(routeP,f->getNormal())){
		      // 	      int s = ind.size();
		      // 	      for(auto i=0; i<s; i++)
		      // 	    	plot.SaveVectorData({{routeP[ind[i]]->getX(),routeP[ind[(i+1)%s]]->getX()-routeP[ind[i]]->getX()}},{{"lw",".1"},{"lc","'magenta'"},{"notitle",""}});				
		      // 	    }
		      // 	  }
		      // 	}
			
		      // for(const auto& f:network::takeIfIntersect(obj2.Faces)){
		      //   for(const auto& routeP:f->getPointsCutFaces()){
		      //     plot.SaveData(getData(routeP),{{"loop",""},{"ps",".5"},{"w","lp"},{"lc","'blue'"},{"lw",".5"},{"notitle",""}});			    
		      //     for(const auto& ind: triangulate(routeP,f->getNormal())){
		      //       int s = ind.size();
		      //       for(auto i=0; i<s; i++)
		      //     	plot.SaveVectorData({{routeP[ind[i]]->getX(),routeP[ind[(i+1)%s]]->getX()-routeP[ind[i]]->getX()}},{{"lw",".1"},{"lc","'blue'"},{"notitle",""}});
		      //     }
		      //   }
		      // }
			
		      // }

		      // GET FACES BEHIND PENETRATING FACES
		      // ===================
		      // {
		      // 	int i=0;
		      // 	for(const auto& p:obj.Points){
		      // 	  if(i++!=0) break;			  

			  
		      // 	  auto fs = p->getFaces();

			  
		      // 	  plot.SaveVectorData(getVectorData(fs),{{"lc",std::to_string(i)},
		      // 						 {"lw",std::to_string(i%2+.5)},
		      // 						 {"title",std::to_string(i)}});


		      // 	  for(const auto& f:network::takeIfIntersect(fs)){
		      // 	    for(const auto& routeP:f->getPointsCutFaces()){
		      // 	      plot.SaveData(getData(routeP),{{"loop",""},{"ps",".5"},{"w","lp"},{"lc","'blue'"},{"lw",".5"},{"notitle",""}});			    
		      // 	      for(const auto& ind: triangulate(routeP,f->getNormal())){
		      // 		int s = ind.size();
		      // 		for(auto i=0; i<s; i++)
		      // 		  plot.SaveVectorData({{routeP[ind[i]]->getX(),routeP[ind[(i+1)%s]]->getX()-routeP[ind[i]]->getX()}},{{"lw",".1"},{"lc","'blue'"},{"notitle",""}});
		      // 	      }
		      // 	    }
		      // 	  }			  
		      // 	}			
		      // }

		      // CALCULATE IG IGn USING BEM 
		      //=====================
		      // using map_P_Vd = std::map<netP*, V_double>;
		      // using map_F_P_Vd = std::map<netF*, map_P_Vd>;
		      // using map_P_P_Vd = std::map<netP*, map_P_Vd>;		      
		      // using map_P_F_P_Vd = std::map<netP*, map_F_P_Vd>;
		      		      
		      // /// CHECK OBJECTS		      
		      // V_netLp passedL;
		      // V_netLp checkedL;
		      // network::takeNearest(obj3D::takeInsideOfBounds(obj.Points, bounds),{0,0,0})->getNeighbors_UpToIntersection_recursive(checkedL,passedL);
		      // auto Faces = extractFaces(checkedL);
		      // auto Points = extractPoints(checkedL);
		      // GNUPLOT plot_;
		      // plot_.Set({{"key",""}});
		      // plot_.SaveVectorData(getVectorData(Faces),{{"lc","'red'"},{"title","Faces"}});
		      // plot_.SaveData(getData(Points),{{"title","Points"}});
		      // plot_.SaveVectorData(getVectorData(passedL),{{"lc","'red'"},{"title","passedL"}});
		      // plot_.SaveVectorData(getVectorData(checkedL),{{"lc","'blue'"},{"title","checkedL"}});
		      // plot_.plot3d();
		      // std::cin.ignore();
		     

		      // /// CALCULATING IG IGn
		      // VV_double gw = GaussianQuadratureWeights(5,0.,1.);
		      // map_P_P_Vd P_P_IGIGn;		      
		      // map_P_Vd P_IGIGn;
		      // map_P_Vd init_P_IGIGn;		      
		      
		      // for(const auto& p:Points)
		      //   init_P_IGIGn[p] = {0.,0.};		      

		      // for(const auto& a:Points){			
		      // 	P_IGIGn = init_P_IGIGn;
		      // 	for(const auto& f:Faces)
		      // 	  for(const auto& [p, igign]:BEM::IGIGn_(f,a,gw))
		      // 	    P_IGIGn[p] += igign;
						
		      //   P_P_IGIGn[a] = P_IGIGn;
		      // }
		      
		      // VV_double IG, IGn;
		      // for(const auto& tmp:P_P_IGIGn){
		      // 	V_double ig, ign;		      
		      // 	for(const auto& igign:tmp.second)
		      // 	  {
		      // 	    ig.emplace_back(igign.second[0]);
		      // 	    ign.emplace_back(igign.second[1]);
		      // 	  }
		      // 	IG.emplace_back(ig);
		      // 	IGn.emplace_back(ign);
		      // }

		      // GNUPLOT mat;
		      // mat.SaveMatrixData(IG,{{"w","image"}});
		      // mat.MatrixPlot(); 
		      // std::cin.ignore();
		      
		      // GNUPLOT mat2;
		      // mat2.SaveMatrixData(IGn,{{"w","image"}});
		      // mat2.MatrixPlot(); 
		      // std::cin.ignore();

		      // CALCULATE IG IGn USING BEM 
		      //=====================
		      using map_P_Vd = std::map<netP*, V_double>;
		      using map_F_P_Vd = std::map<netF*, map_P_Vd>;
		      using map_P_P_Vd = std::map<netP*, map_P_Vd>;		      
		      using map_P_F_P_Vd = std::map<netP*, map_F_P_Vd>;

		      //------
		      V_netFp Faces;
		      V_netPp Points;		      
		      GNUPLOT plot_;
		      plot_.Set({{"key",""}});
		      for(const auto& OBJ:{obj, obj2})
		      {
			V_netLp passedL;
			V_netLp checkedL;
			network::takeNearest(obj3D::takeInsideOfBounds(OBJ.Points, bounds),{0,0,0})->getNeighbors_UpToIntersection_recursive(checkedL,passedL);
			auto faces = extractFaces(checkedL);
			auto points = extractPoints(checkedL);
			plot_.SaveVectorData(getVectorData(faces),{{"lc","'red'"},{"title","Faces"}});
			plot_.SaveData(getData(points),{{"title","Points"}});
			plot_.SaveVectorData(getVectorData(passedL),{{"lc","'red'"},{"title","passedL"}});
			plot_.SaveVectorData(getVectorData(checkedL),{{"lc","'blue'"},{"title","checkedL"}});

			Faces.insert(Faces.begin(),faces.begin(),faces.end());
			Points.insert(Points.begin(),points.begin(),points.end());			
		      }
		      plot_.plot3d();
		      std::cin.ignore();

		      
		      /// CALCULATING IG IGn
		      VV_double gw = GaussianQuadratureWeights(5,0.,1.);
		      map_P_P_Vd P_P_IGIGn;		      
		      map_P_Vd P_IGIGn;
		      map_P_Vd init_P_IGIGn;		      
		      
		      for(const auto& p:Points)
		        init_P_IGIGn[p] = {0.,0.};		      

		      for(const auto& a:Points){			
		      	P_IGIGn = init_P_IGIGn;
		      	for(const auto& f:Faces)
		      	  for(const auto& [p, igign]:BEM::IGIGn_(f,a,gw))
		      	    P_IGIGn[p] += igign;
						
		        P_P_IGIGn[a] = P_IGIGn;
		      }
		      
		      
		      VV_double matOfKnowns, matOfUnknowns;
		      for(const auto& [origin, p_igign]:P_P_IGIGn){
		      	V_double ig, ign;
		      	for(const auto& igign:p_igign)
		      	  {
		      	    ig.emplace_back(igign.second[0]);
		      	    ign.emplace_back(igign.second[1]);
		      	  }

			if(origin->getNetwork() == &obj){
			  Print("true",Green);
			  matOfKnowns.emplace_back(ig);
			  matOfUnknowns.emplace_back(ign);
			}else if(origin->getNetwork() == &obj2){
			  Print("false",Blue);			  
			  matOfKnowns.emplace_back(ign);
			  matOfUnknowns.emplace_back(ig);
			}else{
			  Print("ERROR",Red);
			}
			
		      }

		      
		      GNUPLOT mat;
		      mat.SaveMatrixData(Inverse(matOfKnowns),{{"w","image"}});
		      mat.MatrixPlot(); 
		      std::cin.ignore();
		      
		      GNUPLOT mat2;
		      mat2.SaveMatrixData(matOfUnknowns,{{"w","image"}});
		      mat2.MatrixPlot(); 
		      std::cin.ignore();

		      /////////////////////////////////////
		      
		      // map_P_Vd P_IGIGn;		      
		      // for(const auto& [pp, f_p_igign]:P_F_P_IGIGn)
		      // 	for(const auto& [f, p_igign]:f_p_igign)			
		      // 	  for(const auto& [p, igign]:p_igign)
		      // 	    P_IGIGn[p] = {0.,0.};

		      		      
		      // for(const auto& pfpV:P_F_P_IGIGn)
		      // 	for(const auto& fpV:pfpV.second)			
		      // 	  for(const auto& [p, igign]:fpV.second){
		      // 	    if(!std::isnan(igign[0])){
		      // 	      P_IGIGn[p] += igign;
		      // 	      Print(igign,red);			    
		      // 	      Print(P_IGIGn[p],green);
		      // 	    }
		      // 	  }

		      // for(const auto& [p,igign]:P_IGIGn){
		      // 	Print(p,red);
		      // 	Print(igign,green);
		      // }
		      
		      
		    };
	
	save(obj,water);
      }
      
      plot.plot3d();
      std::cin.ignore();
      plot.Clear();      
    }
  };
  ~test1_7(){Print(CHECK+" test1_7",Red);};
};
