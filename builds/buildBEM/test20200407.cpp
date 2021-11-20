
#include "GNUPLOT.hpp"
#include "fundamental.hpp"

using namespace std;

  //   *
  //   |
  // *-*-*
  //   |
  //   *

class Bbasis{
public:
  std::vector<double> xyz;
  double EPS=1E-15;
  int K;
  std::vector<double> q;//knot
  double itvl=1.;//interval  

  std::vector<double> C;
  std::vector<Bbasis*> VecBbasis;
  std::map<Bbasis*,Bbasis*> LinedBbasis;
  std::map<Bbasis*,Bbasis*> BifurcatedBbasis;  
  std::map<Bbasis*,std::vector<double>> Mapq;

  void Bifurcate(Bbasis* a, Bbasis* b){//counter clockwise
    if(VecBbasis.size()==4){
      this->BifurcatedBbasis[a]=b;
      this->BifurcatedBbasis[b]=LinedBbasis[a];
      this->BifurcatedBbasis[LinedBbasis[a]]=LinedBbasis[b];
      this->BifurcatedBbasis[LinedBbasis[b]]=a;      
    }
    else{
      this->BifurcatedBbasis[a]=b;
    }
  };
  
  void Bifurcate(std::vector<Bbasis*> a){
    for(auto i=0; i<a.size()-1; i++)      
      Bifurcate(a[i],a[i+1]);
    
    // if(a.size()==VecBbasis.size()) // a includes all vertex around this node
    //   if(LinedBbasis(*a.begin()) != *a.rbegin())
    // 	Bifurcate(*a.rbegin(),*a.begin());
  };
  
  void Line(Bbasis* a, Bbasis* b){
      this->LinedBbasis[a]=b;
      this->LinedBbasis[b]=a;
  };

  void PointConnect(Bbasis* b_IN){
    for(const auto& b: VecBbasis){
      if(b==b_IN)
        return;
    }
    this->VecBbasis.push_back(b_IN);
  };
  
  Bbasis(const int K_IN): K(K_IN), VecBbasis(0), xyz(3), C(3){
    q.resize(K+2);
    for(auto i=0; i<K+2; i++){
	q[i] = EPS*i - EPS*(K+1)/2.;
      };
  };
  
  double basis0(int k, double xi, const std::vector<double>& q, int i, int iK, int i1, int i1K){
    if(k == 1)
      return (q[i] <= xi && xi < q[iK] ? (xi-q[i])/(q[iK]-q[i]) : 0.) + (q[i1]<=xi && xi<q[i1K] ? (q[i1K]-xi)/(q[i1K]-q[i1]) : 0.);
    return (xi-q[i])/(q[iK]-q[i])*basis0(K-1, xi, q, i, iK-1, i1, i1K-1) + (q[i1K]-xi)/(q[i1K]-q[i1])*basis0(K-1, xi, q, i+1, iK, i1+1, i1K);
  };  
  double basis(Bbasis* a_IN, double xi){
    if(Mapq.count(a_IN)==0)
      return 0.;
    std::vector<double> q = Mapq[a_IN];
    int i=0, iK=K, i1=1, i1K=K+1;
    if(K == 1)
      return (q[i] <= xi && xi < q[iK] ? (xi-q[i])/(q[iK]-q[i]) : 0.) + (q[i1]<=xi && xi<q[i1K] ? (q[i1K]-xi)/(q[i1K]-q[i1]) : 0.);
    return (xi-q[i])/(q[iK]-q[i])*basis0(K-1, xi, q, i, iK-1, i1, i1K-1) + (q[i1K]-xi)/(q[i1K]-q[i1])*basis0(K-1, xi, q, i+1, iK, i1+1, i1K);
  };
  
  // 0 - 1 - 2;
  int checkEdge(Bbasis* a_IN){
    //    -1         0         1         2
    // (*a_IN) -- (*this) --- (*) --- (*edge)
    int ret = -1, maxloop = 10, i=0;
    Bbasis* b = this;
    Bbasis* a = a_IN;
    Bbasis* tmp;
    while(true)
      {
	if(b->VecBbasis.size() <= i){//is accessing point is stored?
	  break;
	}else if(b->VecBbasis[i] == a){//is connected to this poinr?
	  ret++;
	  if(b->LinedBbasis.count(a)==0)//is connected line exist?
	    break;
	  tmp = b;
	  b = b->LinedBbasis[a];
	  a = tmp;
	  i=0;
	}else{
	  i++;
	}
      }
    return ret;
  }
  
  std::vector<double> RetKnots(Bbasis* a_IN){
    std::vector<double> qBK((int)(K/2.+1.),0);
    std::vector<double> qFT((int)(K/2.+1.),0);
    std::vector<double> qCT(K%2,0);
    for(auto i=0; i<qBK.size(); i++){
  	qBK[i] = i + 1. - 0.5 * ((K+1)%2);
  	qFT[i] = i + 1. - 0.5 * ((K+1)%2);      
      }
    std::vector<double> q = Join(Join(-Reverse(qBK),qCT),qFT);    
    //   BK                            FT
    // ------ (*a_IN) ------ (*this) -------> *

    //                             q           q           q       q
    //  nFT                   -1         0           1         2
    //            (*edge) --- (*) --- (*a_IN) --- (*this) --- (*) --- (*edge)    
    int nFT = this->checkEdge(a_IN);

    //  nBK          2         1         0          -1
    //            (*edge) --- (*) --- (*a_IN) --- (*this) --- (*) --- (*edge)    
    int nBK = a_IN->checkEdge(this);
    
    if(nFT==-1 || nBK==-1)
      std::cout << Red << "ERROR: "<< __FUNCTION__ << " @ " << __LINE__ << reset << std::endl;
    
    for(auto& v: qFT)
      if(v >= nFT){
	v = (double)nFT;
      }
    //                             q           q           q       q   
    //  nBK+1        3         2         1           0
    //            (*edge) --- (*) --- (*a_IN) --- (*this) --- (*) --- (*edge)    
    for(auto& v: qBK)
      if(v >= nBK + 1){
	v = (double)nBK + 1;
      }
    
    return Join(Join(-Reverse(qBK),qCT),qFT);
  }

  void StoreKnots(){
    for(const auto b:VecBbasis)
      Mapq[b] = RetKnots(b);
  }

  void SetPosition(const std::vector<double> xyzIN){xyz=xyzIN;};  
  Bbasis* ForwardBbasis(Bbasis* b){
    for(auto tmp:b->VecBbasis){
      if(tmp == b){
	return tmp;
      };
    };
    return this;
  };
};


void PointConnect(Bbasis* b_IN, Bbasis* a_IN){
  b_IN->PointConnect(a_IN);
  a_IN->PointConnect(b_IN);  
};


int main(){

  int K = 2;
   
  std::vector<std::vector<double>> xyz={{0., 3., 0.},
					{1., 3., 1.5},
					{2., 3., 1.},
					{3., 3., 0.},
					{0., 2., 0.},
					{1., 2., 1.},
					{2., 2., 1.},
					{3., 2., 0.},
					{1., 0., 0.},
					{1.5, 1., .5},
					{2., 0., 0.}};

  std::vector<Bbasis*> B;
  
  for(auto i=0; i<11; i++){
    B.push_back(new Bbasis(K));
    (*B.rbegin())->SetPosition(xyz[i]);   
  }

  PointConnect(B[0],B[1]);
  PointConnect(B[1],B[2]);
  PointConnect(B[2],B[3]);
  
  PointConnect(B[4],B[5]);
  PointConnect(B[5],B[6]);
  PointConnect(B[6],B[7]);

  PointConnect(B[8],B[9]);
  PointConnect(B[9],B[10]);

  PointConnect(B[0],B[4]);
  PointConnect(B[1],B[5]);
  PointConnect(B[2],B[6]);
  PointConnect(B[3],B[7]);

  PointConnect(B[4],B[8]);
  PointConnect(B[5],B[9]);
  PointConnect(B[6],B[9]);
  PointConnect(B[7],B[10]);
  
  PointConnect(B[8],B[10]);
  
  // std::vector<std::vector<std::vector<double>>> vec;

  // for(const auto b:B)
  //   {
  //     for(const auto a:b->VecBbasis)
  // 	{
  // 	  vec.push_back({b->xyz, a->xyz - b->xyz});
  // 	}
  //   }
  
  GNUPLOT plot;
  plot.Set({{"hidden3d",""}});
  // plot.VectorPlot3D(vec,{{"lc","1"}});
  // std::cin.ignore();
  
  B[1]->Line(B[0],B[2]);
  B[2]->Line(B[1],B[3]);

  B[5]->Line(B[4],B[6]);
  B[6]->Line(B[5],B[7]);
  
  B[4]->Line(B[0],B[8]);
  B[5]->Line(B[1],B[9]);
  B[6]->Line(B[2],B[9]);
  B[7]->Line(B[3],B[10]);

  B[9]->Line(B[5],B[10]);
  B[9]->Line(B[6],B[8]);    
  B[8]->Line(B[4],B[10]);
  B[10]->Line(B[8],B[7]);



  B[6]->Bifurcate(B[7],B[2]);
  B[5]->Bifurcate(B[6],B[1]);
  B[9]->Bifurcate(B[6],B[5]);
  B[8]->Bifurcate(B[9],B[4]);

  B[0]->Bifurcate(B[4],B[1]);
  B[1]->Bifurcate({B[0],B[5],B[2]});
  B[2]->Bifurcate({B[1],B[6],B[3]});
  B[3]->Bifurcate({B[2],B[7]});
  B[4]->Bifurcate({B[8],B[5],B[0]});
  B[8]->Bifurcate({B[10],B[9],B[4]});
  B[10]->Bifurcate({B[7],B[9],B[8]});
  B[7]->Bifurcate({B[3],B[6],B[10]});
  
  // vec.clear();
  // plot.Clear();
  // for(const auto b:B)
  //   {
  //     for(const auto a: b->LinedBbasis)
  // 	vec.push_back({b->xyz, a.second->xyz - b->xyz});
  //   }
  // plot.VectorPlot3D(vec,{{"lc","1"}});
  // std::cin.ignore(); 
  
  cout << B[0]->checkEdge(B[1]) << endl;
  cout << B[0]->RetKnots(B[1]) << endl;

  cout << B[0]->checkEdge(B[4]) << endl;
  cout << B[0]->RetKnots(B[4]) << endl;

  cout << B[1]->checkEdge(B[2]) << endl;
  cout << B[1]->RetKnots(B[2]) << endl;

  cout << B[4]->checkEdge(B[10]) << endl;
  cout << B[4]->RetKnots(B[10]) << endl;

  for(auto b:B)
    b->StoreKnots();

  // {
  
  // std::vector<std::vector<int>> IND{
  // 				    {0,1,4},
  // 				    {6,2,7},
  // 				    {6,5,2},
  // 				    {6,7,9},
  // 				    {4,5,0},
  // 				    {4,5,8},
  // 				    {5,4,1},
  // 				    {5,4,9},
  // 				    {5,6,1},
  // 				    {1,5,0},
  // 				    {1,5,2},
  // 				    {2,1,6},
  // 				    {2,6,3},
  // 				    {3,2,7},
  // 				    {7,3,6},
  // 				    {7,10,6},
  // 				    {9,8,10},
  // 				    {9,6,10},
  // 				    {9,5,8},				    
  // 				    {8,9,4},
  // 				    {8,9,10},
  // 				    {10,7,9},
  // 				    {10,9,8}};
  // for(auto ind:IND)
  //   {
  //     Bbasis* fromB0 = B[ind[1]];
  //     Bbasis* fromB1 = B[ind[2]];
  //     Bbasis* toB = B[ind[0]];
  //     {
  // 	std::vector<std::vector<std::vector<double>>> VEC;
  // 	for(auto i=0; i<26; i++){
  // 	  std::vector<std::vector<double>> vec;
  // 	  for(auto j=0; j<26; j++){
  // 	    double xi=-i/25., eta=-j/25.;
  // 	    std::vector<double> xyz = ((toB->xyz - fromB0->xyz)*xi + (toB->xyz - fromB1->xyz)*eta) + toB->xyz;
  // 	    std::vector<double> plus = {0.,0.,(toB->basis(fromB0,xi))*(toB->basis(fromB1,eta))};
  // 	    vec.push_back( std::vector<double>{xyz[0],xyz[1],0.} + plus);
  // 	  }
  // 	  VEC.push_back(vec);
  // 	}
  // 	plot.SaveSplotData(VEC,{{"w","l"},{"lc","2"},{"hidden3d",""}});
  //     }
  //   }    
  // }
  

  for(auto toB:B)      
  {
    for(auto v:toB->VecBbasis)      
      {
	Bbasis* fromB0 = v;
	if(toB->BifurcatedBbasis.count(v)!=0)
	  {
	    Bbasis* fromB1 = toB->BifurcatedBbasis[v];	
	    std::vector<std::vector<std::vector<double>>> VEC;
	      for(auto i=0; i<21; i++){
		std::vector<std::vector<double>> vec;
		for(auto j=0; j<21; j++){
		  double xi=-i/12., eta=-j/12.;
		  std::vector<double> xyz = ((toB->xyz - fromB0->xyz)*xi + (toB->xyz - fromB1->xyz)*eta) + toB->xyz;
		  std::vector<double> plus = {0.,0.,(toB->basis(fromB0,xi))*(toB->basis(fromB1,eta))};
		  vec.push_back( std::vector<double>{xyz[0],xyz[1],0.} + plus);
		}
		VEC.push_back(vec);
	      }
	      plot.SaveSplotData(VEC,{{"w","l"},{"hidden3d",""}});
	    }
	}
    }
    
  plot.Splot();
  std::cin.ignore();
  
  // vec.clear();
  // for(auto l=-100; l<101; l++)
  //   vec.push_back({l/50.,B[1]->basis(B[0],l/50.)});  
  // plot.SaveData(vec,{{"w","l"}});  

  // vec.clear();  
  // for(auto l=-100; l<101; l++)
  //   vec.push_back({1.-l/50.,B[1]->basis(B[2],l/50.)});  
  // plot.SaveData(vec,{{"w","l"}});


  // vec.clear();  
  // for(auto l=-100; l<101; l++)
  //   vec.push_back({2.-l/50.,B[2]->basis(B[3],l/50.)});  
  // plot.SaveData(vec,{{"w","l"}});
  
  // plot.Plot2D();  
  // std::cin.ignore();

  // // cout << b.q << endl;

  // // std::vector<std::vector<std::vector<double>>> vector3{
  // // 							{{1.,1.,1.},{1.,2.,2.}},
  // // 							{{1.,1.,1.},{1.,2.,2.}}
  // // };
  
  // // GNUPLOT plot;
  // // plot.Plot3D({xyz},{{"lc","1"}});

  // std::cin.ignore();


  // for(auto i=0; i<10; i++)
  //   {
  //     B[i]->Setq();
  //     cout << "i: " << B[i]->checkEdge(2) << endl;
  //     cout << B[i]->q << endl;;
  //   }
  
  return 0;

}

