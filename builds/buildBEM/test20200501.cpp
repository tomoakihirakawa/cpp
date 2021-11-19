
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


  
  return 0;

}

