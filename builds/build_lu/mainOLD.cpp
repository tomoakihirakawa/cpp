
#include "fundamental_vectors.hpp"
#include <iostream>
#include <chrono>

struct myludcmp{
  int n;
  VV_d lu;
  std::vector<int> indx;
  double d;
  VV_d aref;

  myludcmp(const VV_d &a) : n(a.size()), lu(a), aref(a), indx(n) {
    const double TINY=1.0e-40;
    int i, imax, j, k;
    double big,temp;
    V_d vv(n);
    d=1.0;
    
    for(const auto& row:lu){
      big=0.0;
      for(const auto& e:row)
        if((temp=std::abs(e)) > big)
          big=temp;     
      if(big == 0.0)
        throw("Singular matrix in LUdcmp");
      vv[i]=1.0/big;
    }

    //row
    for (k=0;k<n;k++) {
      big=0.0;
      for (i=k;i<n;i++) {
	temp=vv[i]*std::abs(lu[i][k]);
	if (temp > big) {
	  big=temp;
	  imax=i;
	}
      }
      
      if (k != imax) {
        std::swap(lu[k],lu[imax]);
	d = - d;
	vv[imax]=vv[k];
      }
      
      indx[k]=imax;
      
      if (lu[k][k] == 0.0)
        lu[k][k]=TINY;

      for (i=k+1;i<n;i++) {
	temp=lu[i][k] /= lu[k][k];
	for (j=k+1;j<n;j++)
	  lu[i][j] -= temp*lu[k][j];
      }
    }
    
  };
  
  void solve(const V_d &b, V_d &x)
  {
    int i,ii=0,ip,j;
    double sum;
    if (b.size() != n || x.size() != n)
      throw("solve bad sizes");

    x=b;
    
    for (i=0;i<n;i++) {
      ip=indx[i];
      sum=x[ip];
      x[ip]=x[i];

      if (ii != 0)
	for (j=ii-1;j<i;j++)
          sum -= lu[i][j]*x[j];
      else if (sum != 0.0)
	ii=i+1;
      x[i]=sum;
    }
    for (i=n-1;i>=0;i--) {
      sum=x[i];
      for (j=i+1;j<n;j++)
        sum -= lu[i][j]*x[j];
      x[i]=sum/lu[i][i];
    }
    
  };

  void solve(const VV_d &b,VV_d &x)
  {
    int i,j,m=b[0].size();
    if (b.size() != n || x.size() != n || b[0].size() != x.size())
      throw("solve bad sizes");
    V_d xx(n);
    for (j=0;j<m;j++) {
      for (i=0;i<n;i++)
        xx[i] = b[i][j];
      solve(xx,xx);
      for (i=0;i<n;i++)
        x[i][j] = xx[i];
    }
  };
  
  void inverse(VV_d &ainv){
    int i=0,j=0;
    ainv.resize(n,V_d(n,0));
    for (i=0;i<n;i++) {
      for (j=0;j<n;j++) ainv[i][j] = 0.;
      ainv[i][i] = 1.;
    }
    solve(ainv,ainv);
  };

  VV_d inverse(){    
    VV_d ainv(n,V_d(n,0));
    int i,j;
    for (i=0;i<n;i++) {
      for (j=0;j<n;j++) ainv[i][j] = 0.;
      ainv[i][i] = 1.;
    }
    solve(ainv,ainv);
    return ainv;
  };

  double det()
  {
    double dd = d;
    for (int i=0;i<n;i++)
      dd *= lu[i][i];
    return dd;
  };
  
  void mprove(V_d &b, V_d &x)
  {
    int i,j;
    V_d r(n);
    for (i=0;i<n;i++) {
      double sdp = -b[i];
      for (j=0;j<n;j++)
	sdp += (double)aref[i][j] * (double)x[j];
      r[i]=sdp;
    }
    solve(r,r);
    for (i=0;i<n;i++) x[i] -= r[i];
  };
};


int main(){

  using V_d = std::vector<double>;
  using VV_d = std::vector<std::vector<double>>;

  VV_d w = {{1.,1.,2.,3.},
            {4.,4.,9.,3.},
            {4.,5.,3.,3.},
            {1.,3.,1.,3.}};

  auto start = std::chrono::high_resolution_clock::now();

  myludcmp lu(w);    
  std::cout << lu.inverse() << std::endl;  

  auto finish = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed = finish - start;

  std::cout << "Elapsed time: " << elapsed.count() << " s\n";
  
}
