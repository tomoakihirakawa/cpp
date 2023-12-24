#include <iostream>
// #include <csdlib>

using namespace std;

class array{
  int *p;
  int size;
public:
  array(int sz){
    p=new int[sz];
    if(!p) exit(1);
    this->size = sz;
    cout << "use default constructor" << endl;
  }
  ~array(){delete [] p;};
  void put(int i, int j){
    if(i>=0 && i<size) p[i] = j;    
  }

  array(const array &a);
  
  int get(int i){
    return p[i];
  }
};
/*初期化の際だけ呼び出される．単なる代入では，利用されないので注意*/
array::array(const array &a){
  int i;
  
  size=a.size;
  p=new int[a.size];
  if(!p) exit(1);
  for(i=0; i<a.size; i++)
    p[i]=a.p[i];
  cout << "use copy constructor" << endl;
}

int main(){

  array num(10);
  cout << "num=" << &num << endl;
  
  int i;

  for(i=0; i<10; i++)
    num.put(i,i);

  for(i=9; i>=0; i--) cout << num.get(i) << endl;
  cout << endl;

  /*よく紹介される，コピーコンストラクタを使った方法*/
  array x = num/*これをコピー*/;
  cout << "x=" << &x << endl;//numとは違うメモリに保存されていることがわかる
  
  for(i=0; i<10; i++) cout << x.get(i) << endl;
  cout << endl;

  /*これでもコピーコンストラクタが呼び出される*/
  array* y = new array(x/*これをコピー*/);
  cout << "y=" << y << endl;//xとは違うメモリに保存されていることがわかる
  
  for(i=0; i<10; i++) cout << y->get(i) << endl;
  cout << endl;

  /*以下は，デフォルトコンストラクタが呼び出される*/
  array* z = new array(10);
  cout << "z=" << z << endl;
  
  for(i=0; i<10; i++) cout << z->get(i) << endl;
  cout << endl;

}
