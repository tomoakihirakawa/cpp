
//まとめ:どのクラスとして初期化されているのかを意識すべし
//
//クラス宣言文では，オブジェクトのクラスを完全に把握し確定する<-------
//A a1 = a;
//baseクラスのaをコピーしても，a1はAクラスとして初期化される
//
//base a2 = a1;
//Aクラスのa1をコピーしても，a2はbaseクラスとして初期化される
//
//
//ポインターの初期化は，ポインターの指すオブジェクトの初期化に従うので，
//ポインター宣言文では，オブジェクトのクラスを把握できない<-------
//base *a = new A;
//base *aとしても，aが指しているものは，Aクラスである

#include <cstdio>
#include <iostream>

#define test3

#ifdef test1
//-----------------
class base
{
public:
    base(){};
    virtual double operator()() = 0; //オーバーライド必須
};
//-----------------
class A : public base
{
public:
    A(){};
    double operator()() override { return 10; };
};
//-----------------
class B : public base
{
public:
    B(){};
    double operator()() override { return 10; };
};
//-----------------
int main()
{
    A a;
    std::cout << a() << std::endl;

    B b;
    std::cout << b() << std::endl;
};
#endif

#ifdef test2
//baseクラスのポインタとして，`new A`を作成しても，
//作成されたポインタは，クラスAの関数が呼び出してくれる．
//-----------------
class base
{
public:
    base(){};
    virtual double operator()() { return 0; };
};
//-----------------
class A : public base
{
public:
    A(){};
    double operator()() override { return 10; };
};
//-----------------
class B : public base
{
public:
    B(){};
    double operator()() override { return 20; };
};
//-----------------
int main()
{
    base *aa = new base;
    std::cout << (*aa)() << std::endl;

    base *a = new A;
    std::cout << (*a)() << std::endl; //派生クラスAでオーバーライドされた関数が呼び出される

    base *b = new B;
    std::cout << (*b)() << std::endl;
};
#endif

#ifdef test3
//baseクラスのポインタとして，`new A`を作成しても，
//作成されたポインタは，クラスAの関数が呼び出してくれる．
//-----------------
class base
{
public:
    base() { std::cout << "base コンストラクタ" << std::endl; };
    ~base() { std::cout << "base デストラクタ" << std::endl; };
    virtual double operator()() { return 0; };
};
//-----------------
class A : public base
{
public:
    A() { std::cout << "A コンストラクタ" << std::endl; };
    ~A() { std::cout << "A デストラクタ" << std::endl; };
    A(const base &) { std::cout << "A コピーコンストラクタ" << std::endl; };
    double operator()() override { return 10; };
};
//-----------------
class B : public base
{
public:
    B() { std::cout << "B コンストラクタ" << std::endl; };
    ~B() { std::cout << "B デストラクタ" << std::endl; };
    B(const base &) { std::cout << "B コピーコンストラクタ" << std::endl; };
    double operator()() override { return 20; };
};
//-----------------
int main()
{
    base a;
    B b;

    A a1 = a; //オブジェクトa1の初期化時にaをコピーしているので，コピーコンストラクタが呼び出される
    std::cout << a1() << std::endl;
    std::cout << "-------------" << std::endl;
    base a2 = b;
    std::cout << a2() << std::endl;
    std::cout << "-------------" << std::endl;
    b = a;
    std::cout << b() << std::endl;
    std::cout << "-------------" << std::endl;
    //

};
#endif