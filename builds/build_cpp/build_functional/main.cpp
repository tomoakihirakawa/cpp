
#include <iostream>
#include <functional>
#include <vector>

/* ------------------------------------------------------ */
class human
{
public:
    int age;
    human(const int ageIN) : age(ageIN){};
};
/* ------------------------------------------------------ */
class magic
{
public:
    std::vector<human *> humans;

    magic() : humans(0){};

    void add(human *h)
    {
        this->humans.push_back(h);
    };

    void operator()(std::function<void(human *)> magic_function)
    {
        for (auto &h : humans)
        {
            magic_function(h);
        }
    };
};
/* ------------------------------------------------------ */
int main()
{

    auto magic_function = [](human *h)
    {
        h->age = 0;
    };

    auto taro = new human(10);
    auto jiro = new human(9);

    std::cout << "taro.age , " << taro->age << std::endl;
    std::cout << "jiro.age , " << jiro->age << std::endl;
    magic m;
    m.add(taro);
    m.add(jiro);
    m(magic_function);

    std::cout << "taro.age , " << taro->age << std::endl;
    std::cout << "jiro.age , " << jiro->age << std::endl;
}