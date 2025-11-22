#include "fundamental.hpp"

int main()
{
    VV_d VV = {{RandomReal({0., 1.}), 0, 0},
               {RandomReal({0., 1.}), RandomReal({0., 1.}), RandomReal({0., 1.})},
               {RandomReal({0., 1.}), RandomReal({0., 1.}), RandomReal({0., 1.})}};
    auto M = Orthogonalize(VV);
    std::cout << M << std::endl;
    for (auto i = 0; i < M.size(); ++i)
        for (auto j = 0; j < M.size(); ++j)
        {
            std::cout << Dot(M[i], M[j]) << std::endl;
        }
};
