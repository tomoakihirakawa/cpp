#include <algorithm>
#include <array>
#include <iostream>
#include <string>
#include <type_traits>
// #include "basic.hpp"
#include "basic_IO.hpp"
#include "basic_arithmetic_array_operations.hpp"
#include "lib_measurement.hpp"

// /Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/build_array/array演算が正しく行われているか.nb

int main() {
   auto rand = 0.71568;
   auto irand = 8;
   std::array<std::array<double, 3>, 3> a{{{0.0800026, 0.0427302, 0.183538}, {0.538665, 0.917936, 0.424192}, {0.512684, 0.0260359, 0.872296}}};
   std::array<std::array<double, 3>, 3> b{{{0.305438, 0.174232, 0.839517}, {0.113538, 0.46904, 0.682904}, {0.659404, 0.104577, 0.269401}}};
   std::array<double, 3> row = {0.0827081, 0.644763, 0.833808};
   std::cout << std::setw(30) << "rand = " << (rand) << std::endl;
   std::cout << std::setw(30) << 0.71568 << std::endl;
   std::cout << std::setw(30) << "irand = " << (irand) << std::endl;
   std::cout << std::setw(30) << 8 << std::endl;
   std::cout << std::setw(30) << "a = " << (a) << std::endl;
   std::cout << std::setw(30) << std::array<std::array<double, 3>, 3>{{{0.0800026, 0.0427302, 0.183538}, {0.538665, 0.917936, 0.424192}, {0.512684, 0.0260359, 0.872296}}} << std::endl;
   std::cout << std::setw(30) << "b = " << (b) << std::endl;
   std::cout << std::setw(30) << std::array<std::array<double, 3>, 3>{{{0.305438, 0.174232, 0.839517}, {0.113538, 0.46904, 0.682904}, {0.659404, 0.104577, 0.269401}}} << std::endl;
   std::cout << std::setw(30) << "a/irand = " << (a / irand) << std::endl;
   std::cout << std::setw(30) << std::array<std::array<double, 3>, 3>{{{0.0100003, 0.00534128, 0.0229422}, {0.0673332, 0.114742, 0.053024}, {0.0640855, 0.00325448, 0.109037}}} << std::endl;
   std::cout << std::setw(30) << "row = " << (row) << std::endl;
   std::cout << std::setw(30) << std::array<double, 3>{{0.0827081, 0.644763, 0.833808}} << std::endl;
   std::cout << std::setw(30) << "row/irand = " << (row / irand) << std::endl;
   std::cout << std::setw(30) << std::array<double, 3>{{0.0103385, 0.0805954, 0.104226}} << std::endl;
   std::cout << "negate" << std::endl;
   std::cout << std::setw(30) << "-a = " << (-a) << std::endl;
   std::cout << std::setw(30) << std::array<std::array<double, 3>, 3>{{{-0.0800026, -0.0427302, -0.183538}, {-0.538665, -0.917936, -0.424192}, {-0.512684, -0.0260359, -0.872296}}} << std::endl;
   std::cout << std::setw(30) << "-b = " << (-b) << std::endl;
   std::cout << std::setw(30) << std::array<std::array<double, 3>, 3>{{{-0.305438, -0.174232, -0.839517}, {-0.113538, -0.46904, -0.682904}, {-0.659404, -0.104577, -0.269401}}} << std::endl;
   std::cout << "+" << std::endl;
   std::cout << std::setw(30) << "a+b = " << (a + b) << std::endl;
   std::cout << std::setw(30) << std::array<std::array<double, 3>, 3>{{{0.38544, 0.216962, 1.02305}, {0.652204, 1.38698, 1.1071}, {1.17209, 0.130612, 1.1417}}} << std::endl;
   std::cout << std::setw(30) << "a-b = " << (a - b) << std::endl;
   std::cout << std::setw(30) << std::array<std::array<double, 3>, 3>{{{-0.225435, -0.131502, -0.655979}, {0.425127, 0.448896, -0.258712}, {-0.14672, -0.0785408, 0.602895}}} << std::endl;
   std::cout << std::setw(30) << "a+row = " << (a + row) << std::endl;
   std::cout << std::setw(30) << std::array<std::array<double, 3>, 3>{{{0.162711, 0.125438, 0.266246}, {1.18343, 1.5627, 1.06896}, {1.34649, 0.859844, 1.7061}}} << std::endl;
   std::cout << std::setw(30) << "row+a = " << (row + a) << std::endl;
   std::cout << std::setw(30) << std::array<std::array<double, 3>, 3>{{{0.162711, 0.125438, 0.266246}, {1.18343, 1.5627, 1.06896}, {1.34649, 0.859844, 1.7061}}} << std::endl;
   std::cout << std::setw(30) << "a+irand = " << (a + irand) << std::endl;
   std::cout << std::setw(30) << std::array<std::array<double, 3>, 3>{{{8.08, 8.04273, 8.18354}, {8.53867, 8.91794, 8.42419}, {8.51268, 8.02604, 8.8723}}} << std::endl;
   std::cout << std::setw(30) << "row+irand = " << (row + irand) << std::endl;
   std::cout << std::setw(30) << std::array<double, 3>{{8.08271, 8.64476, 8.83381}} << std::endl;
   std::cout << "-" << std::endl;
   std::cout << std::setw(30) << "a-row = " << (a - row) << std::endl;
   std::cout << std::setw(30) << std::array<std::array<double, 3>, 3>{{{-0.00270548, -0.0399779, 0.10083}, {-0.106098, 0.273173, -0.220571}, {-0.321124, -0.807772, 0.0384882}}} << std::endl;
   std::cout << std::setw(30) << "row-a = " << (row - a) << std::endl;
   std::cout << std::setw(30) << std::array<std::array<double, 3>, 3>{{{0.00270548, 0.0399779, -0.10083}, {0.106098, -0.273173, 0.220571}, {0.321124, 0.807772, -0.0384882}}} << std::endl;
   std::cout << std::setw(30) << "a+rand = " << (a + rand) << std::endl;
   std::cout << std::setw(30) << std::array<std::array<double, 3>, 3>{{{0.795682, 0.75841, 0.899217}, {1.25435, 1.63362, 1.13987}, {1.22836, 0.741716, 1.58798}}} << std::endl;
   std::cout << std::setw(30) << "rand+a = " << (rand + a) << std::endl;
   std::cout << std::setw(30) << std::array<std::array<double, 3>, 3>{{{0.795682, 0.75841, 0.899217}, {1.25435, 1.63362, 1.13987}, {1.22836, 0.741716, 1.58798}}} << std::endl;
   std::cout << std::setw(30) << "a+irand = " << (a + irand) << std::endl;
   std::cout << std::setw(30) << std::array<std::array<double, 3>, 3>{{{8.08, 8.04273, 8.18354}, {8.53867, 8.91794, 8.42419}, {8.51268, 8.02604, 8.8723}}} << std::endl;
   std::cout << std::setw(30) << "row+irand = " << (row + irand) << std::endl;
   std::cout << std::setw(30) << std::array<double, 3>{{8.08271, 8.64476, 8.83381}} << std::endl;
   std::cout << "*" << std::endl;
   std::cout << std::setw(30) << "row*rand = " << (row * rand) << std::endl;
   std::cout << std::setw(30) << std::array<double, 3>{{0.0591925, 0.461444, 0.596739}} << std::endl;
   std::cout << std::setw(30) << "rand*row = " << (rand * row) << std::endl;
   std::cout << std::setw(30) << std::array<double, 3>{{0.0591925, 0.461444, 0.596739}} << std::endl;
   std::cout << std::setw(30) << "a*rand = " << (a * rand) << std::endl;
   std::cout << std::setw(30) << std::array<std::array<double, 3>, 3>{{{0.0572562, 0.0305811, 0.131354}, {0.385512, 0.656948, 0.303586}, {0.366917, 0.0186333, 0.624285}}} << std::endl;
   std::cout << std::setw(30) << "rand*a = " << (rand * a) << std::endl;
   std::cout << std::setw(30) << std::array<std::array<double, 3>, 3>{{{0.0572562, 0.0305811, 0.131354}, {0.385512, 0.656948, 0.303586}, {0.366917, 0.0186333, 0.624285}}} << std::endl;
   std::cout << std::setw(30) << "a*irand = " << (a * irand) << std::endl;
   std::cout << std::setw(30) << std::array<std::array<double, 3>, 3>{{{0.640021, 0.341842, 1.4683}, {4.30932, 7.34349, 3.39354}, {4.10147, 0.208287, 6.97837}}} << std::endl;
   std::cout << std::setw(30) << "row*irand = " << (row * irand) << std::endl;
   std::cout << std::setw(30) << std::array<double, 3>{{0.661664, 5.1581, 6.67046}} << std::endl;
   std::cout << "/" << std::endl;
   std::cout << std::setw(30) << "row/rand = " << (row / rand) << std::endl;
   std::cout << std::setw(30) << std::array<double, 3>{{0.115566, 0.90091, 1.16506}} << std::endl;
   std::cout << std::setw(30) << "rand/row = " << (rand / row) << std::endl;
   std::cout << std::setw(30) << std::array<double, 3>{{8.65308, 1.10999, 0.858327}} << std::endl;
   std::cout << std::setw(30) << "a/rand = " << (a / rand) << std::endl;
   std::cout << std::setw(30) << std::array<std::array<double, 3>, 3>{{{0.111785, 0.0597058, 0.256452}, {0.752663, 1.28261, 0.592712}, {0.716359, 0.0363792, 1.21884}}} << std::endl;
   std::cout << std::setw(30) << "rand/a = " << (rand / a) << std::endl;
   std::cout << std::setw(30) << std::array<std::array<double, 3>, 3>{{{8.94571, 16.7488, 3.89936}, {1.32862, 0.779662, 1.68716}, {1.39595, 27.4882, 0.820455}}} << std::endl;
   std::cout << std::setw(30) << "a/irand = " << (a / irand) << std::endl;
   std::cout << std::setw(30) << std::array<std::array<double, 3>, 3>{{{0.0100003, 0.00534128, 0.0229422}, {0.0673332, 0.114742, 0.053024}, {0.0640855, 0.00325448, 0.109037}}} << std::endl;
   std::cout << std::setw(30) << "row/irand = " << (row / irand) << std::endl;
   std::cout << std::setw(30) << std::array<double, 3>{{0.0103385, 0.0805954, 0.104226}} << std::endl;
   std::cout << "+=, -=, *=, /=" << std::endl;
   std::cout << std::setw(30) << "a+=b = " << (a += b) << std::endl;
   std::cout << std::setw(30) << std::array<std::array<double, 3>, 3>{{{0.38544, 0.216962, 1.02305}, {0.652204, 1.38698, 1.1071}, {1.17209, 0.130612, 1.1417}}} << std::endl;
   std::cout << std::setw(30) << "a-=b = " << (a -= b) << std::endl;
   std::cout << std::setw(30) << std::array<std::array<double, 3>, 3>{{{0.0800026, 0.0427302, 0.183538}, {0.538665, 0.917936, 0.424192}, {0.512684, 0.0260359, 0.872296}}} << std::endl;
   std::cout << std::setw(30) << "a*=b = " << (a *= b) << std::endl;
   std::cout << std::setw(30) << std::array<std::array<double, 3>, 3>{{{0.0244358, 0.00744498, 0.154083}, {0.0611591, 0.430548, 0.289682}, {0.338066, 0.00272274, 0.234997}}} << std::endl;
   std::cout << std::setw(30) << "a*=row = " << (a *= row) << std::endl;
   std::cout << std::setw(30) << std::array<std::array<double, 3>, 3>{{{0.00202104, 0.00061576, 0.0127439}, {0.0394331, 0.277602, 0.186777}, {0.281882, 0.00227024, 0.195942}}} << std::endl;
   std::cout << std::setw(30) << "a/=b = " << (a /= b) << std::endl;
   std::cout << std::setw(30) << std::array<std::array<double, 3>, 3>{{{0.00661686, 0.00353413, 0.01518}, {0.347312, 0.591851, 0.273503}, {0.42748, 0.0217089, 0.727327}}} << std::endl;
   std::cout << std::setw(30) << "a/=row = " << (a /= row) << std::endl;
   std::cout << std::setw(30) << std::array<std::array<double, 3>, 3>{{{0.0800026, 0.0427302, 0.183538}, {0.538665, 0.917936, 0.424192}, {0.512684, 0.0260359, 0.872296}}} << std::endl;
   std::cout << std::setw(30) << "a/=b = " << (a /= b) << std::endl;
   std::cout << std::setw(30) << std::array<std::array<double, 3>, 3>{{{0.261928, 0.245248, 0.218623}, {4.74435, 1.95705, 0.621159}, {0.777496, 0.248964, 3.23791}}} << std::endl;
   std::cout << std::setw(30) << "a/=row = " << (a /= row) << std::endl;
   std::cout << std::setw(30) << std::array<std::array<double, 3>, 3>{{{3.16689, 2.96523, 2.64331}, {7.35829, 3.03531, 0.963391}, {0.932465, 0.298587, 3.88328}}} << std::endl;
}