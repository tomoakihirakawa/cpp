# !/bin/zh

g++-10 -std=c++17 -Ofast -fopenmp -x c++-header -o fundamental_exception.hpp.gch fundamental_exception.hpp
g++-10 -std=c++17 -Ofast -fopenmp -x c++-header -o basic_arithmetic_vector_operations.hpp.gch basic_arithmetic_vector_operations.hpp
g++-10 -std=c++17 -Ofast -fopenmp -x c++-header -o fundamental_constants.hpp.gch fundamental_constants.hpp
g++-10 -std=c++17 -Ofast -fopenmp -x c++-header -o fundamental_vectors.hpp.gch fundamental_vectors.hpp
g++-10 -std=c++17 -Ofast -Winvalid-pch -fopenmp -H -x c++-header -o fundamental.hpp.gch fundamental.hpp
# g++-10 -std=c++17 -Ofast -fopenmp -x c++-header -o fundamental_geometry.hpp.gch fundamental_geometry.hpp 
# g++-10 -std=c++17 -Ofast -fopenmp -x c++-header -o Network.hpp.gch Network.hpp
g++-10 -std=c++17 -Ofast -fopenmp -x c++-header -o GNUPLOT.hpp.gch GNUPLOT.hpp
