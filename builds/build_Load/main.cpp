#include "fundamental.hpp"

int main()
{

  std::string file = "/Users/tomoaki/Downloads/GEBCO_24_Jan_2022_f73708a40215/gebco_2021_n40.64308941364288_s39.06588077545166_w139.46481049060822_e140.11575043201447.asc";
  using VVV_d = std::vector<std::vector<std::vector<double>>>;

  std::vector<V_s> read_line;
  Load(file, read_line, {"    ", "   ", "  ", " "});
  std::cout << stod(read_line) << std::endl;
}
