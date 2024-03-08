#include <array>
#include "Network.hpp"
#include "lib_multipole_expansion.hpp"

/*DOC_EXTRACT 1_1_Bessel_function

## ベッセル関数



*/

int main() {

   auto net = new Network("./bunny.obj", "bunny");
   net->makeBucketFaces(net->getScale() / 10.);
   net->makeBucketPoints(net->getScale() / 10.);

   /*



   */

   auto bucket_double = copyPartition<double>(net->BucketPoints);
   auto bucket_tdd = copyPartition<Tddd>(net->BucketPoints);

   for (const auto& cell : net->BucketPoints) {
      auto g = gradG(X, A);
      auto error = std::log10(std::abs(1 - Gapx(n, X, A, center) / G(X, A)));
      auto v = Norm(g - gradGapx(n, X, A, center)) / Norm(g);
   }
}