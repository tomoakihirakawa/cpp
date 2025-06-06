#ifndef lib_measurement_H
#define lib_measurement_H

#include <chrono>
#include <tuple>
#include <vector>

// struct Timer {
//    // std::chrono::time_point<std::chrono::high_resolution_clock> start;
//    //! 一つ目は，前回からの経過時間，二つ目は作成した時間からの経過時間
//    std::chrono::time_point<std::chrono::system_clock> start;
//    std::chrono::time_point<std::chrono::system_clock> lasttime;
//    Timer() : start(std::chrono::system_clock::now()), lasttime(start){};
//    std::array<double, 2> operator()() {
//       auto tmp = this->lasttime;
//       this->lasttime = std::chrono::system_clock::now();
//       // return {this->lasttime - tmp, this->lasttime - this->start};
//       return {(std::chrono::duration<double>(this->lasttime - tmp)).count(),
//               (std::chrono::duration<double>(this->lasttime - this->start)).count()};
//    };
// };

/* ------------------------------------------------------ */

struct TimeWatch {
   // std::chrono::time_point<std::chrono::high_resolution_clock> start;
   std::chrono::time_point<std::chrono::system_clock> start;
   std::chrono::time_point<std::chrono::system_clock> lasttime;
   TimeWatch() {
      this->start = std::chrono::system_clock::now();
      this->lasttime = start;
   };

   void reset() {
      this->start = std::chrono::system_clock::now();
      this->lasttime = start;
   };

   //! 一つ目は，前回からの経過時間，二つ目は作成した時間からの経過時間
   // std::vector<std::chrono::duration<double>> get()
   std::array<double, 2> get() {
      auto tmp = this->lasttime;
      this->lasttime = std::chrono::system_clock::now();
      // return {this->lasttime - tmp, this->lasttime - this->start};
      return {(std::chrono::duration<double>(this->lasttime - tmp)).count(),
              (std::chrono::duration<double>(this->lasttime - this->start)).count()};
   };
   std::array<double, 2> operator()() {
      auto tmp = this->lasttime;
      this->lasttime = std::chrono::system_clock::now();
      // return {this->lasttime - tmp, this->lasttime - this->start};
      return {(std::chrono::duration<double>(this->lasttime - tmp)).count(),
              (std::chrono::duration<double>(this->lasttime - this->start)).count()};
   };
};

#endif