#ifndef basic_exception_H
#define basic_exception_H

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include "basic_constants.hpp"

class error_message : public std::exception {
  public:
   std::string filename;
   std::string name;
   int line;
   std::string message;

   error_message(const std::string &filename_,
                 const std::string &name_,
                 const int line_,
                 const std::string &message_)
       : filename(filename_), name(name_), line(line_), message(message_) {
      std::cout << nani() << std::endl;
   };

   std::string nani() const {
      std::ostringstream ss;
      ss << colorOff << "----------------------------------------------------\n"
         << colorOff << "    FILE: " << red << filename << ":" << Green << line << "\n"
         << colorOff << "FUNCTION: " << Magenta << name << "\n"
         << colorOff << " message: " << Red << message << colorOff << "\n"
         << colorOff << "----------------------------------------------------\n";
      return ss.str();
   };

   void print() const {
      std::cout << nani() << std::endl;
   };
};

std::string message(const std::string &filename_,
                    const std::string &name_,
                    const int line_,
                    const std::string &message_) {
   std::ostringstream ss;
   ss << colorOff << Green << " LINE: " << line_ << colorOff << Magenta << " FUNCTION: " << name_ << colorOff << blue << " FILE: " << filename_ << colorOff << red << " : " << message_ << colorOff;
   return ss.str();
};

#endif
