#ifndef basic_exception_H
#define basic_exception_H

#include <iostream>
#include <string>
#include <stdexcept>
#include "basic_constants.hpp"

class error_message : public std::exception
{
public:
  std::string filename;
  std::string name;
  int line;
  std::string message;
  error_message(const std::string &filename_,
                const std::string &name_,
                const int line_,
                const std::string &message_) : filename(filename_),
                                               name(name_),
                                               line(line_),
                                               message(message_)
  {
    std::cout << nani() << std::endl;
  };

  std::string nani() const
  {
    return colorOff + "----------------------------------------------------\n" +
           colorOff + "    LINE: " + Green + std::to_string(this->line) + "\n" +
           colorOff + "FUNCTION: " + Magenta + this->name + "\n" +
           colorOff + "    FILE: " + red + this->filename + "\n" +
           colorOff + " message: " + Red + this->message + colorOff + "\n" +
           colorOff + "----------------------------------------------------\n";
  };

  void print() const
  {
    std::cout << nani() << std::endl;
  };
};

std::string message(const std::string &filename_,
                    const std::string &name_,
                    const int line_,
                    const std::string &message_)
{
  return colorOff + green + " LINE: " + std::to_string(line_) + colorOff + magenta + " FUNCTION: " + name_ + colorOff + blue + " FILE: " + filename_ + colorOff + red + " : " + message_ + colorOff;
};

#endif