#ifndef INCL_fundamental_exception
#define INCL_fundamental_exception

#include <iostream>
#include <string>
#include <stdexcept>
#include "fundamental_constants.hpp"

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
    return reset + "----------------------------------------------------\n" +
           reset + "    LINE: " + Green + std::to_string(this->line) + "\n" +
           reset + "FUNCTION: " + Magenta + this->name + "\n" +
           reset + "    FILE: " + red + this->filename + "\n" +
           reset + " message: " + Red + this->message + reset + "\n" +
           reset + "----------------------------------------------------\n";
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
  return reset + green + " LINE: " + std::to_string(line_) + reset + magenta + " FUNCTION: " + name_ + reset + blue + " FILE: " + filename_ + reset + red + " : " + message_ + reset;
};

#endif