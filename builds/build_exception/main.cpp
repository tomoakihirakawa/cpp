#include <iostream>
#include <string>
#include "fundamental_constants.hpp"

class error_message{
 public:
  std::string filename;
  std::string name;
  int line;
  std::string message;
  error_message(const std::string& filename_, const std::string& name_, const int line_, const std::string message_):
    filename(filename_),
    name(name_),
    line(line_),
    message(message_)
    {};

  std::string what() const {
    return
      reset + "-------------------------------------------------------\n" +
      reset + "    line: " + Green + std::to_string(this->line) + "\n"  +
      reset + "function: " + Red + name + "\n" +
      reset + "    file: " + red + filename + "\n" +
      reset + " message: " + Red + message + reset ;
  };  
};


void func(int i)
{
  if(i==0){

    throw(error_message(__FILE__,__PRETTY_FUNCTION__,__LINE__,""));
    
  }else if(i==2){

    throw(error_message(__FILE__,__PRETTY_FUNCTION__,__LINE__,""));
    
  }else{
    
    throw(10);
    
  }
  
}


int main(){

  try{
    
    try{
      func(0);    
    }

    catch (error_message e){
      std::cerr << e.what() << std::endl;
      throw(error_message(__FILE__,__PRETTY_FUNCTION__,__LINE__,""));    
    }

    catch (int e) {
      std::cerr << e << std::endl;
    }

  }
  catch (error_message e){
    std::cerr << e.what() << std::endl;
    throw(error_message(__FILE__,__PRETTY_FUNCTION__,__LINE__,""));    
  }

  
  
}
