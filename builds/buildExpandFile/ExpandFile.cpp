#include "fundamental.hpp"

int main(int argc, char** argv){
  std::string filename_org = argv[1];  
  std::string filename = argv[2];
  std::cout << Blue << argv[1] << " -> " << argv[2] << reset << std::endl;
  std::ifstream in_strm;
  in_strm.open(filename_org,std::ios::in);
  //////////////
  if(!in_strm){
    std::cout << Magenta << "file :  " << filename_org  << reset << std::endl;    
    std::cout << Magenta << " can not be opened " << reset << std::endl;    

    if (in_strm.bad())
      std::cout << "I/O error while reading\n";
    else if (in_strm.eof())
      std::cout << "End of file reached successfully\n";
    else if (in_strm.fail())
      std::cout << "Non-integer data encountered\n";

    abort();
  }
  std::ofstream out_strm;
  out_strm.open(filename,std::ios::out);
  if(!out_strm){
    std::cout << Magenta << "file :  " << filename  << reset << std::endl;    
    std::cout << Magenta << " can not be opened " << reset << std::endl;


    if (out_strm.bad())
      std::cout << "I/O error while reading\n";
    else if (out_strm.eof())
      std::cout << "End of file reached successfully\n";
    else if (out_strm.fail())
      std::cout << "Non-integer data encountered\n";

    abort();
  }

  
  ///////////////
  int counter = 0;
  while(!in_strm.eof())// loop until the end of file 
    {
      std::string read;
      std::getline(in_strm,read);

      // if(read.find("Expand(")!=std::string::npos && read.find(")")!=std::string::npos){
      //   std::string file = StringSplit(read,{"Expand(\"","\")"})[0];

      if(!read.empty() && read.find("ExpandFile(\"")!=std::string::npos && read.find("\")")!=std::string::npos){
        std::string file = StringSplit(read,{"ExpandFile(\"","\")"})[0];
      
	std::ifstream ifs(file);
	if(!ifs)
	  {
	    Print(counter,Red);
	    std::cout << Green << "in filename :  " << filename  << reset << std::endl;    
	    std::cout << Green << "try to open file :  " << std::endl;
	    std::cout << Green << file  << reset << std::endl;


	    std::cout << green << "while reading :  " << std::endl;
	    std::cout << Green << read  << reset << std::endl;
	    std::cout << green << StringSplit(read,{"ExpandFile(\"","\")"})  << reset << std::endl;	    
	    
	    std::cout << Green << " can not be opened " << reset << std::endl;
	    abort();
	  }	
	std::string content( (std::istreambuf_iterator<char>(ifs) ),
			     (std::istreambuf_iterator<char>()    ) );
	
	out_strm << content << std::endl;
	out_strm << "filename: " << file << std::endl;
      }else{
	out_strm << read << std::endl;	  
      }
      counter++;
    }
  in_strm.close();
  out_strm.close();  
  
  return 0;
}
