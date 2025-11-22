#include <fstream>
#include <sstream>
#include <string>
#include "basic.hpp"

/*DOC_EXTRACT JSON

# JSONクラス

\insert{basic::JSON}

*/

int main() {
   // sample.jsonを文字として読み込み表示
   {
      std::string filename = "./sample.json";
      std::ifstream istrm(filename);

      std::string lines;
      if (!istrm.is_open())
         std::cout << "can not open " << filename << std::endl;
      else {
         std::cout << "----------------------------------------" << std::endl;
         std::cout << magenta << "sample.jsonを文字として読み込み表示" << colorOff << std::endl;
         std::string line;
         while (!istrm.eof()) {
            getline(istrm, line);
            lines += line;
            std::cout << line << std::endl;
         }
         std::cout << "----------------------------------------" << std::endl;
         auto json = parseJSON(lines);
         for (auto &[a, b] : json)
            std::cout << a << ", " << b << std::endl;
         std::cout << "----------------------------------------" << std::endl;
      }
   }

   {
      std::cout << red << "1. 文字列からJSONをコンストラクト" << colorOff << std::endl;
      JSON json("./sample.json");
      std::cout << json["translate"] << std::endl;
      std::cout << stod(json["price"]) << std::endl;
      std::cout << stod(json["translate"]) << std::endl;
      std::cout << stod(json["scale"]) << std::endl;
      std::cout << stod(json["rotation"]) << std::endl;
      std::cout << stob(json["extra_cheese"]) << std::endl;
      std::cout << stob(json["delivery"]) << std::endl;
      std::cout << stob(json["phone"]) << std::endl;
      std::cout << stob(json["mooring"]) << std::endl;
   }

   {
      std::cout << red << "2. ifstreamでJSONをコンストラクト" << colorOff << std::endl;
      JSON json(std::ifstream("./sample.json"));
      json["price"] = {"10."};
      std::cout << json["translate"] << std::endl;
      std::cout << stod(json["price"]) << std::endl;
      std::cout << stod(json["translate"]) << std::endl;
      std::cout << stod(json["scale"]) << std::endl;
      std::cout << stod(json["rotation"]) << std::endl;
      std::cout << stob(json["extra_cheese"]) << std::endl;
      std::cout << stob(json["delivery"]) << std::endl;
      std::cout << stob(json["phone"]) << std::endl;

      std::ofstream os("./output.json");
      os << json;
      os.close();
   }

   JSONoutput jsonout;

   jsonout.push("time", 1);
   jsonout.push("time", 1);
   jsonout.push("time", 1);
   jsonout.push("time", 1);
   jsonout.push("volume", 1);
   jsonout.push("volume", 1);
   jsonout.push("volume", 1);
   std::ofstream os("./output2.json");
   jsonout.output(os);
   os.close();

   /* ------------------------------------------------------ */

   // 内容表示方法
   std::cout << "内容表示方法" << std::endl;
   JSON json(std::ifstream("./output.json"));
   std::cout << ToString(json) << std::endl;
   // for (auto &[key, value] : json())
   //    std::cout << key << ", " << value << std::endl;
}