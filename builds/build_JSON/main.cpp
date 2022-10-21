#include <fstream>
#include <sstream>
#include <string>
#include "basic.hpp"

// 2021/03/31JSONとして導入した
//  class myJSON {
//  	using V_S = std::vector<std::string>;
//  	std::map<std::string, V_S> map_S_S;

//    public:
// 	myJSON(const std::string& str_IN) : map_S_S() {
// 		V_S SEP = {" ", "\t", "\n", "\r"};
// 		std::string trimed = StringTrim(str_IN, SEP);
// 		std::string L = "[", R = "]", Col = ":", Cam = ",", Lb = "{", Rb = "}";
// 		bool is_array = false;
// 		std::string array = "", value = "";
// 		V_S vstr;
// 		for (auto it = trimed.begin(); it != trimed.end(); it++) {
// 			//array catcher
// 			if (is_array || L.find(*it) != std::string::npos) {
// 				if ((R.find(*it) != std::string::npos)) {
// 					is_array = false;
// 					// std::cout << blue << array << std::endl;
// 					vstr.emplace_back(array);
// 					array.clear();
// 				} else {
// 					if (is_array)
// 						array.push_back(*it);
// 					is_array = true;
// 				}
// 			} else {
// 				if ((Lb.find(*it) != std::string::npos) ||
// 				    (Cam.find(*it) != std::string::npos) ||
// 				    (Col.find(*it) != std::string::npos) ||
// 				    (Rb.find(*it) != std::string::npos)) {
// 					if (!value.empty()) {
// 						// std::cout << red << value << std::endl;
// 						vstr.emplace_back(value);
// 						value.clear();
// 					}
// 				} else {
// 					value.push_back(*it);
// 				}
// 			}
// 		}
// 		for (auto i = 0; i < vstr.size() - 1; i += 2) {
// 			this->map_S_S[StringTrim(vstr[i], {"\""})] = StringSplit(StringTrim(vstr[i + 1], {"\""}), {","});
// 		}
// 	};

// 	std::map<std::string, V_S> operator()() const {
// 		return this->map_S_S;
// 	};
// };

// int main()
// {
// 	//sample.jsonを文字として読み込み，表見
// 	{
// 		std::string filename = "./sample.json";
// 		std::ifstream istrm(filename);

// 		std::string lines;
// 		if (!istrm.is_open())
// 			std::cout << "can not open " << filename << std::endl;
// 		else
// 		{
// 			std::cout << "----------------------------------------" << std::endl;
// 			std::string line;
// 			while (!istrm.eof())
// 			{
// 				getline(istrm, line);
// 				lines += line;
// 				std::cout << line << std::endl;
// 			}
// 			std::cout << "----------------------------------------" << std::endl;
// 			JSON json(lines);
// 			for (auto &[a, b] : json())
// 				std::cout << a << ", " << b << std::endl;
// 			std::cout << "----------------------------------------" << std::endl;
// 		}
// 	}

// 	{
// 		std::cout << magenta << "stringでJSONをコンストラクト" << reset << std::endl;
// 		std::ifstream istrm("./sample.json");
// 		if (!istrm.is_open())
// 		{
// 			std::cout << "can not open" << std::endl;
// 		}
// 		else
// 		{
// 			std::stringstream ss;
// 			ss << istrm.rdbuf();
// 			JSON json(ss.str());
// 			auto parsed = json();
// 			std::cout << parsed["translate"] << std::endl;
// 			std::cout << stod(parsed["translate"]) << std::endl;
// 			std::cout << stod(parsed["scale"]) << std::endl;
// 			std::cout << stod(parsed["rotation"]) << std::endl;
// 		}
// 	}

// 	{
// 		std::cout << red << "ifstreamでJSONをコンストラクト" << reset << std::endl;
// 		std::ifstream istrm("./sample.json");
// 		JSON json(istrm);
// 		auto parsed = json();
// 		std::cout << parsed["translate"] << std::endl;
// 		std::cout << stod(parsed["translate"]) << std::endl;
// 		std::cout << stod(parsed["scale"]) << std::endl;
// 		std::cout << stod(parsed["rotation"]) << std::endl;
// 	}
// }

// 2021/09/16
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
      std::cout << magenta << "1. filenameでJSONをコンストラクト" << colorOff << std::endl;
      JSON json("./sample.json");
      std::cout << json["translate"] << std::endl;
      std::cout << stod(json["price"]) << std::endl;
      std::cout << stod(json["translate"]) << std::endl;
      std::cout << stod(json["scale"]) << std::endl;
      std::cout << stod(json["rotation"]) << std::endl;
      std::cout << stob(json["extra_cheese"]) << std::endl;
      std::cout << stob(json["delivery"]) << std::endl;
      std::cout << stob(json["phone"]) << std::endl;
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
   JSON json(std::ifstream("./sample.json"));
   for (auto &[key, value] : json())
      std::cout << key << ", " << value << std::endl;
}