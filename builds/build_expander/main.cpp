#include "fundamental.hpp"

using V_s = std::vector<std::string>;
using VV_s = std::vector<std::vector<std::string>>;

std::string StringConcat(const V_s &vs)
{
	std::string ret = "";
	for (const auto &s : vs)
		ret += s;
	return ret;
};

bool StringContainsQ(const std::string read, const V_s &vs)
{
	if (read.empty())
		return false;
	for (const auto &s : vs)
		if (read.find(s) == std::string::npos)
			return false;
	return true;
};

int main(int argc, char **argv)
{
	std::string filename_org = argv[1];
	std::string filename = argv[2];

	if (filename_org.compare(filename) == 0)
	{
		std::cout << _Red << "input and output are both " << filename << std::endl;
		std::cout << _Red << "Please change name" << filename << std::endl;		
		return 0;
	}

	std::cout << Blue << argv[1] << " -> " << argv[2] << reset << std::endl;
	std::ifstream in_strm;
	in_strm.open(filename_org, std::ios::in);
	//////////////
	if (!in_strm)
	{
		std::cout << Magenta << "file :  " << filename_org << reset << std::endl;
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
	out_strm.open(filename, std::ios::out);
	if (!out_strm)
	{
		std::cout << Magenta << "file :  " << filename << reset << std::endl;
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
	while (!in_strm.eof()) // loop until the end of file
	{
		std::string read;
		std::getline(in_strm, read);
		if (StringContainsQ(read, {"ExpandPart(", ",", ")"}))
		{
			read = StringTrim(read, {" "});
			V_s file_part = StringSplit(read, {"ExpandPart(", ",", ")"});
			std::string file = file_part[0];
			std::string part = file_part[1];

			Print("file:" + file, red);
			Print("part:" + part, red);

			std::ifstream ifs(file); //open
			if (!ifs)
			{
				Print(counter, Red);
				Print(std::string("in filename :  ") + filename, Green);
				Print(std::string("try to open :  ") + file, Green);
				Print(std::string("while reading :  ") + read, green);
				Print(std::string("while reading :  ") + StringConcat(StringSplit(read, {"ExpandPart(", ",", ")"})), green);
				Print(std::string("can not be opened "), Green);
				abort();
			}

			std::string content((std::istreambuf_iterator<char>(ifs)),
								(std::istreambuf_iterator<char>()));

			V_s read_part = StringSplit(content, {part});

			if (read_part.size() < 2)
				break;

			// int i=0;
			// for(const auto& r:read_part){
			//   Print(i++,Red);
			//   std::cout << r << std::endl;
			// }

			out_strm << read_part[1] << std::endl;
			out_strm << "//filename: " << file << std::endl;
			out_strm << "//part: " << part << std::endl;
		}
		else if (StringContainsQ(read, {"Expand(", ",", ")"}))
		{
			read = StringTrim(read, {" "});
			V_s file_part = StringSplit(read, {"Expand(", ",", ")"});
			std::string file = file_part[0];
			std::string part = file_part[1];

			Print("file:" + file, red);
			Print("part:" + part, red);

			std::ifstream ifs(file); //open
			if (!ifs)
			{
				Print(counter, Red);
				Print(std::string("in filename :  ") + filename, Green);
				Print(std::string("try to open :  ") + file, Green);
				Print(std::string("while reading :  ") + read, green);
				Print(std::string("while reading :  ") + StringConcat(StringSplit(read, {"Expand(", ",", ")"})), green);
				Print(std::string("can not be opened "), Green);
				abort();
			}

			std::string content((std::istreambuf_iterator<char>(ifs)),
								(std::istreambuf_iterator<char>()));

			V_s read_part = StringSplit(content, {part});

			if (read_part.size() < 2)
				break;

			// int i=0;
			// for(const auto& r:read_part){
			//   Print(i++,Red);
			//   std::cout << r << std::endl;
			// }

			out_strm << read_part[1] << std::endl;
			// out_strm << "//filename: " << file << std::endl;
			// out_strm << "//part: " << part << std::endl;
		}
		else if (StringContainsQ(read, {"Expand(", ")"}))
		{
			read = StringTrim(read, {" "});
			V_s file_part = StringSplit(read, {"Expand(", ")"});
			std::string file = file_part[0];

			Print("file:" + file, red);

			std::ifstream ifs(file);
			if (!ifs)
			{
				Print(counter, Red);
				Print(std::string("in filename :  ") + filename, Green);
				Print(std::string("try to open :  ") + file, Green);
				Print(std::string("while reading :  ") + read, green);
				Print(std::string("while reading :  ") + StringConcat(StringSplit(read, {"Expand(", ")"})), green);
				Print(std::string("can not be opened "), Green);
				abort();
			}
			std::string content((std::istreambuf_iterator<char>(ifs)),
								(std::istreambuf_iterator<char>()));

			out_strm << content << std::endl;
		}
		else if (StringContainsQ(read, {"ExpandFile(", ")"}))
		{
			read = StringTrim(read, {" "});
			auto file_part = StringSplit(read, {"ExpandFile(", ")"});
			auto file = file_part[0];
			auto part = file_part[1];

			Print("file:" + file, red);

			std::ifstream ifs(file);
			if (!ifs)
			{
				Print(counter, Red);
				Print(std::string("in filename :  ") + filename, Green);
				Print(std::string("try to open :  ") + file, Green);
				Print(std::string("while reading :  ") + read, green);
				Print(std::string("while reading :  ") + StringConcat(StringSplit(read, {"ExpandFile(", ")"})), green);
				Print(std::string("can not be opened "), Green);
				abort();
			}
			std::string content((std::istreambuf_iterator<char>(ifs)),
								(std::istreambuf_iterator<char>()));

			out_strm << content << std::endl;
			out_strm << "//filename: " << file << std::endl;
		}
		else
		{
			out_strm << read << std::endl;
		}
		// Print(counter);
		if (counter > 100000)
			break;
		counter++;
	}

	in_strm.close();
	out_strm.close();

	return 0;
}
