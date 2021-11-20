#include "fundamental.hpp"

int main(int argc, char **argv)
{
	std::string filename_org = argv[1];
	std::string filename = argv[2];
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
		read = StringTrim(read, {" "});//ここを加えた
		if (!read.empty() && read.find("ExpandPart(") != std::string::npos && read.find(",") != std::string::npos && read.find(")") != std::string::npos)
		{
			std::vector<std::string> file_part = StringSplit(read, {"ExpandPart(", ",", ")"});
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
				Print(std::string("while reading :  ") + StringSplit(read, {"ExpandPart(", ",", ")"}), green);
				Print(std::string("can not be opened "), Green);
				abort();
			}

			std::string content((std::istreambuf_iterator<char>(ifs)),
								(std::istreambuf_iterator<char>()));

			std::vector<std::string> read_part = StringSplit(content, {part});

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
		else if (!read.empty() && read.find("Expand(") != std::string::npos && read.find(",") != std::string::npos && read.find(")") != std::string::npos)
		{
			std::vector<std::string> file_part = StringSplit(read), {"Expand(", ",", ")"});
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
				Print(std::string("while reading :  ") + StringSplit(read, {"Expand(", ",", ")"}), green);
				Print(std::string("can not be opened "), Green);
				abort();
			}

			std::string content((std::istreambuf_iterator<char>(ifs)),
								(std::istreambuf_iterator<char>()));

			std::vector<std::string> read_part = StringSplit(content, {part});

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
		else if (!read.empty() && read.find("Expand(") != std::string::npos && read.find(")") != std::string::npos)
		{
			std::vector<std::string> file_part = StringSplit(read, {"Expand(", ")"});
			std::string file = file_part[0];

			Print("file:" + file, red);

			std::ifstream ifs(file);
			if (!ifs)
			{
				Print(counter, Red);
				Print(std::string("in filename :  ") + filename, Green);
				Print(std::string("try to open :  ") + file, Green);
				Print(std::string("while reading :  ") + read, green);
				Print(std::string("while reading :  ") + StringSplit(read, {"Expand(", ")"}), green);
				Print(std::string("can not be opened "), Green);
				abort();
			}
			std::string content((std::istreambuf_iterator<char>(ifs)),
								(std::istreambuf_iterator<char>()));

			out_strm << content << std::endl;
		}
		else if (!read.empty() && read.find("ExpandFile(") != std::string::npos && read.find(")") != std::string::npos)
		{
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
				Print(std::string("while reading :  ") + StringSplit(read, {"ExpandFile(", ")"}), green);
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
