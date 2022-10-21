#include <algorithm>
#include <iostream>
#include <vector>

int main()
{
	std::vector<int> vec = {3, 1, 4, 6, 5};

	// ソート：二分探索はソートされた配列に対して適用する
	sort(vec.begin(), vec.end()); // 1, 3, 4, 5, 6

	for (auto i = 0; i < 10; i++)
	{
		vec.insert(std::lower_bound(vec.cbegin(), vec.cend(), i), i);
		for (auto &v : vec)
			std::cout << v << " ,";
		std::cout << "\n";
	}

	// 0 - 10 の整数が vec に含まれているか検索
	for (int key = 0; key <= 10; ++key)
	{
		if (binary_search(vec.begin(), vec.end(), key))
		{
			std::cout << key << ": "
					  << "found" << std::endl;
		}
		else
		{
			std::cout << key << ": "
					  << "not found" << std::endl;
		}
	}

	//チェックその２

	// vec = {3, 1, 4, 6, 5, 0, -1};
	vec = {};
	// ソート：二分探索はソートされた配列に対して適用する
	sort(vec.begin(), vec.end()); // 1, 3, 4, 5, 6

	std::vector<int>::iterator it;

	for (auto i = 0; i < 10; i++)
	{
		if ((it = std::lower_bound(vec.begin(), vec.end(), i)) != vec.end() || vec.empty())
			vec.insert(it, i);
		for (auto &v : vec)
			std::cout << v << " ,";
		std::cout << "\n";
	}

	// 0 - 10 の整数が vec に含まれているか検索
	for (int key = 0; key <= 20; ++key)
	{
		it = std::lower_bound(vec.begin(), vec.end(), key);
		if (it != vec.end())
		{
			std::cout << key << ": "
					  << "found" << std::endl;
		}
		else
		{
			std::cout << key << ": "
					  << "not found" << std::endl;
		}
	}

	return 0;
}
