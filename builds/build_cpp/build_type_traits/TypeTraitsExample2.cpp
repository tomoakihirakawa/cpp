// TypeTraitsExample2.cpp
#include <iostream>
#include <type_traits>

template <typename T>
void printType(T&& input) {
   using ValueType = typename std::remove_reference<T>::type;

   std::cout << "Input: " << input << std::endl;
   std::cout << "Type of input: " << typeid(input).name() << std::endl;
   std::cout << "Type of ValueType: " << typeid(ValueType).name() << std::endl;
}

int main() {
   int x = 42;
   int& xRef = x;
   int&& xRvalueRef = 43;

   std::cout << "Passing an lvalue reference: " << std::endl;
   printType(xRef);

   std::cout << "\nPassing an rvalue reference: " << std::endl;
   printType(std::move(xRvalueRef));

   return 0;
}
