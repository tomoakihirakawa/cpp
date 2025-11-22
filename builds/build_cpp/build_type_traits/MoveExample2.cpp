#include <iostream>
#include <string>

class Resource {
  public:
   Resource(const std::string& data) : data_(data) {
      std::cout << "Resource created: " << data_ << std::endl;
   }

   Resource(const Resource& other) : data_(other.data_) {
      std::cout << "Resource copied: " << data_ << std::endl;
   }

   Resource(Resource&& other) : data_(std::move(other.data_)) {
      std::cout << "Resource moved: " << data_ << std::endl;
   }

   ~Resource() {
      std::cout << "Resource destroyed: " << data_ << std::endl;
   }

  private:
   std::string data_;
};

void use_resource(Resource res) {
   std::cout << "Using resource: " << std::endl;
}

int main() {
   Resource res1("Resource1");
   Resource res2(res1);             // Copy constructor is called
   Resource res3(std::move(res1));  // Move constructor is called

   use_resource(res2);             // Copy constructor is called
   use_resource(std::move(res3));  // Move constructor is called

   return 0;
}
