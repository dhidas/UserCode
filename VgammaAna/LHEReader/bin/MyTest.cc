#include <iostream>


int MyTest ()
{
  std::cout << "Hello world!" << std::endl;
  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 1) {
    std::cerr << "Usage: " << argv[0] << " " << std::endl;
    return 1;
  }

  MyTest();

  return 0;
}
