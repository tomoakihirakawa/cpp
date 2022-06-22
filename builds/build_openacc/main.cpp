#include <iostream>
#include <vector>
#include <tuple>
#include <string>

// Hello_World.c
void Print_Hello_World()
{
#pragma acc kernels
  for (int i = 0; i < 1000; i++)
  {
    printf("Hello World!\n");
  }
};

int main()
{

  Print_Hello_World();
};