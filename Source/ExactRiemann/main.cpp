#include "EulerSystem.h"

int main()
{
  EulerSystem E;

  E.run();

  E.outputFile("results/data.dat");
}
