#include "EulerSystem.h"

int main()
{
  EulerSystem E;

  E.run();

  E.outputFile("RiemannExactResults/Test5/Test5_400cells.dat");
}
