#include "EulerSystem.h"

int main()
{
  EulerSystem E;

  E.run();

  E.outputFile("RiemannExactResults/Test1/Test1_100cells.dat");
}
