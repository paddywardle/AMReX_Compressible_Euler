#include "EulerSystem.h"

int main()
{
  EulerSystem E;

  E.run();

  E.outputFile("RiemannExactResults/Test3/Test3_100cells.dat");
}
