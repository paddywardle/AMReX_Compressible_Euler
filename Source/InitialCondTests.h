#include <array>

typedef std::array<double,4> var_array;

class InitialCondTests
{
 public:

  InitialCondTests(){};
  
  var_array Test1D_1(double);

  var_array Test1D_2(double);

  var_array Test1D_3(double);

  var_array Test1D_4(double);

  var_array Test1D_5(double);

  var_array Test1D_1Y(double);

  var_array Test1D_2Y(double);

  var_array Test1D_3Y(double);

  var_array Test1D_4Y(double);

  var_array Test1D_5Y(double);

  var_array Test1_diag(double, double);

  var_array Test2_diag(double, double);

  var_array Test3_diag(double, double);

  var_array Test4_diag(double, double);

  var_array Test5_diag(double, double);
  
  /*
  var_array Test2D_1(double);

  var_array Test2D_2(double);

  var_array Test2D_3(double);

  var_array Test2D_4(double);

  var_array Test2D_5(double);
  */

};
