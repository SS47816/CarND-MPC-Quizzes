#include <vector>
#include "Eigen-3.3/Eigen/QR"
#include "helpers.h"
#include "matplotlibcpp.h"
#include "MPC.h"

namespace plt = matplotlibcpp;

using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

int main() {
  MPC mpc;
  int iters = 50;

  VectorXd ptsx(12);
  VectorXd ptsy(12);
  // ptsx << -100, 100;
  // ptsy << -1, -1;

  // ptsx << -2.02222, -1.72357, -1.41589, -1.09661, -0.763399, -0.414111, -0.0468359,
  // 		0.340131, 0.748285, 1.17891, 1.17891, 1.17891;

  // ptsx << -0.0435319, -0.125675, -0.208346, -0.291873, -0.376647, -0.463111, -0.551755,
		// -0.643102, -0.737702, -0.83612, -0.83612, -0.83612;

  ptsx << -1.89952, -1.59594, -1.28631, -0.96791, -0.638175, -0.294751,
			0.0645291, 0.44164, 0.838364, 1.25629, 1.69682, 1.69682;

  ptsy << 0.00573502, 0.0125499, 0.0213456, 0.0328391, 0.0475759, 0.0659421, 
  			0.0881765, 0.114384, 0.144545, 0.178533, 0.21612, 0.21612;
 

	auto coeffs = polyfit(ptsx, ptsy, 3);

  /**
   * TODO: fit a polynomial to the above x and y coordinates
   * The polynomial is fitted to a straight line 
   * so a ploynomial with order of 1 is sufficient
   */

  // Start state
  // NOTE: free feel to play around with these
  double x = 0;
  double y = 0;
  double psi = 0.1;
  double v = 0.2;
  
  /**
   * TODO: calculate the cross track error
   * The cross track error is calculated by evaluating 
   * at polynomial at (x, f(x)) and subtracting y.
   */
  double cte = -polyeval(coeffs, 0);
  
  /**
   * TODO: calculate the orientation error
   * Due to the sign starting at 0, the orientation error is -f'(x)
   * derivative of coeffs[0] + coeffs[1] * x -> coeffs[1]
   */
  double epsi = -atan(coeffs[1]);

  std::cout << "Coefficients: " << std::endl;
  for (int i = 0; i < coeffs.size(); i++)
  {
    std::cout << coeffs(i) << std::endl;
  }
  std::cout << "CTE: " << cte << " EPSI: " << epsi << std::endl;

  VectorXd state(6);
  state << x, y, psi, v, cte, epsi;

  vector<double> x_vals = {state[0]};
  vector<double> y_vals = {state[1]};
  vector<double> psi_vals = {state[2]};
  vector<double> v_vals = {state[3]};
  vector<double> cte_vals = {state[4]};
  vector<double> epsi_vals = {state[5]};
  vector<double> delta_vals = {};
  vector<double> a_vals = {};

  for (size_t i = 0; i < iters; ++i) {
    cout << "Iteration " << i << endl;

    auto vars = mpc.Solve(state, coeffs);

    x_vals.push_back(vars[0]);
    y_vals.push_back(vars[1]);
    psi_vals.push_back(vars[2]);
    v_vals.push_back(vars[3]);
    cte_vals.push_back(vars[4]);
    epsi_vals.push_back(vars[5]);

    delta_vals.push_back(vars[6]);
    a_vals.push_back(vars[7]);

    state << vars[0], vars[1], vars[2], vars[3], vars[4], vars[5];
    cout << "x = " << vars[0] << endl;
    cout << "y = " << vars[1] << endl;
    cout << "psi = " << vars[2] << endl;
    cout << "v = " << vars[3] << endl;
    cout << "cte = " << vars[4] << endl;
    cout << "epsi = " << vars[5] << endl;
    cout << "delta = " << vars[6] << endl;
    cout << "a = " << vars[7] << endl;
    cout << endl;
  }

  // Plot values
  // NOTE: feel free to play around with this.
  // It's useful for debugging!
  plt::subplot(6, 1, 1);
  plt::title("Track");
  plt::plot(x_vals, y_vals);
  
  plt::subplot(6, 1, 2);
  plt::title("CTE");
  plt::plot(cte_vals);

  plt::subplot(6, 1, 3);
  plt::title("EPSI");
  plt::plot(epsi_vals);

  plt::subplot(6, 1, 4);
  plt::title("Delta (Radians)");
  plt::plot(delta_vals);

  plt::subplot(6, 1, 5);
  plt::title("Velocity");
  plt::plot(v_vals);
  
  plt::subplot(6, 1, 6);
  plt::title("Acceleration (m/s^2)");
  plt::plot(a_vals);

  plt::show();
}