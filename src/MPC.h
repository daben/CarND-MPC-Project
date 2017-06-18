#pragma once

#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

class MPC {
 public:
  MPC();

  virtual ~MPC();

  struct Result {
    double cost;
    double delta;
    double a;
    vector<double> x;
    vector<double> y;
  };
  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuations.
  bool Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs, Result &result);
};
