#pragma once

#include <vector>
#include "Eigen-3.3/Eigen/Core"

class MPC {
 public:
  // Timestep length
  size_t N;
  // Timestemp duration
  double dt;
  
 public:
  MPC(size_t N, double dt);

  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuations.
  std::vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
};
