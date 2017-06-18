#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"
#include "utils.h"

using CppAD::AD;
using Eigen::VectorXd;
namespace ad = CppAD;

// Set the timestep length and duration
size_t N    =   10;
double dt   = 0.08;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
constexpr double Lf   = 2.67;

// Target values
double ref_cte    =     0;
double ref_epsi   =     0;
double ref_v      =   120;

// Loss weights
double Wcte       =  800.0;
double Wepsi      =  300.0;
double Wdelta     = 4000.0;
double Wddelta    =    1e7;
double Wv         =    1.0;
double Wa         =    1.0;
double Wda        =    1.0;

// For tuning the parameters from the commandline
// For example `env ref_v=200 ./mpc` will set the
// target speed to 200... and most probably will
// send the car to the river!
static void init_globals_from_environment() {
  char *value;

  value = getenv("N");
  if (value) N = atoi(value);
  value = getenv("dt");
  if (value) dt = atof(value);
  value = getenv("ref_v");
  if (value) ref_v = atof(value);

  value = getenv("Wcte");
  if (value) Wcte = atof(value);
  value = getenv("Wepsi");
  if (value) Wepsi = atof(value);
  value = getenv("Wv");
  if (value) Wv = atof(value);
  value = getenv("Wdelta");
  if (value) Wdelta = atof(value);
  value = getenv("Wa");
  if (value) Wa = atof(value);
  value = getenv("Wddelta");
  if (value) Wddelta = atof(value);
  value = getenv("Wda");
  if (value) Wda = atof(value);
}

// Variable offsets in the variables vector.
struct variable_offsets {
  size_t begin, end;

  variable_offsets(size_t begin, size_t count)
  : begin(begin), end(begin + count)
  {}

  // Syntactic sugar
  inline operator size_t() const {
    return begin;
  }

  // Ditto
  inline size_t operator()(int t) const {
    assert(t >= 0);
    assert(begin + t < end);
    return begin + t;
  }
};

// The solver takes all the state variables and actuator
// variables in a singular vector. This struct allows
// to handle that.
struct variables {
  // state
  const variable_offsets x, y, psi, v, cte, epsi;
  // actuators (steering and acceleration)
  const variable_offsets delta, a;
  // ctor
  explicit variables(size_t N)
    : x(0, N), y(N, N), psi(2*N, N), v(3*N, N), cte(4*N, N), epsi(5*N, N)
    , delta(6*N, N-1), a(7*N-1, N-1)
  {
    assert(count() == 8*N - 2);
  }

  inline size_t count() const { return a.end; }
  inline size_t state_dim() const { return 6; }
  inline size_t actuators_dim() const { return 2; }
};


// Cost and Constraints
class FG_eval {
 public:
  // variables offsets
  const variables var;
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs)
    : var(N), coeffs(coeffs) {}

  using ADdouble = AD<double>;
  typedef CPPAD_TESTVECTOR(ADdouble) ADvector;

  // `fg` a vector of the cost constraints,
  // `vars` is a vector of variable values (state & actuators)
  void operator()(ADvector& fg, const ADvector& vars) {

    // Cost
    // ======================================================

    // The cost is stored is the first element of `fg`.
    // Any additions to the cost should be added to `fg[0]`.
    fg[0] = 0;

    // Minimize the state cost
    for (int t = 0; t < N; t++) {
      fg[0] +=  Wcte  * ad::pow(vars[var.cte(t) ] - ref_cte,  2);
      fg[0] +=  Wepsi * ad::pow(vars[var.epsi(t)] - ref_epsi, 2);
      fg[0] +=  Wv    * ad::pow(vars[var.v(t)   ] - ref_v,    2);
    }

    // Minimize the use of actuators
    for (int t = 0; t < N - 1; t++) {
      fg[0] +=  Wdelta * ad::pow(vars[var.delta(t)], 2);
      fg[0] +=  Wa     * ad::pow(vars[var.a(t)    ], 2);
    }

    // Minimize the value gap between sequential actuations.
    for (int t = 0; t < N - 2; t++) {
      fg[0] +=  Wddelta * ad::pow(vars[var.delta(t+1)] - vars[var.delta(t)], 2);
      fg[0] +=  Wda     * ad::pow(vars[var.a(t+1)    ] - vars[var.a(t)    ], 2);
    }

    // Constraints
    // ======================================================

    // We add 1 to each of the starting indices due to cost being located at
    // index 0 of `fg`. This bumps up the position of all the other values.
    fg[1 + var.x   ] = vars[var.x   ];
    fg[1 + var.y   ] = vars[var.y   ];
    fg[1 + var.psi ] = vars[var.psi ];
    fg[1 + var.v   ] = vars[var.v   ];
    fg[1 + var.cte ] = vars[var.cte ];
    fg[1 + var.epsi] = vars[var.epsi];

    for (int t = 0; t < N-1; t++) {
      // The state at time t.
      ADdouble x0      = vars[var.x(t)     ];
      ADdouble y0      = vars[var.y(t)     ];
      ADdouble psi0    = vars[var.psi(t)   ];
      ADdouble v0      = vars[var.v(t)     ];
      ADdouble cte0    = vars[var.cte(t)   ];
      ADdouble epsi0   = vars[var.epsi(t)  ];

      // The state at time t+1.
      ADdouble x1      = vars[var.x(t+1)   ];
      ADdouble y1      = vars[var.y(t+1)   ];
      ADdouble psi1    = vars[var.psi(t+1) ];
      ADdouble v1      = vars[var.v(t+1)   ];
      ADdouble cte1    = vars[var.cte(t+1) ];
      ADdouble epsi1   = vars[var.epsi(t+1)];

      // Only consider the actuation at time t.
      ADdouble delta0  = vars[var.delta(t) ];
      ADdouble a0      = vars[var.a(t)     ];

      // f(x[t])
      ADdouble f0      = polyeval(coeffs, x0);
      // desired side
      ADdouble psides0 = ad::atan(coeffs[1] +
                                  (2 * coeffs[2]
                                   + 3 * coeffs[3] * x0) * x0);

      // Constraints for the kinematic model:
      fg[2 + var.x(t)   ] =    x1 - (x0 + v0 * ad::cos(psi0) * dt);
      fg[2 + var.y(t)   ] =    y1 - (y0 + v0 * ad::sin(psi0) * dt);
      fg[2 + var.v(t)   ] =    v1 - (v0 + a0 * dt);
      // Note that we will have to invert the sign of delta0 in the actuator
      fg[2 + var.psi(t) ] =  psi1 - (psi0 + v0 / Lf * delta0 * dt);
      fg[2 + var.epsi(t)] = epsi1 - ((psi0 - psides0) + v0 / Lf * delta0 * dt);
      fg[2 + var.cte(t) ] =  cte1 - ((f0 - y0) + v0 * ad::sin(epsi0) * dt);
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {
  init_globals_from_environment();
  // Log the set of parameters for reference.
  printf("N=%d  dt=%.2f  ref_v=%.0f  Wcte=%.1f  Wepsi=%.1f  Wv=%.1f  Wdelta=%.1f  Wa=%.1f  Wddelta=%.1f  Wda=%.1f\n",
         (int)N, dt, ref_v, Wcte, Wepsi, Wv, Wdelta, Wa, Wddelta, Wda);
}

MPC::~MPC() {}

bool MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs, Result &result) {
  typedef CPPAD_TESTVECTOR(double) Dvector;

  // Variables offsets
  const variables var(N);

  const size_t n_vars = var.count();
  const size_t n_constraints = N * var.state_dim();

  // Explicit names for clarity
  const double
    x    = state[0],
    y    = state[1],
    psi  = state[2],
    v    = state[3],
    cte  = state[4],
    epsi = state[5];

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) vars[i] = 0;

  // Set the initial variable values
  vars[var.x   ] = x;
  vars[var.y   ] = y;
  vars[var.psi ] = psi;
  vars[var.v   ] = v;
  vars[var.cte ] = cte;
  vars[var.epsi] = epsi;

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);

  // Set lower and upper limits for variables.

  // Set all non-actuators upper and lowerlimits
  // to the max negative and positive values.
  for (int i = 0; i < var.delta.begin; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = +1.0e19;
  }

  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  for (int i = var.delta.begin; i < var.delta.end; i++) {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] =  0.436332;
  }

  // Acceleration/decceleration upper and lower limits.
  for (int i = var.a.begin; i < var.a.end; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] =  1.0;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  // Initial state is constant
  constraints_lowerbound[var.x   ] = x;
  constraints_lowerbound[var.y   ] = y;
  constraints_lowerbound[var.psi ] = psi;
  constraints_lowerbound[var.v   ] = v;
  constraints_lowerbound[var.cte ] = cte;
  constraints_lowerbound[var.epsi] = epsi;

  constraints_upperbound[var.x   ] = x;
  constraints_upperbound[var.y   ] = y;
  constraints_upperbound[var.psi ] = psi;
  constraints_upperbound[var.v   ] = v;
  constraints_upperbound[var.cte ] = cte;
  constraints_upperbound[var.epsi] = epsi;

  // object that computes the objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  ad::ipopt::solve_result<Dvector> solution;

  // solve the problem
  ad::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  bool ok = solution.status == ad::ipopt::solve_result<Dvector>::success;

  // Return the result with the cost, the actuator values and the predicted x, y
  result.cost = solution.obj_value;
  result.delta = solution.x[var.delta(0)];
  result.a = solution.x[var.a(0)];
  // next x, y coordinates
  result.x.clear();
  result.y.clear();
  for (int t = 1; t < N; t++) {
    result.x.push_back(solution.x[var.x(t)]);
    result.y.push_back(solution.x[var.y(t)]);
  }

  return ok;
}
