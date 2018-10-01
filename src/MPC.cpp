#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// Set the timestep length and duration
size_t N = 10;
double dt = 0.1;

// Flags for debugging
bool debug_mode = false;

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
const double Lf = 2.67;

// Reference Values
double ref_cte = 0;
double ref_epsi = 0;
double ref_v = 100;

// Constraint indicies
size_t x_start = 0;
size_t y_start = x_start + N;
size_t psi_start = y_start + N;
size_t v_start = psi_start + N;
size_t cte_start = v_start + N;
size_t epsi_start = cte_start + N;
size_t delta_start = epsi_start + N;
size_t a_start = delta_start + N - 1;


class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.

    if (debug_mode) { cout << "FG COST INIT START " << endl; };

    // Initialise Cost Function
    fg[0] = 0;

    // Reference State Costs
    // Adjust scalar values to penalise cost (high value > higher importance constraint)
    for (int i = 0; i < N; i++) {

      fg[0] += 2000 * CppAD::pow(vars[cte_start + i] - ref_cte, 2);
      fg[0] += 2000 * CppAD::pow(vars[epsi_start + i] - ref_epsi, 2);
      fg[0] += CppAD::pow(vars[v_start + i] - ref_v, 2);

    }

    // Steering and Acceleration Costs
    for (int i = 0; i < N - 1; i++) {

      fg[0] += 5 * CppAD::pow(vars[delta_start + i], 2);
      fg[0] += 5 * CppAD::pow(vars[a_start + i], 2);

    }

    // Change in Steering and Acceleration Costs
    for (int i = 0; i < N - 2; i++) {

      fg[0] += 200 * CppAD::pow(vars[delta_start + i + 1] - vars[delta_start + i], 2);
      fg[0] += 10 * CppAD::pow(vars[a_start + i + 1] - vars[a_start + i], 2);


    }

    if (debug_mode) { cout << "FG CONSTRAINT START" << endl; };

    // Setup Initial Constraints using initial conditions
    fg[x_start + 1] = vars[x_start];
    fg[y_start + 1] = vars[y_start];
    fg[psi_start + 1] = vars[psi_start];
    fg[v_start + 1] = vars[v_start];
    fg[cte_start + 1] = vars[cte_start];
    fg[epsi_start + 1] = vars[epsi_start];

    // Setup the rest of the constraints
    for (int i = 0; i < N - 1; i++) {

      // State at time t+1
      AD<double> x1 = vars[x_start + i + 1];
      AD<double> y1 = vars[y_start + i + 1];
      AD<double> psi1 = vars[psi_start + i + 1];
      AD<double> v1 = vars[v_start + i + 1];
      AD<double> cte1 = vars[cte_start + i + 1];
      AD<double> epsi1 = vars[epsi_start + i + 1];

      // State at time t
      AD<double> x0 = vars[x_start + i];
      AD<double> y0 = vars[y_start + i];
      AD<double> psi0 = vars[psi_start + i];
      AD<double> v0 = vars[v_start + i];
      AD<double> cte0 = vars[cte_start + i];
      AD<double> epsi0 = vars[epsi_start + i];

      // Actuation at time t.
      AD<double> delta0 = vars[delta_start + i];
      AD<double> a0 = vars[a_start + i];

      // Polynomial and Psi values (Desired state)
      AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * x0 * x0 + coeffs[3] * x0 * x0 * x0;
      AD<double> psides0 = CppAD::atan(3 * coeffs[3] * x0 * x0 + 2 * coeffs[2] * x0 + coeffs[1]);

      // Set Model Constraints to zero based on equations:

      // Recall the equations for the model:
      // x_[t] = x[t-1] + v[t-1] * cos(psi[t-1]) * dt
      // y_[t] = y[t-1] + v[t-1] * sin(psi[t-1]) * dt
      // psi_[t] = psi[t-1] + v[t-1] / Lf * delta[t-1] * dt
      // v_[t] = v[t-1] + a[t-1] * dt
      // cte[t] = f(x[t-1]) - y[t-1] + v[t-1] * sin(epsi[t-1]) * dt
      // epsi[t] = psi[t] - psides[t-1] + v[t-1] * delta[t-1] / Lf * dt
      fg[x_start + i + 2]    = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[y_start + i + 2]    = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[psi_start + i + 2]  = psi1 - (psi0 - v0 * delta0 / Lf * dt);
      fg[v_start + i + 2]    = v1 - (v0 + a0 * dt);
      fg[cte_start + i + 2]  = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
      fg[epsi_start + i + 2] = epsi1 - ((psi0 - psides0) - v0 * delta0 / Lf * dt);

    }



  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  //size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  // Initialise state
  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double v = state[3];
  double cte = state[4];
  double epsi = state[5];

  // Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  // 4 * 10 + 2 * 9
  size_t n_vars =  N * 6 + 2 * (N - 1);
  
  // Set the number of constraints
  size_t n_constraints = N * 6;


  if (debug_mode) { cout << "VAR INIT START" << endl; };


  // Initial value of the independent variables.
  // These should be 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  
  // Set lower and upper limits for variables.

  // Set all non-actuator upper and lowerlimits 
  // to the max negative and positive values
  for (int i = 0; i < delta_start; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  // Set upper and lower limits of delta (steering angle) to -25 to 25
  // degrees (in radians)
  for (int i = delta_start; i < a_start; i++) {
    vars_lowerbound[i] =  -0.436332 * Lf; // -1.0
    vars_upperbound[i] =  0.436332 * Lf;   // 1.0
  }

  // Set upper and lower limits to accelerations to -1 to 1.
  for (int i = a_start; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }

  if (debug_mode) { cout << "VAR INIT END" << endl; };

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  // Set constraints upper and lower bound to initial state
  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[psi_start] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[epsi_start] = epsi;
  
  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[psi_start] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = epsi;
  
  if (debug_mode) { cout << "CONSTRAINT INIT END" << endl; };

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  if (debug_mode) { cout << "FG EVAL END" << endl; }

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


  if (debug_mode) { cout << "OPTIONS END" << endl; };


  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  if (debug_mode) { cout << "START SOLVER" << endl; };

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);


  if (debug_mode) { cout << "END SOLVER" << endl; };

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.


  if (debug_mode) { cout << "Result Start" << endl; };

  vector<double> result;

  result.push_back(solution.x[delta_start]);
  result.push_back(solution.x[a_start]);

  for (int i = 0; i < N - 1; i++) {

    result.push_back(solution.x[x_start + i + 1]);
    result.push_back(solution.x[y_start + i + 1]);

  }

  if (debug_mode) { cout << "Result End" << endl; };

  return result;
}
