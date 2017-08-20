#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

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

    // Initialize the cost to zero
    fg[0] = 0;

    // Minimize the cte, velocity and epsilon errors
    for (int t = 0; t < N - 1; t++){
    	fg[0] += 2000 * CppAD::pow(vars[cte_start + t], 2);
    	fg[0] += 8000 * CppAD::pow(vars[epsi_start + t], 2);
    	fg[0] += CppAD::pow(vars[v_start + t] - ref_v, 2);
    }

    // Minimize the use of actuators.
    for (int t = 0; t < N - 1; t++) {
      fg[0] += CppAD::pow(vars[delta_start + t], 2);
      fg[0] += CppAD::pow(vars[a_start + t], 2);
    }

    // Minimize the value gap between sequential actuations.
    for (int t = 0; t < N - 2; t++) {
      fg[0] += 3200000 * CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
      fg[0] += CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
    }

    // Initialize model
    fg[1 + x_start] = vars[x_start];
		fg[1 + y_start] = vars[y_start];
		fg[1 + psi_start] = vars[psi_start];
		fg[1 + v_start] = vars[v_start];
		fg[1 + cte_start] = vars[cte_start];
		fg[1 + epsi_start] = vars[epsi_start];

		// Initialize the constraints
		for (int t = 1; t < N; t++) {
		  // The state at time t+1 .
		  AD<double> x1 = vars[x_start + t];
		  AD<double> y1 = vars[y_start + t];
		  AD<double> psi1 = vars[psi_start + t];
		  AD<double> v1 = vars[v_start + t];
		  AD<double> cte1 = vars[cte_start + t];
		  AD<double> epsi1 = vars[epsi_start + t];

		  // The state at time t.
		  AD<double> x0 = vars[x_start + t - 1];
		  AD<double> y0 = vars[y_start + t - 1];
		  AD<double> psi0 = vars[psi_start + t - 1];
		  AD<double> v0 = vars[v_start + t - 1];
		  AD<double> cte0 = vars[cte_start + t - 1];
		  AD<double> epsi0 = vars[epsi_start + t - 1];

		  // Only consider the actuation at time t.
		  AD<double> delta0 = vars[delta_start + t - 1];
		  AD<double> a0 = vars[a_start + t - 1];

		  // Polynomial calculations
		  AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * x0 * x0 + coeffs[3] * x0 * x0 * x0;
		  AD<double> psides0 = CppAD::atan(coeffs[1] + coeffs[2] * x0 + coeffs[3] * x0 * x0);

		  // Here's `x` to get you started.
		  // The idea here is to constraint this value to be 0.
		  //
		  // Recall the equations for the model:
		  // x_[t] = x[t-1] + v[t-1] * cos(psi[t-1]) * dt
		  // y_[t] = y[t-1] + v[t-1] * sin(psi[t-1]) * dt
		  // psi_[t] = psi[t-1] + v[t-1] / Lf * delta[t-1] * dt
		  // v_[t] = v[t-1] + a[t-1] * dt
		  // cte[t] = f(x[t-1]) - y[t-1] + v[t-1] * sin(epsi[t-1]) * dt
		  // epsi[t] = psi[t] - psides[t-1] + v[t-1] * delta[t-1] / Lf * dt
		  fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
		  fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
		  fg[1 + psi_start + t] = psi1 - (psi0 + v0 * delta0 / Lf * dt);
		  fg[1 + v_start + t] = v1 - (v0 + a0 * dt);
		  fg[1 + cte_start + t] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
		  fg[1 + epsi_start + t] = epsi1 - ((psi0 - psides0) + v0 * delta0 / Lf * dt);
		}

  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {
	// Vectors to show the green spline in the simulator
	mpc_x.resize(N-1);
	mpc_y.resize(N-1);
}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  size_t n_vars = 6*N + 2*(N-1);
  // TODO: Set the number of constraints
  size_t n_constraints = 6*N;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // TODO: Set lower and upper limits for variables.

  // Bounds for state variables
  for (int i = 0; i < delta_start; i++){
  	vars_lowerbound[i] = -1.0e-19;
  	vars_upperbound[i] = 1.0e19;
  }

  // Bounds for control signal variables
  for (int i = delta_start; i < a_start; i++){
  	vars_lowerbound[i] = -0.436332; // -25 degrees
  	vars_upperbound[i] = 0.436332;  //  25 degrees
  }
  for (int i = a_start; i < n_vars; i++){
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

  // object that computes objective and constraints
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
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.

  // This adds values to plot into the simulator
  for (int i = 1; i < N; i++){
  	mpc_x[i] = solution.x[i];

  }
  return {solution.x[delta_start], solution.y[a_start]};
}
