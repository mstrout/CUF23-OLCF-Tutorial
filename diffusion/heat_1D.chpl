/*
  A 1D finite-difference heat/diffusion equation solver

  Computation is local (single compute node) and runs
  the kernel in parallel using a `forall` loop

  Values of the `config const` variables can be modified in
  the command line (e.g., `./heat_1D --nt=100`)
*/

// declare configurable constants with default values
config const xLen = 2.0,    // length of the grid in x
             nx = 31,       // number of grid points in x
             nt = 50,       // number of time steps
             sigma = 0.25,  // stability parameter
             nu = 0.05;     // viscosity

// define non-configurable constants
const dx : real = xLen / (nx - 1),       // grid spacing in x
      dt : real = sigma * dx / nu;       // time step size

// define a domain and subdomain to describe the grid and its interior
const indices = {0..<nx},
      indicesInner = {1..<nx-1}; // equivalent to `= indices.expand(-1)`

// define an array over the above domain
var u : [indices] real;

// set up initial conditions
u = 1.0;
u[(0.5 / dx):int..<(1.0 / dx + 1):int] = 2;

// create a temporary copy of 'u' to store the previous time step
var un = u;

// iterate for 'nt' time steps
for 1..nt {
  // swap the arrays to prepare for the next time step
  u <=> un;

  // update the solution over the interior of the domain in parallel
  forall i in indicesInner do
    u[i] = un[i] + nu * dt / dx**2 *
            (un[i-1] - 2 * un[i] + un[i+1]);
}

// print final results
const mean = (+ reduce u) / u.size,
      stdDev = sqrt((+ reduce (u - mean)**2) / u.size);
writeln("mean: ", mean, " stdDev: ", stdDev);