/*
  A 1D finite-difference heat/diffusion equation solver

  Computation is local (single compute node) and runs
  the kernel in parallel using a `forall` loop

  Values of the `config const` variables can be modified in
  the command line (e.g., `./heat_1D --nt=100`)
*/

// declare configurable constants with default values
config const nx = 4096,     // number of grid points in x
             nt = 50,       // number of time steps
             alpha = 0.25;  // diffusion constant

// define a domain and subdomain to describe the grid and its interior
const indices = {0..<nx},
      indicesInner = {1..<nx-1}; // equivalent to `= indices.expand(-1)`

// define an array over the above domain
var u : [indices] real;

// set up initial conditions
u = 1.0;
u[nx/4..nx/2] = 2.0;

// create a temporary copy of 'u' to store the previous time step
var un = u;

// iterate for 'nt' time steps
for 1..nt {
  // swap the arrays to prepare for the next time step
  u <=> un;

  // update the solution over the interior of the domain in parallel
  forall i in indicesInner do
    u[i] = un[i] + alpha * (un[i-1] - 2 * un[i] + un[i+1]);
}

// print final results
const mean = (+ reduce u) / u.size,
      stdDev = sqrt((+ reduce (u - mean)**2) / u.size);
writeln("mean: ", mean, " stdDev: ", stdDev);
