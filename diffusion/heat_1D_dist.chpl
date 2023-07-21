/*
  A distributed 1D finite-difference heat/diffusion equation solver

  Computation is executed over a 1D distributed array.
  The array distribution is managed by the `Block` distribution.
  Values are shared between neighboring locales using implicit
  communication. The `forall` loop manages task creation and
  synchronization across and within locales.

  Values of the `config const` variables can be modified in
  the command line (e.g., `./heat_1D --nt=100`)
*/

use BlockDist;

// declare configurable constants with default values
config const nx = 4096,     // number of grid points in x
             nt = 50,       // number of time steps
             alpha = 0.25;  // diffusion constant

// define distributed domains to describe the grid and its interior
const indices = Block.createDomain({0..<nx}),
      indicesInner = indices[{1..<nx-1}];

// define an array over the above domain
var u : [indices] real;

// set up initial conditions
u = 1.0;
u[nx/4..nx/2] = 2;

// create a temporary copy of 'u' to store the previous time step
var un = u;

// iterate for 'nt' time steps
for 1..nt {
  // swap arrays to prepare for next time step
  u <=> un;

  // update the solution over the interior of the domain in parallel
  forall i in indicesInner do
    u[i] = un[i] + alpha * (un[i-1] - 2 * un[i] + un[i+1]);
}

// print final results
const mean = (+ reduce u) / u.size,
      stdDev = sqrt((+ reduce (u - mean)**2) / u.size);
writeln("mean: ", mean, " stdDev: ", stdDev);
