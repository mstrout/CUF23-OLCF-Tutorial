/*
  A distributed 2D finite-difference heat/diffusion equation solver

  Computation is executed over a 2D distributed array.
  The array distribution and sharing of "fluff"/"halo"
  regions are managed by the `Stencil` distribution. The
  `forall` loop manages task creation and synchronization
  across and within locales.

  Values of the `config const` variables can be modified in
  the command line (e.g., `./heat_2D_dist --nt=100`)
*/

import StencilDist.Stencil,
       Time.Timer;

// create a stopwatch to time kernel execution
var t = new Timer();

// declare configurable constants with default values
config const nx = 4096,     // number of grid points in x
             ny = 4096,     // number of grid points in y
             nt = 50,       // number of time steps
             alpha = 0.25;  // diffusion constant

// define a distributed 2D domain and subdomain to describe the grid and its interior
const indices = {0..<nx, 0..<ny},
      indicesInner = indices.expand(-1),
      INDICES = indices dmapped Stencil(indicesInner, fluff=(1,1)),
      INDICES_INNER = INDICES[{1..<nx-1, 1..<ny-1}];

// define a distributed 2D array over the above domain
var u: [INDICES] real;

// set up initial conditions
u = 1.0;
u[nx/4..nx/2, ny/4..ny/2] = 2.0;

// create a temporary copy of 'u' to store the previous time step
var un = u;

// iterate for 'nt' time steps
t.start();
for 1..nt {
  // swap arrays to prepare for next time step
  u <=> un;

  // swap halo regions between neighboring tasks
  un.updateFluff();

  // compute the FD kernel in parallel
  forall (i, j) in INDICES_INNER do
    u[i, j] = un[i, j] + alpha *
      (un[i, j-1] + un[i-1, j] + un[i+1, j] + un[i, j+1] - 4 * un[i, j]);
}

// print final results
const mean = (+ reduce u) / u.size,
      stdDev = sqrt((+ reduce (u - mean)**2) / u.size);
t.stop();
writeln("mean: ", mean, " stdDev: ", stdDev);
writeln("time: ", t.elapsed(), " (sec)");
