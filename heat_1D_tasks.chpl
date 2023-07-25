/*
  A 1D finite difference heat/diffusion equation solver

  Computation is local (single compute node). Tasks are spawned
  manually using a `coforall` loop and synchronization is
  managed via `sync` variables in a global array of halo
  cells.

  Values of the `config const` variables can be modified in
  the command line (e.g., `./heat_1D --nt=100`)
*/

import Time.Timer;

// create a stopwatch to time kernel execution
var t = new Timer();

// declare configurable constants with default values
config const nx = 4096,     // number of grid points in x
             nt = 50,       // number of time steps
             alpha = 0.25;  // diffusion constant

// set the number of tasks equal to the number of available cores
const nTasks = here.maxTaskPar,
      npt = nx / nTasks;

// define a domain to describe the grid
const indices = {0..<nx},
      indicesInner = {1..<(nx-1)};

// define an array over the above domain
var u : [indices] real;

// set up initial conditions
u = 1.0;
u[nx/4..nx/2] = 2.0;

// define array of halo cells for each side of each task
var halos : [0..1, 0..<nTasks] sync real;
param LEFT = 0, RIGHT = 1;

// run the simulation across tasks
t.start();
coforall tid in 0..<nTasks do work(tid);

// print final results
const mean = (+ reduce u) / u.size,
      stdDev = sqrt((+ reduce (u - mean)**2) / u.size);
t.stop();
writeln("mean: ", mean, " stdDev: ", stdDev);
writeln("time: ", t.elapsed(), " (sec)");

proc work(tid: int) {
  // define region of the global array owned by this task
  const lo = tid * npt,
        hi = min((tid + 1) * npt, nx);

  const taskIndices = {lo..<hi},
        taskIndicesBuffered = taskIndices.expand(1),
        taskIndicesComp = taskIndices[indicesInner];

  // declare local array and load values from global array
  var uLocal1, uLocal2: [taskIndicesBuffered] real = 1.0;
  uLocal1[taskIndices] = u[taskIndices];
  uLocal2 = uLocal1;

  // iterate for 'nt' time steps
  for 1..nt {
    // write results from previous iteration to halo cells
    if tid != 0        then halos[RIGHT, tid-1].writeEF(uLocal2[taskIndices.low]);
    if tid != nTasks-1 then halos[LEFT, tid+1].writeEF(uLocal2[taskIndices.high]);

    // swap local arrays
    uLocal1 <=> uLocal2;

    // load values from halo cells into local array's borders
    if tid != 0        then uLocal1[taskIndicesBuffered.low] = halos[LEFT, tid].readFE();
    if tid != nTasks-1 then uLocal1[taskIndicesBuffered.high] = halos[RIGHT, tid].readFE();

    // compute the FD kernel in parallel
    foreach i in taskIndicesComp do
      uLocal2[i] = uLocal1[i] + alpha *
        (uLocal1[i-1] - 2 * uLocal1[i] + uLocal1[i+1]);
  }

  // store this task's results in global array
  uLocal1 <=> uLocal2;
  u[taskIndices] = uLocal1[taskIndices];
}
