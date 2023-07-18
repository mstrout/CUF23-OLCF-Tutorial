/*
  A distributed 2D finite-difference heat/diffusion equation solver

  Computation is executed over a 2D distributed array.
  The array distribution is managed by the `Block` distribution.
  Tasks are spawned manually with a `coforall` loop and synchronization
  is done manually using a `barrier`. Halo regions are shared across
  locales manually via halo/buffer arrays.

  Values of the `config const` variables can be modified in
  the command line (e.g., `./heat_2D_buffers --nt=100`)
*/

import BlockDist.Block,
       Collectives.barrier;

use CommDiagnostics;
config param runCommDiag = false;

// declare configurable constants with default values
config const xLen = 2.0,    // length of the grid in x
             yLen = 2.0,    // length of the grid in y
             nx = 31,       // number of grid points in x
             ny = 31,       // number of grid points in y
             nt = 50,       // number of time steps
             sigma = 0.25,  // stability parameter
             nu = 0.05;     // viscosity

// define non-configurable constants
const dx : real = xLen / (nx - 1),       // grid spacing in x
      dy : real = yLen / (ny - 1),       // grid spacing in y
      dt : real = sigma * dx * dy / nu;  // time step size

// define distributed domains and block-distributed array
const indices = {0..<nx, 0..<ny},
      indicesInner = indices.expand(-1),
      INDICES = Block.createDomain(indices);
var u: [INDICES] real;

// apply initial conditions
u = 1.0;
u[
  (0.5 / dx):int..<(1.0 / dx + 1):int,
  (0.5 / dy):int..<(1.0 / dy + 1):int
] = 2;

// type for creating a "skyline" array of halos
record haloArray {
  var d: domain(1);
  var v: [d] real;

  proc init() do this.d = {0..0};
  proc init(r: range(int, boundKind.both, strideKind.one)) do
    this.d = {r};
}

// set up array of halos over same distribution as 'u.targetLocales'
var TL_DOM = Block.createDomain(u.targetLocales().domain);
var haloArrays: [TL_DOM] [0..<4] haloArray;

// numbers for indexing into local array edges
param N = 0, S = 1, E = 2, W = 3;

// number of tasks per dimension based on Block distribution's decomposition
const tidXMax = u.targetLocales().dim(0).high,
      tidYMax = u.targetLocales().dim(1).high;

// barrier for one task per locale
var b = new barrier(numLocales);

proc main() {
  if runCommDiag then startCommDiagnostics();

  // spawn one task for each locale
  coforall (loc, (tidX, tidY)) in zip(u.targetLocales(), u.targetLocales().domain) {
    // run initialization and computation on the task for this locale
    on loc {
      // initialize halo arrays
      haloArrays[tidX, tidY][N] = new haloArray(u.localSubdomain().dim(1).expand(1));
      haloArrays[tidX, tidY][S] = new haloArray(u.localSubdomain().dim(1).expand(1));
      haloArrays[tidX, tidY][E] = new haloArray(u.localSubdomain().dim(0).expand(1));
      haloArrays[tidX, tidY][W] = new haloArray(u.localSubdomain().dim(0).expand(1));

      // synchronize across tasks
      b.barrier();

      // run the portion of the FD computation owned by this task
      work(tidX, tidY);
    }
  }

  if runCommDiag {
    stopCommDiagnostics();
    printCommDiagnosticsTable();
  }

  // print final results
  const mean = (+ reduce u) / u.size,
        stdDev = sqrt((+ reduce (u - mean)**2) / u.size);
  writeln("mean: ", mean, " stdDev: ", stdDev);
}

proc work(tidX: int, tidY: int) {
  // define domains to describe the indices owned by this task
  const localIndices = u.localSubdomain(here),
        localIndicesBuffered = localIndices.expand(1),
        localIndicesInner = localIndices[indicesInner];

	// declare two local arrays for computation on this locale
  var uLocal1, uLocal2: [localIndicesBuffered] real = 1;

  // populate local arrays with initial conditions from global array
  uLocal1[localIndices] = u[localIndices];
  uLocal2 = uLocal1;

  // define constants for indexing into edges of local array
  const WW = localIndicesBuffered.dim(1).low,
        EE = localIndicesBuffered.dim(1).high,
        SS = localIndicesBuffered.dim(0).high,
        NN = localIndicesBuffered.dim(0).low;

  // iterate for 'nt' time steps
  for 1..nt {
    // store results from last iteration in neighboring task's halos
    if tidY > 0       then haloArrays[tidX, tidY-1][E].v = uLocal2[.., WW+1];
    if tidY < tidYMax then haloArrays[tidX, tidY+1][W].v = uLocal2[.., EE-1];
    if tidX > 0       then haloArrays[tidX-1, tidY][S].v = uLocal2[NN+1, ..];
    if tidX < tidXMax then haloArrays[tidX+1, tidY][N].v = uLocal2[SS-1, ..];

    // swap local arrays
    b.barrier();
    uLocal1 <=> uLocal2;

    // populate edges of local array from halo arrays
    if tidY > 0       then uLocal1[.., WW] = haloArrays[tidX, tidY][W].v;
    if tidY < tidYMax then uLocal1[.., EE] = haloArrays[tidX, tidY][E].v;
    if tidX > 0       then uLocal1[NN, ..] = haloArrays[tidX, tidY][N].v;
    if tidX < tidXMax then uLocal1[SS, ..] = haloArrays[tidX, tidY][S].v;

    // compute the FD kernel in parallel
    forall (i, j) in localIndicesInner do
      uLocal2[i, j] = uLocal1[i, j] +
        nu * dt / dy**2 *
          (uLocal1[i-1, j] - 2 * uLocal1[i, j] + uLocal1[i+1, j]) +
        nu * dt / dx**2 *
          (uLocal1[i, j-1] - 2 * uLocal1[i, j] + uLocal1[i, j+1]);

    b.barrier();
  }

  // store results in global array
  uLocal1 <=> uLocal2;
  u[localIndices] = uLocal1[localIndices];
}
