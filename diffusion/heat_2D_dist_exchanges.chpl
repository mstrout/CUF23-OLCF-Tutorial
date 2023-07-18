/*
  A distributed 2D finite-difference heat/diffusion equation solver

  Computation is executed over a 2D distributed array.
  The array distribution is managed by the `Block` distribution.
  Tasks are spawned manually with a `coforall` loop and synchronization
  is done manually using a `barrier`. Halo regions are shared across
  locales manually via direct assignment between locally owned arrays.

  Values of the `config const` variables can be modified in
  the command line (e.g., `./heat_2D_exchanges --nt=100`)
*/

import BlockDist.Block,
       Collectives.barrier;

// compile with `-sRunCommDiag=true` to see comm diagnostics
use CommDiagnostics;
config param RunCommDiag = false;

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

// number of tasks per dimension based on Block distributions decomposition
const tidXMax = u.targetLocales().dim(0).high,
      tidYMax = u.targetLocales().dim(1).high;

// barrier for one task per locale
var b = new barrier(u.targetLocales().size);

// a type for creating a "skyline" array of local arrays
record localArray {
  var d: domain(2);
  var v: [d] real;

  proc init() do this.d = {0..0, 0..0};
  proc init(d: domain(2)) do this.d = d;
}

// set up an arrays of local arrays over same distribution as 'u.targetLocales'
var LOCALE_DOM = Block.createDomain(u.targetLocales().domain);
var uTaskLocal, unTaskLocal: [LOCALE_DOM] localArray;

proc main() {
  if RunCommDiag then startCommDiagnostics();

  // spawn one task for each locale
  coforall (loc, (tidX, tidY)) in zip(u.targetLocales(), LOCALE_DOM) {
    // run initialization and computation on the task for this locale
    on loc do work(tidX, tidY);
  }

  if RunCommDiag {
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
  const myGlobalIndices = u.localSubdomain(here),
        localIndicesBuffered = myGlobalIndices.expand(1),
        localIndicesInner = myGlobalIndices[indicesInner];

  // initialize this tasks local arrays using indices from `Block` dist.
  uTaskLocal[tidX, tidY] = new localArray(localIndicesBuffered);
  unTaskLocal[tidX, tidY] = new localArray(localIndicesBuffered);

  // get a reference to this task's local arrays
  ref uLocal = uTaskLocal[tidX, tidY].v;
  ref unLocal = unTaskLocal[tidX, tidY].v;

  // copy initial conditions from global array
  uLocal = 1;
  uLocal[myGlobalIndices] = u[myGlobalIndices];
  unLocal = uLocal;

  // define constants for indexing into edges of local array
  const NN = localIndicesBuffered.dim(1).low,
        SS = localIndicesBuffered.dim(1).high,
        EE = localIndicesBuffered.dim(0).high,
        WW = localIndicesBuffered.dim(0).low;

  b.barrier();

  // run FD computation
  for 1..nt {
    // store results from last iteration in neighboring task's halos
    if tidX > 0       then unTaskLocal[tidX-1, tidY].v[EE, ..] = unLocal[WW+1, ..];
    if tidX < tidXMax then unTaskLocal[tidX+1, tidY].v[WW, ..] = unLocal[EE-1, ..];
    if tidY > 0       then unTaskLocal[tidX, tidY-1].v[.., SS] = unLocal[.., NN+1];
    if tidY < tidYMax then unTaskLocal[tidX, tidY+1].v[.., NN] = unLocal[.., SS-1];

    // swap all local arrays
    b.barrier();
    uLocal <=> unLocal;

    // compute the FD kernel in parallel
    foreach (i, j) in localIndicesInner do
      unLocal[i, j] = uLocal[i, j] +
              nu * dt / dy**2 *
                (uLocal[i-1, j] - 2 * uLocal[i, j] + uLocal[i+1, j]) +
              nu * dt / dx**2 *
                (uLocal[i, j-1] - 2 * uLocal[i, j] + uLocal[i, j+1]);

    b.barrier();
  }

  // store results in global array
  uLocal <=> unLocal;
  u[myGlobalIndices] = uLocal[myGlobalIndices];
}
