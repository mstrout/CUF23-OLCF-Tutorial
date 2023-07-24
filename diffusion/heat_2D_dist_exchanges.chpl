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
       Collectives.barrier,
       Time.Timer;

// create a stopwatch to time kernel execution
var t = new Timer();

// compile with `-sRunCommDiag=true` to see comm diagnostics
use CommDiagnostics;
config param RunCommDiag = false;

// declare configurable constants with default values
config const nx = 4096,     // number of grid points in x
             ny = 4096,     // number of grid points in y
             nt = 50,       // number of time steps
             alpha = 0.25;  // diffusion constant

// define distributed domains and block-distributed array
const indices = {0..<nx, 0..<ny},
      indicesInner = indices.expand(-1),
      INDICES = Block.createDomain(indices);
var u: [INDICES] real;

// apply initial conditions
u = 1.0;
u[nx/4..nx/2, ny/4..ny/2] = 2.0;

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
  proc init(dGlobal: domain(2)) do this.d = dGlobal;
}

// set up an arrays of local arrays over same distribution as 'u.targetLocales'
var LOCALE_DOM = Block.createDomain(u.targetLocales().domain);
var uTaskLocal, unTaskLocal: [LOCALE_DOM] localArray;

proc main() {
  if RunCommDiag then startCommDiagnostics();

  // spawn one task for each locale
  t.start();
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
  t.stop();
  writeln("mean: ", mean, " stdDev: ", stdDev);
  writeln("time: ", t.elapsed(), " (sec)");
}

proc work(tidX: int, tidY: int) {
  // define domains to describe the indices owned by this task
  const myGlobalIndices = u.localSubdomain(here),
        localIndicesBuffered = myGlobalIndices.expand(1),
        localIndicesInner = myGlobalIndices[indicesInner];

  // initialize this tasks local arrays using indices from `Block` dist.
  uTaskLocal[tidX, tidY] = new localArray(localIndicesBuffered);
  unTaskLocal[tidX, tidY] = new localArray(localIndicesBuffered);

  // copy initial conditions from global array
  uTaskLocal[tidX, tidY].v = 1;
  uTaskLocal[tidX, tidY].v[myGlobalIndices] = u[myGlobalIndices];
  unTaskLocal[tidX, tidY].v = uTaskLocal[tidX, tidY].v;

  // define constants for indexing into edges of local array
  const N = myGlobalIndices.dim(1).low,
        S = myGlobalIndices.dim(1).high,
        E = myGlobalIndices.dim(0).high,
        W = myGlobalIndices.dim(0).low;

  b.barrier();

  // define constants for indexing into halo regions of neighboring arrays
  const NN = if tidY < tidYMax then unTaskLocal[tidX, tidY+1].d.dim(1).low else 0,
        SS = if tidY > 0       then unTaskLocal[tidX, tidY-1].d.dim(1).high else 0,
        WW = if tidX < tidXMax then unTaskLocal[tidX+1, tidY].d.dim(0).low else 0,
        EE = if tidX > 0       then unTaskLocal[tidX-1, tidY].d.dim(0).high else 0;

  // run FD computation
  for 1..nt {
    // store results from last iteration in neighboring task's halos
    if tidX > 0       then unTaskLocal[tidX-1, tidY].v[EE, ..] = unTaskLocal[tidX, tidY].v[W, ..];
    if tidX < tidXMax then unTaskLocal[tidX+1, tidY].v[WW, ..] = unTaskLocal[tidX, tidY].v[E, ..];
    if tidY > 0       then unTaskLocal[tidX, tidY-1].v[.., SS] = unTaskLocal[tidX, tidY].v[.., N];
    if tidY < tidYMax then unTaskLocal[tidX, tidY+1].v[.., NN] = unTaskLocal[tidX, tidY].v[.., S];

    // swap all local arrays
    b.barrier();
    uTaskLocal[tidX, tidY].v <=> unTaskLocal[tidX, tidY].v;

    // compute the FD kernel in parallel
    foreach (i, j) in localIndicesInner do
      unTaskLocal[tidX, tidY].v[i, j] = uTaskLocal[tidX, tidY].v[i, j] + alpha * (
          uTaskLocal[tidX, tidY].v[i-1, j] + uTaskLocal[tidX, tidY].v[i+1, j] +
          uTaskLocal[tidX, tidY].v[i, j-1] + uTaskLocal[tidX, tidY].v[i, j+1] -
          4 * uTaskLocal[tidX, tidY].v[i, j]
        );

    b.barrier();
  }

  // store results in global array
  uTaskLocal[tidX, tidY].v <=> unTaskLocal[tidX, tidY].v;
  u[myGlobalIndices] = uTaskLocal[tidX, tidY].v[myGlobalIndices];
}
