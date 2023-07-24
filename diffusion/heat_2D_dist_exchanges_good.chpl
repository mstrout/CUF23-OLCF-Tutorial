/*
  A distributed 2D finite-difference heat/diffusion equation solver

  Computation is executed over a 2D distributed array.
  The array distribution is managed by the `Block` distribution.
  Tasks are spawned manually with a `coforall` loop and synchronization
  is done manually using a `barrier`. Halo regions are shared across
  locales manually via direct assignment between arrays.

  Values of the `config const` variables can be modified in
  the command line (e.g., `./heat_2D_exchanges --nt=100`)
*/

import BlockDist.Block,
       Collectives.barrier,
       Time.Timer;

// create a stopwatch to time kernel execution
var t = new Timer();

use CommDiagnostics;
config param runCommDiag = false;

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

// class to facilitate sharing of local edges between locales
enum Edge { N, E, S, W }
class LocalArrayPair {
  // the set of global indices that this locale owns
  //  used to index into the global array
  var globalIndices: domain(2) = {0..0, 0..0};

  // indices over which local arrays are defined
  //  same as 'globalIndices' set with a halo region along the edge
  var indices: domain(2) = {0..0, 0..0};

  // indices over which the kernel is executed
  //  can differ from 'globalIndices' for locales along the border of the global domain
  var compIndices: domain(2) = {0..0, 0..0};

  var u: [indices] real;
  var un: [indices] real;

  proc init() {
    this.globalIndices = {0..0, 0..0};
    this.indices = {0..0, 0..0};
    this.compIndices = {0..0, 0..0};
  }

  proc init(myGlobalIndices: domain(2), globalInnerAll: domain(2)) {
    this.globalIndices = myGlobalIndices;
    this.indices = myGlobalIndices.expand(1);
    this.compIndices = myGlobalIndices[globalInnerAll];
  }

  proc ref copyInitialConditions(const ref u: [] real) {
    this.u = 1.0;
    this.u[this.globalIndices] = u[globalIndices];
    this.un = this.u;
  }

   proc fillBuffer(edge: Edge, values: [] real) {
    select edge {
      when Edge.N do this.un[.., this.indices.dim(1).low] = values;
      when Edge.S do this.un[.., this.indices.dim(1).high] = values;
      when Edge.E do this.un[this.indices.dim(0).high, ..] = values;
      when Edge.W do this.un[this.indices.dim(0).low, ..] = values;
    }
  }

  proc getEdge(edge: Edge) {
    select edge {
      when Edge.N do return this.un[.., this.indices.dim(1).low+1];
      when Edge.S do return this.un[.., this.indices.dim(1).high-1];
      when Edge.E do return this.un[this.indices.dim(0).high-1, ..];
      when Edge.W do return this.un[this.indices.dim(0).low+1, ..];
    }
    return this.un[this.indices.dim(0).low, ..]; // never actually returned
  }

  proc swap() do this.u <=> this.un;
}

// set up an array of local array pairs over same distribution as 'u.targetLocales'
var TL_DOM = Block.createDomain(u.targetLocales().domain);
var localArrays = [d in TL_DOM] new LocalArrayPair();

proc main() {
  if runCommDiag then startCommDiagnostics();

  // execute the FD computation with one task per locale
  t.start();
  coforall (loc, (tidX, tidY)) in zip(u.targetLocales(), u.targetLocales().domain) do on loc {
    // initialize local arrays owned by this task
    localArrays[tidX, tidY] = new LocalArrayPair(u.localSubdomain(here), indicesInner);
    localArrays[tidX, tidY].copyInitialConditions(u);

    // synchronize across tasks
    b.barrier();

    // run the portion of the FD computation owned by this task
    work(tidX, tidY);
  }

  if runCommDiag {
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
  // get a pointer to this task's local array-pair
  var uLocal: borrowed LocalArrayPair = localArrays[tidX, tidY];

  // run FD computation
  for 1..nt {
    // store results from last iteration in neighboring task's halos
    if tidX > 0       then localArrays[tidX-1, tidY].fillBuffer(Edge.E, uLocal.getEdge(Edge.W));
    if tidX < tidXMax then localArrays[tidX+1, tidY].fillBuffer(Edge.W, uLocal.getEdge(Edge.E));
    if tidY > 0       then localArrays[tidX, tidY-1].fillBuffer(Edge.S, uLocal.getEdge(Edge.N));
    if tidY < tidYMax then localArrays[tidX, tidY+1].fillBuffer(Edge.N, uLocal.getEdge(Edge.S));

    // swap local 'u' and 'un' arrays
    b.barrier();
    uLocal.swap();

    // compute the FD kernel in parallel
    foreach (i, j) in uLocal.compIndices do
      uLocal.un[i, j] = uLocal.u[i, j] + alpha * (
          uLocal.u[i-1, j] + uLocal.u[i+1, j] +
          uLocal.u[i, j-1] + uLocal.u[i, j+1] +
          4 * uLocal.u[i, j]
        );

    b.barrier();
  }

  // store results in global array
  uLocal.swap();
  u[uLocal.globalIndices] = uLocal.u[uLocal.globalIndices];
}
