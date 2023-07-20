/*
  A distributed 2D finite-difference heat/diffusion equation solver

  Computation is executed over a 2D distributed array.
  The array distribution is managed by the `Block` distribution.
  Tasks are spawned manually with a `coforall` loop and synchronization
  is done manually using a `barrier`. Halo regions are shared across
  locales manually via direct assignment between locally owned arrays.

  Similar to `heat_2D_dist_exchanges.chpl`, except local array definition
  and halo exchanges are facilitated by `localArrayPair` type.

  Values of the `config const` variables can be modified in
  the command line (e.g., `./heat_2D_exchanges --nt=100`)
*/

import BlockDist.Block,
       Collectives.barrier;

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

// buffer edge enum: North, East, South, West
enum Edge { N, E, S, W }

// a type to store a tasks local arrays and facilitate sharing of
//  halo regions between neighboring tasks
record localArrayPair {
  var dummy = false;

  // the set of global indices that this locale owns
  //  used to index into the global array
  var globalIndices: domain(2);

  // indices over which local arrays are defined
  //  same as 'globalIndices' set with a halo regions along edges
  var indices: domain(2);

  // indices over which the kernel is executed
  //  can differ from 'globalIndices' for locales with
  //  indices on the border of the global domain
  var compIndices: domain(2);

  var u: [indices] real;
  var un: [indices] real;

  // default initializer
  proc init(dummy = false) {
    this.dummy = dummy;
    this.globalIndices = {0..0, 0..0};
    this.indices = {0..0, 0..0};
    this.compIndices = {0..0, 0..0};
  }

  proc init(myGlobalIndices: domain(2), globalInnerAll: domain(2)) {
    this.dummy = false;
    this.globalIndices = myGlobalIndices;
    this.indices = myGlobalIndices.expand(1);
    this.compIndices = myGlobalIndices[globalInnerAll];
  }

  proc ref copyInitialConditions(const ref uGlobal: [] real) {
    this.u = 1.0;
    this.u[this.globalIndices] = uGlobal[globalIndices];
    this.un = this.u;
  }

  proc fillHalo(edge: Edge, const ref values: [] real) {
    if !dummy {
      select edge {
        when Edge.N do this.un[.., this.indices.dim(1).low] = values;
        when Edge.S do this.un[.., this.indices.dim(1).high] = values;
        when Edge.E do this.un[this.indices.dim(0).high, ..] = values;
        when Edge.W do this.un[this.indices.dim(0).low, ..] = values;
      }
    }
  }

  // `[ ]` access
  proc this(edge: Edge) {
    select edge {
      when Edge.N do return this.un[.., this.indices.dim(1).low+1];
      when Edge.S do return this.un[.., this.indices.dim(1).high-1];
      when Edge.E do return this.un[this.indices.dim(0).high-1, ..];
      when Edge.W do return this.un[this.indices.dim(0).low+1, ..];
    }
    return this.un[this.indices.dim(0).low, ..]; // never actually returned
  }
}

// set up an array of local arrays over same distribution as 'u.targetLocales'
var LOCALE_DOM = Block.createDomain(u.targetLocales().domain),
    uTaskLocal : [LOCALE_DOM] localArrayPair;

proc main() {
  if RunCommDiag then startCommDiagnostics();

  // spawn one task for each locale
  coforall (loc, (tidX, tidY)) in zip(u.targetLocales(), LOCALE_DOM) {
    // run initialization and computation on the task for this locale
    on loc {
      // initialize local array pairs
      uTaskLocal[tidX, tidY] = new localArrayPair(
        u.localSubdomain(here), // indices owned by this locale from `Block` dist
        indicesInner            // global "compIndices"
      );
      uTaskLocal[tidX, tidY].copyInitialConditions(u);

      // synchronize across tasks
      b.barrier();

      // run the portion of the FD computation owned by this task
      work(tidX, tidY);
    }
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
  // array pair whose `fillHalo` call is a no-op
  var uDummy = new localArrayPair(dummy=true);

  // get references to the local array pairs for this task and its neighbors
  ref uLocal = uTaskLocal[tidX, tidY],
      uWest  = if tidX > 0       then uTaskLocal[tidX-1, tidY] else uDummy,
      uEast  = if tidX < tidXMax then uTaskLocal[tidX+1, tidY] else uDummy,
      uNorth = if tidY > 0       then uTaskLocal[tidX, tidY-1] else uDummy,
      uSouth = if tidY < tidYMax then uTaskLocal[tidX, tidY+1] else uDummy;

  // run FD computation
  for 1..nt {
    // store results from last iteration in neighboring task's halos
    uWest.fillHalo(Edge.E, uLocal[Edge.W]);
    uEast.fillHalo(Edge.W, uLocal[Edge.E]);
    uNorth.fillHalo(Edge.S, uLocal[Edge.N]);
    uSouth.fillHalo(Edge.N, uLocal[Edge.S]);

    // swap local arrays
    b.barrier();
    uLocal.u <=> uLocal.un;

    // compute the FD kernel in parallel
    foreach (i, j) in uLocal.compIndices do
      uLocal.un[i, j] = uLocal.u[i, j] + alpha * (
        uLocal.u[i-1, j] + uLocal.u[i+1, j] +
        uLocal.u[i, j-1] + uLocal.u[i, j+1] -
        4 * uLocal.u[i, j]
      );

    b.barrier();
  }

  // store results in global array
  uLocal.u <=> uLocal.un;
  u[uLocal.globalIndices] = uLocal.u[uLocal.globalIndices];
}
