use StencilDist;

// declare configurable constants with default values
config const xLen = 2.0,    // length of the domain in x
             yLen = 2.0,    // length of the domain in y
             nx = 31,       // number of grid points in x
             ny = 31,       // number of grid points in y
             nt = 50,       // number of time steps
             sigma = 0.25,  // CFL condition
             nu = 0.05;     // viscosity

// define non-configurable constants
const dx: real = xLen / (nx - 1),       // grid spacing in x
      dy: real = yLen / (ny - 1),       // grid spacing in y
      dt: real = sigma * dx * dy / nu;  // time step size

// define a distributed 2D domain and subdomain to describe the grid and its interior
const indices = {0..<nx, 0..<ny} dmapped Stencil({1..<nx-1, 1..<ny-1}, fluff=(1,1)),
      indicesInner = indices[{1..<nx-1, 1..<ny-1}];

// define a distributed 2D array over the above domain
var u: [indices] real;

// set up initial conditions
u = 1.0;
u[
  (0.5 / dx):int..<(1.0 / dx + 1):int,
  (0.5 / dy):int..<(1.0 / dy + 1):int
] = 2;

// create a temporary copy of 'u' to store the previous time step
var un = u;

// iterate for 'nt' time steps
for 1..nt {
  // swap arrays to prepare for next time step
  u <=> un;

  // swap halo regions between neighboring tasks
  un.updateFluff();

  // compute the FD kernel in parallel
  forall (i, j) in indicesInner do
    u[i, j] = un[i, j] +
              nu * dt / dy**2 *
                (un[i-1, j] - 2 * un[i, j] + un[i+1, j]) +
              nu * dt / dx**2 *
                (un[i, j-1] - 2 * un[i, j] + un[i, j+1]);
}

// print final results
const mean = (+ reduce u) / u.size,
      stdDev = sqrt((+ reduce (u - mean)**2) / u.size);
writeln("mean: ", mean, " stdDev: ", stdDev);
