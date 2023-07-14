// declare configurable constants with default values
config const xLen = 2.0,    // length of the grid in x
             nx = 31,       // number of grid points in x
             nt = 50,       // number of time steps
             sigma = 0.25,  // stability parameter
             nu = 0.05;     // viscosity

// define non-configurable constants
const dx : real = xLen / (nx - 1),       // grid spacing in x
      dt : real = sigma * dx / nu,       // time step size
      nTasks = here.maxTaskPar,          // number of tasks
      npt = nx / nTasks;                 // number of grid points per task

// define a domain to describe the grid
const indices = {0..<nx};

// define an array over the above domain
var u : [indices] real;

// set up initial conditions
u = 1.0;
u[(0.5 / dx):int..<(1.0 / dx + 1):int] = 2;

// define array of ghost cells for each side of each task
var ghosts : [0..1, 0..<nTasks] sync real;
param LEFT = 0, RIGHT = 1;

proc work(tid: int) {
  // define region of the global array owned by this task
  const lo = tid * npt,
        hi = min((tid + 1) * npt, nx);

  const taskIndices = {lo..<hi},
        taskIndicesInner = taskIndices.expand(-1);

  // declare local array and load values from global array
  var uLocal1, uLocal2: [taskIndices] real;
  uLocal1 = u[taskIndices];

	// iterate for 'nt' time steps
  for 1..nt {
		// write results from previous iteration to ghost cells
		if tid != 0        then ghosts[RIGHT, tid-1].writeEF(uLocal2[taskIndicesInner.low]);
		if tid != nTasks-1 then ghosts[LEFT, tid+1].writeEF(uLocal2[taskIndicesInner.high]);

		// swap local arrays
		uLocal1 <=> uLocal2;

		// load values from ghost cells into local array's borders
		if tid != 0        then uLocal1[taskIndices.low] = ghosts[LEFT, tid].readFE();
		if tid != nTasks-1 then uLocal1[taskIndices.high] = ghosts[RIGHT, tid].readFE();

		// compute the FD kernel in parallel
		foreach i in taskIndicesInner do
			uLocal2[i] = uLocal1[i] + nu * dt / dx**2 *
				(uLocal1[i-1] - 2 * uLocal1[i] + uLocal1[i+1]);
		}

  // store this task's results in global array
	uLocal1 <=> uLocal2;
  u[taskIndices] = uLocal1;
}

// run the simulation across tasks
coforall tid in 0..<nTasks do work(tid);

// print final results
const mean = (+ reduce u) / u.size,
      stdDev = sqrt((+ reduce (u - mean)**2) / u.size);
writeln("mean: ", mean, " stdDev: ", stdDev);
