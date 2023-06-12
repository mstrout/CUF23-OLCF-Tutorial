use IO.FormattedIO;
//import surfPlot.printAndPlot;

// setup simulation parameters
config const xLen : real = 2,
             yLen : real = 2,
             nx = 31,
             ny = 31,
             nt = 50,
             sigma = 0.25,
             nu = 0.05;

const dx : real = xLen / (nx - 1),
      dy : real = yLen / (ny - 1),
      dt : real = sigma * dx * dy / nu;

config param createPlots = true;

// define 2D domain and subdomain
const dom = {0..<nx, 0..<ny};
const domInner : subdomain(dom) = dom.expand(-1);

// define initial conditions
var u : [dom] real = 1;
u[
    (0.5 / dx):int..<(1.0 / dx + 1):int,
    (0.5 / dy):int..<(1.0 / dy + 1):int
] = 2;

/* MMS, 1/4/23, may leave out of initial tutorial example
if createPlots then printAndPlot(
    u, (xLen, yLen), "results/diffusion_ics", "ICs", "density"
);
*/

// run finite difference computation
var un = u;
for n in 0..#nt {
    u <=> un;
    forall (i, j) in domInner {
        u[i, j] = un[i, j] +
                nu * dt / dy**2 *
                    (un[i-1, j] - 2 * un[i, j] + un[i+1, j]) +
                nu * dt / dx**2 *
                    (un[i, j-1] - 2 * un[i, j] + un[i, j+1]);
    }
}

/*
if createPlots then printAndPlot(
    u, (xLen, yLen), "results/diffusion_sol", "Sol", "density"
);
*/

// compute the standard deviation of 'u'
var mean = (+ reduce u) / u.size;
var std_dev = (+ reduce (u - mean)**2)**(0.5) / u.size;
writef("Final std(u): %.6dr\n", std_dev);

