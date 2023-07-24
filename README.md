# CUF23-OLCF-Tutorial Chapel Examples

This README.md describes how to use Chapel during and after the tutorial.
It then lists the example codes that are described in the tutorial slides
and extra ones included in this repository.  The tutorial slides themselves
are also included in this repository.

If you have any questions about this tutorial material, please email
michelle.strout@hpe.com with the subject heading "CUF23-OLCF-Tutorial question".

## During the tutorial

```
wget https://go.lbl.gov/cuf23.tar.gz
tar xzf cuf23.tar.gz
cd cuf23
FIXME: what next?
```

## After the tutorial

After the tutorial, the AWS instances and guest Perlmutter accounts will
go away.  However, you can still work with Chapel by trying examples
through a web browser, with Docker, asking your HPC consultants to install
Chapel on your local HPC resources, or by following the instructions
on the Chapel documentation webpages for installing Chapel yourself.

### Webbrowser
[Attempt this online](https://ato.pxeger.com/run?1=m70sOSOxIDVnwYKlpSVpuhY7y4syS1Jz8jSUPFJzcvJ1FMrzi3JSFJU0rSHyUGUw5QA)

### Docker
Using a container on your laptop
- First, install docker for your machine and then start it up
- Then, the use the below commands
```
 docker pull docker.io/chapel/chapel-gasnet    # takes about 5 minutes
 docker run --rm -v "$PWD":/myapp -w /myapp chapel/chapel-gasnet chpl hello.chpl
 docker run --rm -v "$PWD":/myapp -w /myapp chapel/chapel-gasnet ./hello -nl 1
```

### Frontier
See the commands in the provided `SOURCE_ME.sh` script for setup on Frontier
and change them so are using your account and the correct queue.
Here are examples of those detailed commands.

```
module load ums/default ums014/default
module load chapel

git clone git@github.com:mstrout/CUF23-OLCF-Tutorial.git
cd CUF23-OLCF-Tutorial
chpl hello6-taskpar-dist.chpl

export CHPL_LAUNCHER_ACCOUNT=CSC296
export CHPL_LAUNCHER_WALLTIME=00:10:00
./hello6-taskpar-dist -nl 2
```

### Perlmutter

See the commands in the provided FIXME.bash script for setup on Perlmutter
and change them so are using your account and the correct queue.
Here are examples of those detailed commands.

```
module unload $(module --terse list 2>&1 | grep PrgEnv-)
module load PrgEnv-gnu
module load cray-pmi
module load chapel

# Work around cxi provider bugs that limit memory registration
# AND using a smaller heap size to reduce startup time for examples
export CHPL_RT_MAX_HEAP_SIZE=16GB
export CHPL_LAUNCHER_MEM=unset

# to avoid doing a manual salloc, this isn't working
export SLURM_QOS=interactive # do I actually need this?
export SLURM_CONSTRAINT=cpu
export SLURM_ACCOUNT=<account>
export CHPL_LAUNCHER_WALLTIME=00:10:00
```
Figure out what account you can charge to by using the iris command.

## Examples covered in slides

20-minute intro
- `hello-dist-node-name.chpl` - parallel and distributed 'Hello World' that prints out locale names

- `diffusion/heat_2D.chpl` - shared memory parallel version that runs in locale 0

- `diffusion/heat_2D_dist.chpl` - parallel and distributed version that is the same as

  `heat_2D.chpl` but with distributed arrays

- `heat_2D_dist_exchanges.chpl` - parallel and distributed version that copies from local subarrays into neighbors' halos


90-minute tutorial
- `kmer.chpl` - k-mer counting from bioinformatics, shows how to use map/dictionary

- `hellopar.chpl` - another parallel and distributed 'Hello World', short version of
  hello6-taskpar-dist.chpl

- `basics/basics-*.chpl` - some basic illustrations of parallelism and/or data locality
  *NOTE*: basics-on.chpl assume two or more locales

- `parfilekmer.chpl` - parallel and distributed processing of files using a 1D distributed
  array of filenames

- `image_analysis/main.chpl` - coral reef diversity example by Scott Bachman at NCAR
  shows how to read a single file in parallel
```
  chpl main.chpl --fast
  ./main --in_name=banda_ai --map_type=benthic --window_size=100000
```

- `gpuExample.chpl` - illustrates how to run code either on the current locale in
  parallel or across GPUs if available

- `stream-ep.chpl` - shows how to run in parallel on the CPU and all the GPUs available

## Additional examples

- `writelnExamples.chpl` - you can write out almost any variable in Chapel with `writeln`

- `diffusion/*.chpl` - there are some 1D and other 2D heat diffusion examples in this
   directory.
