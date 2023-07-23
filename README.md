# CUF23-OLCF-Tutorial Chapel Examples

This README.md describes how to use Chapel during and after the tutorial.
It then lists the example codes that are described in the tutorial slides
and extra ones included in this repository.  The tutorial slides themselves
are also included in this repository.

If you have any questions about this tutorial material, please email
michelle.strout@hpe.com with the subject heading "CUF23-OLCF-Tutorial question".

## During the tutorial

wget https://go.lbl.gov/cuf23.tar.gz

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
See the commands in the provided FIXME.bash script for setup on Frontier
and change them so are using your account and the correct queue.
Here are examples of those detailed commands.

module load ums/default ums014/default
module load chapel

git clone git@github.com:mstrout/CUF23-OLCF-Tutorial.git
cd CUF23-OLCF-Tutorial
chpl hello6-taskpar-dist.chpl

export CHPL_LAUNCHER_ACCOUNT=CSC296
export CHPL_LAUNCHER_WALLTIME=00:10:00
./hello6-taskpar-dist -nl 2

### Perlmutter

See the commands in the provided FIXME.bash script for setup on Perlmutter
and change them so are using your account and the correct queue.
Here are examples of those detailed commands.

FIXME

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
FIXME: write and test including basics-distarr.chpl

- `spmd/toy.chpl` and `spmd/toy-SPMD.chpl` - examples showing typical chapel and chapel
  using SPMD
FIXME: write and test

- `parfilekmer.chpl` - parallel and distributed processing of files using a 1D distributed
  array of filenames

- `onClause.chpl` - shows how computation on a locale X can access a variable on a locale Y
FIXME: write and test

- FIXME Scott's code, grab from testing hierarchy

- `gpuExample.chpl` - illustrates how to run code either on the current locale in
  parallel or across GPUs if available

- `stream-ep.chpl` - shows how to run in parallel on the CPU and all the GPUs available

## Additional examples

- `writelnExamples.chpl` - you can write out almost any variable in Chapel with `writeln`

- `diffusion/*.chpl` - there are some 1D and other 2D heat diffusion examples in this
   directory.
