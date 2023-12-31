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
// you might want to read the README.md in this directory
source SOURCE_ME.sh
cd chapel
```

To compile and run any of the example programs, e.g. hello.chpl, do the following command:
```
make run-hello
```

To run the `image_analysis/main.chpl` example.
```
  cd image_analysis
  chpl main.chpl --fast
  ./main --in_name=banda_ai --map_type=benthic --window_size=100000
```

Do some scaling experiments:
```
for ((i=8;i<=64;i*=2)); do echo @$i; env CHPL_RT_NUM_THREADS_PER_LOCALE=$i ./heat_2D -nl 1; done |& tee ll
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


## Examples covered in slides

20-minute intro
- `hello-dist-node-name.chpl` - parallel and distributed 'Hello World' that prints out locale names

- `heat_2D.chpl` - shared memory parallel version that runs in locale 0

- `heat_2D_dist.chpl` - parallel and distributed version that is the same as

  `heat_2D.chpl` but with distributed arrays

- `heat_2D_dist_exchanges.chpl` - parallel and distributed version that copies from local subarrays into neighbors' halos


90-minute tutorial
- `kmer.chpl` - k-mer counting from bioinformatics, shows how to use map/dictionary

- `hellopar.chpl` - another parallel and distributed 'Hello World', short version of
  hello6-taskpar-dist.chpl

- `basics-*.chpl` - some basic illustrations of parallelism and/or data locality
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

- `heat*.chpl` - there are some 1D and other 2D heat diffusion examples in this
   directory.  The `heat_1D_tasks.chpl` example has an example of using `sync` variables.
