# slicer-docker
Environment for executing the slicer trajectory inference algorithm.
Contains [functions](https://github.com/Stuartlab-UCSC/traj-converters) to convert the output to [common formats](https://github.com/Stuartlab-UCSC/traj-formats). The container is hosted on [docker hub](https://hub.docker.com/r/stuartlab/slicer/).


This readme documents the native analysis script of the container, which has a uniform api across all Stuart lab trajectory algorithm dockers. For other possible ways to use the container, see the documentation for the [monocle container](https://github.com/Stuartlab-UCSC/monocle-docker).

For instructions below I suggest creating a shared volume,
`./shared` in the directory you are launching the docker from. This 
volume acts as shared storage. The shared storage is the easiest way to
access input and output from the container.

You can run this command in bash (R needs to be installed) to create the test expression data, exp.tab:

`R -e 'set.seed(13);write.table(replicate(110, rnbinom(500, c(3, 10, 45, 100), .1)),file="shared/exp.tab",col.names=1:110, sep="\t")'`

## <a name="min"></a>Execute the container's native analysis script:
By convention our containers for TI algorithms have a`run_method` script in the image's $PATH. This allows a uniform interface for containers with different algorithms. The pattern is:

`docker run -v $(pwd):/data <container_name> run_method <tab delim exp file> <output.json>`

If the second positional argument is ommitted then a file called output_*method*.json will be produced. Hence, if your expression file 'exp.tab' is in your current working directory:

`docker run -v $(pwd):/data stuartlab/slicer run_method exp.tab`

will produce a 'output_slicer.json' file.
