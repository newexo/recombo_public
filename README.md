# RECOMBO 

## Synopsis
RECOMBO is the working title of the Topological Molecular Biology Recombination Tool Suite. The suite is made up of a number of different applications that share a common code base are designed to work well together. As a basic overview, these applications cover everything from generating and performing topological reconnection on SAP conformations of fixed topology (MMC), removing reidermeister crossings in knot shadows (xinger), computing homfly polynomials from extended gauss codes (homfly), and identifying homfly polynomials with a topology (idknot). RECOMBO also includes an integrated unit testing suite with ok-ish test coverage, and a thoroughly experimental fugacity parameter finder (zAnalyzer). Each application in this suite is described in detail below.

## Installation
This software suite is intended for use on Linux and OSX. Install by following these steps:
1. Download and extract this repository.
2. In a shell command line, navigate inside the repository's `src/` subdirectory. 
3. execute `make all`.
4. Assuming compilation is successful and all unit tests pass, navigate to the `/src/bin` subdirectory to find the executables. 
5. If compilation is not successful, examine the error messages and ensure your compiler and environment are treating the code appropriately.

## MMC
MMC is the largest and most complex piece of software in the RECOMBO suite. It uses an extension of the BFACF algorithm (CMC-BFACF) [Stolz et. al., 2017,Szafron,Orlandini] to generate and sample self-avoiding polygons. Because of its complexity, there are a number of options and run-modes avaliable to the user to best meet your simulation needs. 

#### Usage
For a more detailed explanation of the theory behind Composite Markov Chain Monte Carlo, I refer the reader to [Stolz PhD Thesis]. Here is described a number of the important arguments with examples. Determining which arguments are specifically best for their simulation is left to the reader, though suggestions appear below. (wording)

#### Argument List

`-zmin` the smallest z value set to one of the initial markov chains

`-zmax` the largest z value set to one of the initial markov chains

`-sr` the user-defined minimum acceptaple swap ratio that must be satisfied between adjacent markov chains. If the swap-ratio is not satisfied, a new Markov chain will be added with a z value half-way between them on a linear curve. 

`-s` the swap interval, or the number of BFACF moves between swap attempts in each chain

`-c` the sampling interval, or the number of BFACF moves between samples (or attempted samples)

`-m` initial number of chains initialized on the interval defined by the `-zmin` and `-zmax` arguments. 

The sampling mode dictates how BFACF samples conformations. Basic options:

`-mode a` Analyze only mode. Perform warmup and calibration, measuring and outputting statistics for each chain. 

`-mode s` Sample mode. Sample conformations in at an interval specified by the `-c` option into a binary format file.

`-bfs n` Block-file size. Specifies the number of conformations to save per file, `n`. Useful for splitting up the dataset. With the `-bfs` option unused or if `n` is 1, each conformation will be output in its own CUBE binary formatted `.b` file. NOTE: The binary representation of the CUBE format should be thought of as a linked list. Retrieving a conformation at a particular index requires iterating through all the preceding conformations in the file. Therefore, it is advisable to think about how you will want to organize and analyze the output conformations when specifing the `-bfs` option.

`-mode b` Combination of `a` and `s` mode.

Advanced operating modes work in conjunction with recombo to perform reconnection in mmc:

`-mode r` Reconnection mode. Samples like `s` mode but saves a post reconnection version of each sampled conformation in a separate file. Requires reconnection criteria to also be specified.

`-mode f` Filter sample mode. Samples like `r` mode but saves only conformations that have sites satisfying the defined reconnection criteria.

`-minarc n` Mandatory when in one of the two advanced operating modes. Defines the minimum arclength to identify a reconnection site as the integer `n`.

`-maxarc n` Mandatory when in one of the two advanced operating modes. Defines the maximum arclength to identify a reconnection site as the integer `n`.

`-targetlength n` Optional when in one of the two advanced operating modes. For knots, the conformation must be of length `n` within a defined arclength tolerance. For two-component links, the components must add to `n` within a defined arclength tolerance.

`-recomboParas < >` is used in place of recomboFromFile.

`-recomboParas ds` Optional when in f mode. Do standard (on-lattice) recombination in direct repeat

`-recomboParas da` Optional when in f mode. Do writhe_based virtual (off-lattice) recombination in direct repeat

`-recomboParas dp` Optional when in f mode. Do positive virtual (off-lattice) recombination in direct repeat

`-recomboParas dn` Optional when in f mode. Do negative virtual (off-lattice) recombination in direct repeat 

`-recomboParas dsa` Optional when in f mode. Do a mix of standard and writhe-based recombination in direct repeat

`-recomboParas dsp` Optional when in f mode. Do a mix of standard and positive virtual recombination in direct repeat

`-recomboParas dsn` Optional when in f mode. Do a mix of standard and negative virtual recombination in direct repeat

`-recomboParas is` Optional when in f mode. Do standard recombination in inverted repeat

`-recomboParas ip` Optional when in f mode. Do positive virtual recombination in inverted repeat  (Not implemented)

`-recomboParas in` Optioanl when in f mode. Do positive virtual recombination in inverted repeat  (Not implemented)

`-recomboParas ia` Optional when in f mode. Do writhe-based virtual recombination in inverted repeat  (Not implemented)

`-recomboParas isp` Optioanl when in f mode. Do a mix of standard and positive virtual recombination in inverted repeat  (Not implemented)

`-recomboParas isn` Optional when in f mode. Do a mix of standard and negative virtual recombination in inverted repeat  (Not implemented)

`-recomboParas isa` Optional when in f mode. Do a mix of standard and writhe-based virtual recombination in inverted repeat  (Not implemented)

##### Optional Arguments
`-q` (Optional) the q value used for all markov chains in the mmc process. If unspecified, defaults to q=1.

`-seed` (Optional) Manually specify a seed for the random number generator. If unspecified, uses the system clock and records the seed to stdout (console output).

`+s` (Optional) Supresses progress output. Should be used when running long term jobs on a server that are not being directly monitored during the run.

`-info`(Optional) Outputs information needed to perform batch-mean analysis as described in George Fishman's "Monte Carlo: Concepts, Algorithms, and Applications" into a new file with the .info extension. 

##### Halting Conditions
At least one of these must be specified for the program to terminate. If more than one is specified, mmc will sample conformations until any halting condition is met. 

`-n` number of conformations to sample

`-t` integer number of hours to run before safe termination. 

#### Example Usage (revise)
`mmc initial/3_1 3_1 -zmin 0.20466 -zmax 0.20966 -q 1 -sr .8 -s 5 -n 800000 -c 20000 -m 5 -w 1000000 -mode m -minarc 57 -maxarc 63 -targetlength 120 -seed 42 +s > 3_1_log.txt &`

## recomboFromFile (Not working now. It is integrated to MMC by -recomboParas < >)
recomboFromFile takes an input file or stream of CUBE binary formatted polygons, and attempts to perform reconnection on each one according to user-specified reconnection criteria. It is capable of performing reconnection in direct ~~and inverted~~ repeat. 

## homfly
Homfly takes as an input a file of extended gauss codes (.egc), and outputs the HOMFLY-PT polynomial [x] for each one. 

#### Usage
`homfly` for interactive usage. If you are a mathematician looking to 

`homfly input_file output_file` 

`homfly -- |` or ` homfly -- > output_file` for use in a stream.

## xinger


## idknot
Idknot takes a input file or stream of homfly polynomials and associates each one with a topology, or 'unknown' if unrecognized. It expects an input file of homfly polynomials, one on each line, and can output the identifications in a variety of ways as described below.

#### Usage
Running `idknot` with no other arguments will print out detailed usage information. Some common usage scenarios:

`idknot input_file output_file -nz -u -P` to specify input output files.

`idknot input_file output_file -k` to print conformation identifications line-by-line. 

`idknot -- output_file -nz -u -P` or `idknot -- -nz -u -P > output_file` substitutes the input file for a data stream. 

## Other Software
There are other useful applications that are not included as part of RECOMBO: 

Xinger attempts to remove reidermeister crossings at the extended gauss code level, improving the performance of Homfly. It is avaliable in limited circumstances upon request by emailing a request to: mariel@math.ucdavis.edu. 

Seqconvert is a generally useful application that efficiently converts files filled with self-avoiding polygons between different formats. For example, it can be used to convert CUBE binary formatted polygons to ascii formatted files. It can also perform and output some useful computations for all polygons in a file, such as polygon length and writhe. It is also avaliable in limited circumstances upon request by emailing a request to: mariel@math.ucdavis.edu.

Knotplot is a wonderful knot computation and visualization software developed and maintained by Rob Scharein. In the context of reconnection studies, it is useful for computing the extended gauss codes for each Polygon. Knotplot can be purchased for a small fee at http://knotplot.com/.

---

# Guide
Updated Sep 5, 2020

## Introduction
In this guide, we will walk through setup and all the steps the whole RECOMBO workflow. These steps are: generating randomized self-avoiding polygons, modifying the geometry of existing conformations, identifying knot and link-type, and analysis of the resulting transition data.

For our example, we will start where every project starts, a question: "What linktypes are possible outcomes of topological reconnection on the 5_1 knot, and which transition is most likely to ocurr at random?" If some of those words or concepts are unfamiliar to you, I would recommend pausing to read [this] and [this](https://www.nature.com/articles/s41598-017-12172-2).

## Installation (OSX)
1. Set up the directory structure and pull the repository using Git.
Open the terminal. Start by navigating to the directory where you want to download the repository using the `cd` command. You can create a new directory using the `mkdir` command. 

For example, if we wanted to create a new folder to host the repository called `recombo_public` in our home directory, we might type:
```bash
cd ~
mkdir recombo_public
cd recombo_public
git init
```
2. Download the repository using the `git pull` command, as described in the Git documentation. Continuing the above example, we might use:
```bash
git pull https://github.com/USER/recombo_public master
```
3. Compile and build the software using `make all`. You will see messages from the make process in the terminal. If compilation is not successful, examine the error messages in the terminal window and ensure your compiler and environment are treating the code appropriately.
4. Assuming compilation is successful, navigate to the /src/bin subdirectory to find the executables. Continuing our example, the commands would be:
```bash
cd src/bin
ls .
```
...which should show four executables: `mmc`, `xinger`, `homfly`, and `idknot`.

## Generating randomized self-avoiding polygons
The Multiple Markov Chain BFACF executable (`mmc`) is the most direct path to generating randomized self avoiding polygons of known knot or link-type from scratch. `mmc` is often the starting point for the RECOMBO workflow. Its options and their descriptions are listed above. `mmc` requires a coordinate file specifying an initial SAP conformation to randomize. Initial conformations for many knot and link-types can be found in the initial folder. NAMING SCHEME. There are many parameters that must be chosen, and the interested user should dig into the [RECOMBO appendix] to learn more about the theoretical significance of each one. Most importantly, we must choose the fugacity parameter (`-zmax`: controls average length), the number of conformations to be sampled (`-n`), the number of BFACF moves between samples (`-c`: controls statistical dependence of samples), the number of BFACF moves to take as a warmup to randomize the initial conformation (`-w`: controls the statistical dependence between the first sample and the initial conformation), and the sampling mode (`-mode`: controls what is saved / output from `mmc`). The parameter that controls the number of Markov chains (`-m`) can be set to 1, which makes `mmc` run as an instance of serial BFACF. I recommend the beginner run `mmc` in this mode, since running more than one Markov chain requires a slew of other parameters to be specified (and understood). For the advanced user, these parameters are the lower fugacity parameter (`-zmin`), the swap ratio (`-sr`), the number of steps between attempted swaps (`-s`). I have provided reasonable choices of these parameters below as a guide

To continue our example, we would like to generate a collection of random 5_1 knots. We will return to the project directory, create a directory for our results, and execute mmc like so:
```bash
cd ~/recombo_public
mkdir results
src/bin/mmc initial/5_1 results/5_1 -zmin 0.2000 -zmax 0.2100 -q 1 -sr .8 -s 5 -n 5000 -c 20000 -m 1 -w 100000 -mode s -bfs 5000 +s > 5_1_log.txt
```
When `mmc` finishes, 5000 5_1 knot conformations should have been deposited in a file called `5_1n1.b` in our `results/` sub-directory.

### CMC Mode (Optional)
The `mmc` executable also has the ability to generate conformations and sample as a composite Markov chain. One might be motivated to use this functionality if sampling a wide range of lengths, or if randomization of a particular knot or link-type is too slow in serial BFACF mode. The software orders the BFACF Markov chains by their fugacity parameter, which provides a concept of 'adjacency' between two chains. To run `mmc` in composite Markov chain mode, the `-m` parameter, which specifies the initial number of Markov chains, should be set to an integer value larger than 1. The `-zmin` parameter must also be specified as a value lower than `-zmax`. Any Markov chains beyond 2 will have their fugacity parameter specified automatically by equidistantly chosen points on a linear curve between `-zmin` and `-zmax`. During the warmup, `mmc` will check if the swap-ratio `-sr` (the proportion of attempted conformation exchanges between adjacent chains) is at or above the specified frequency. Failure to swap conformations frequently enough is a consequence of adjacent Markov chains with fugacity parameters too far apart, so new Markov chains will be added as needed with fugacity parameters between them. The user must also choose a number of BFACF moves between attempted conformation swaps with the `-s` option. Small integer values work best, as the goal of MMC is to swap information as frequently as possible within performance limitations.

This is <b>not</b> part of our example, but we might use a commmand like this to generate our same dataset above in composite markov chain mode.
```bash
src/bin/mmc initial/5_1 results/5_1 -zmin 0.2000 -zmax 0.2100 -q 1 -sr .8 -s 5 -n 5000 -c 20000 -m 1 -w 100000 -mode s -bfs 5000 +s > 5_1_log.txt
```

## Performing Topological Reconnection
It may be desirable to modify the geometry of the generated conformations in a controlled way. The RECOMBO suite in its current form was assembled to carry out one such modification, called topological reconnection (see [this] or [this](https://www.nature.com/articles/s41598-017-12172-2) for more details). The `mmc` executable supports topological reconnection at the time of polygon generation through the `-mode f` (filter mode) and `-mode r` (reconnection mode). Filter mode saves only conformations that have sites meeting the arclength specified through both `-minarc` and `-maxarc`, and/or `-targetlength`. Reconnection mode saves conformations on the specified interval like `-mode s`, but also saves a post-reconnection version of any conformations that has a site meeting the specified criteria.

To continue our example, suppose we wanted to generate a larger version of the datasetwe made above, but also perform reconnection in inverted repeat. For our arclength criteria, we will require that at least 6 edges fall on one side of the reconnection site, and we will require the conformation to be exactly 120 edges long.
```bash
`mmc initial/5_1 5_1 -zmin 0.2000 -zmax 0.2100 -q 1 -sr .8 -s 5 -n 50000 -c 20000 -m 1 -w 100000 -mode r -minarc 6 -maxarc 114 -targetlength 120 -seed 42 +s > 3_1_log.txt
```

## Identifying Knot and Link-type
The conformation identification pathway relies on an existing installation of Knotplot with BFACF functionality enabled. Knotplot is commercial software that can be purchased [here](https://knotplot.com/download/). Its creator, Rob scharein, must be contacted via email to provide instructions and a special license which enables the needed special functionality. Once knotplot is installed, the we must locate the stand-alone knotplot executable included in the installation (default location is `/Applications/KnotPlot/utilities/knotplot`), and copy/paste it into the `/knotplot` subdirectory of your local repository.

The first step in the identification pathway is to use knotplot to carry out a sequence of critical tasks. First, the conformation is loaded into knotplot's BFACF engine. Then, a number of BFACF moves is performed allowing only moves that maintain or reduce the size of the conformation. This has the effect of reducing the conformation's compplexity, which . The conformation is taken off the lattice and embedded in R3, and an energy minimization algorithm straightens the conformation out to disambiguate the crossings in the projection. Finally, a projection of the conformation is taken onto the plane and the extended gauss code is computed. We will carry out all these tasks in sequence using a knotplot script. It is neccisary to know the number of components in our conformations of interest, as knot and two-component link conformations require different treatment. 

To continue our example, the following knotplot script would be suitable for identifying the post-reconnection conformations we generated in the previous section:
```
alias idknot "bfacf load; bfacf prob 1 1 0; bfacf step 20000; ago 200; centre; align axes; gauss"
alias idlink "save |tmp; sequence +; load combine |tmp"
gauss noblank yes

gauss open 5_1.egc
sequence open 5_1_after.b
until seqend "idlink; idknot; sequence +"
sequence close
gauss close
```
To execute this script, we should make sure we are in the main directory of our local repository, then invoke knotplot like so:
```bash
cd ~/recombo_public
knotplot/knotplot 
```

The above knotplot script and another example for knotted conformations have been included in the `/knotplot` subdirectory of this repository. I recommend contacting Rob Scharein with questions about KnotPlot's built-in scripting functionality.

## Pipelining

## Analysis of statistical dependence in sampled conformations

---

## Contributors (ordered chronologically)
Rob Schaerine, Reuben Brasher, Robert Stolz, Michelle Flanner, Zihao Zhu, Diwen Lu

Repository currently maintained by Robert Stolz. 

## Funding


## License
We need to put this under a GPL or similar license.
