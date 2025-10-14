---
ospool:
    path: software_examples/bioinformatics/tutorial-fastqc/README.md
---

# Bioinformatics Tutorial: Quality Assessment of Data with FastQC

The first step of most biofinformatic analyses is to assess the quality
of the data you have recieved. In this example, we are working with real
DNA sequencing data from a research project studying E. coli. We will
use a common software,
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), to
assess the quality of the data.

Before we begin, let's download the tutorial's Github repository into
our current working directory using the `git clone` command.

``` shell
git clone https://github.com/osg-htc/tutorial-fastqc.git
```

Let us make sure we are in our `tutorial-fastqc` directory by printing
our working directory:

``` shell
cd ~/tutorial-fastqc
```

    [Errno 2] No such file or directory: '/Users/useradmin/tutorial-fastqc'
    /Users/useradmin/Downloads/tutorial-fastqc-main

``` shell
pwd
```

We should see `/home/<username>/tutorial-fastqc`.

## Workload Components

Before thinking about how to run a list of jobs, let's bring the
components of our workload (data and software) onto this computer.

### Data

For the data, we will be using a dataset originally prepared by [Data
Carpentry](https://datacarpentry.org/). This data includes both the
genome of Escherichia coli (E. coli) and paired-end RNA sequencing reads
obtained from a study carried out by Blount et al. published in
[PNAS](http://www.pnas.org/content/105/23/7899). Additional information
about how the data was modified in preparation for this analysis can be
found on the [Data Carpentry's workshop
website](https://datacarpentry.org/wrangling-genomics/aio.html).

We have a script called `download_data.sh` that will download our
bioinformatic data. Let's go ahead and run this script to download our
data.

``` shell
./data/download_data.sh
```

Our sequencing data files, all ending in .fastq, can now be found in a
folder called `data`.

### Software

The first step of this analysis uses a quality control tool called
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
It's fairly simple to install - just download and unzip. However,
because we want to run this analysis in a distributed way, we've gone
through the steps of making a software container with FastQC installed.
This makes the workload more flexible in where it can run.

It's sometimes possible to find an existing container that has what you
need (for example, the [State Public Health Bioinformatics
Community](https://hub.docker.com/u/staphb)), but for this tutorial,
we'll use a container built by OSPool facilitators specifically for this
tutorial.

``` shell
./software/download_software.sh
```

``` shell
ls software/
```

> [!TIP]
> If you wanted to replicate the container build, you could do so by using this definition file:
> ``` shell
>cat software/fastqc.def
>```
> ```shell
>Bootstrap: docker
>From: ubuntu:24.04
>
>%post
>export DEBIAN_FRONTEND=noninteractive 
>apt update && apt install -yy default-jre unzip perl wget
>wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && \
>  mv fastqc_v0.12.1.zip /opt && \
>  cd /opt && \
>  unzip fastqc_v0.12.1.zip && \
>  rm fastqc_v0.12.1.zip
>
>%environment
>export PATH=/opt/FastQC:$PATH
>```
> And then running this command:
> ``` shell
>apptainer build fastqc.sif software/fastqc.def
>```

## Building Our List of Jobs

To run the fastqc program on one sample the command is:

    fastqc <sample_file>

But we have a list of sample files:

``` shell
ls data/*.fastq
```

So we want to run the `fastqc` command for each of these samples. Our
list of jobs will be based on the list of samples -- we want to submit
one job per sample. To do this, we need to make two things:

-   a list of our samples
-   a "template" for the jobs we want to run.

### Create A List

Creating a list of samples is fairly straightforward -- with some shell
commands piped together, we can create a file that has all of the sample
names:

``` shell
ls data/*fastq | cut -d '/' -f 2 | cut -f1 -d "." > list_of_samples.txt
```

Briefly, this command lists any files ending with "fastq" in our ./data/
directory, cuts anything before (and including) the first "/", cuts
everything after the "." character and saves it to
`list_of_samples.txt`.

``` shell
cat list_of_samples.txt
```

### Job Template

The job template will be the HTCondor submit file. To start out, we're
going to write a submit file that submits a list of one (samples), for
testing.

``` shell
cat fastqc.submit
```
```shell
# HTCondor Submit File: fastqc.submit

# Provide our executable and arguments
executable = /opt/FastQC/fastqc
arguments = SRR2584863_1.trim.sub.fastq

# Provide the container for our software
universe    = container
container_image = software/fastqc.sif

# List files that need to be transferred to the job
transfer_input_files = data/SRR2584863_1.trim.sub.fastq
should_transfer_files = YES
transfer_executable = false

# Tell HTCondor to transfer output to our /results directory
transfer_output_files = SRR2584863_1.trim.sub_fastqc.html
transfer_output_remaps = "SRR2584863_1.trim.sub_fastqc.html = results/SRR2584863_1.trim.sub_fastqc.html"

# Track job information
log = logs/fastqc.log
output = logs/fastqc.out
error = logs/fastqc.err

# Resource Requests
request_cpus = 1
request_memory = 1GB
request_disk = 2GB

# Tell HTCondor to run our job once:
queue 1
```

Some highlights from the submit file:

-   The command `fastqc <sample_file>` has been mapped into the
    `executable` and `arguments` options of the submit file.
-   The data (reads) we need is listed in `transfer_input_files`.
-   The software container image with fastqc installed is listed in
    `container_image`
-   Computational resources needed by the job are indicated with
    `request_*` options.
-   We're being careful to set up a nice organizational structure from
    the start, with our logs and error files in a `logs` folder and
    moving the output files (an html file) to a `results` folder using
    the `transfer_output_remaps` option.

We are now ready to submit our test job!

``` shell
condor_submit fastqc.submit
```

We can check on the status of our job in HTCondor's queue using:

``` shell
condor_q
```

We told HTCondor to store our FastQC output files in the results
directory. Let's take a look at our scientific results:

``` shell
ls -lh results/
```

It's always good practice to look at our standard error, standard out,
and HTCondor log files to catch unexpected output:

``` shell
ls -lh logs/
```

### Scaling Up to a List of Jobs

We can now combine our template and our list of samples to generate a
list of jobs! See what our new submit file looks like:

``` shell
cat many-fastqc.submit
```

Two changes have turned our previous submit file into something that can
submit many jobs at once:

-   We have incorporated our list, `sample_list.txt` in the `queue`
    option at the end of the file. There are different ways to "queue"
    items from a list; we've chosen `queue <variable> from <list.txt>`
    as a good all-purpose option.
-   Wherever our job template has a value that will be unique for every
    job (the sample id), we have replaced the value from our first
    submit file with a variable, `$(sample)`, which was defined in the
    queue statement.

One way to think about this file is as an inverted for-loop for
submitting jobs - where the for statement
`queue sample from sample_list.txt` is at the bottom of the file and the
rest of the file above the for statement is the body of the loop.

We're now ready to submit our list of jobs!

``` shell
condor_submit many-fastqc.submit
```
``` shell
# HTCondor Submit File: many-fastqc.submit

# Provide our executable and arguments
executable = /opt/FastQC/fastqc
arguments = $(sample).trim.sub.fastq

# Provide the container for our software
universe    = container
container_image = software/fastqc.sif

# List files that need to be transferred to the job
transfer_input_files = data/$(sample).trim.sub.fastq
should_transfer_files = YES
transfer_executable = false

# Tell HTCondor to transfer output to our /results directory
transfer_output_files = $(sample).trim.sub_fastqc.html
transfer_output_remaps = "$(sample).trim.sub_fastqc.html = results/$(sample).trim.sub_fastqc.html"

# Track job information
log = logs/$(sample).fastqc.log
output = logs/$(sample).fastqc.out
error = logs/$(sample).fastqc.err

# Resource Requests
request_cpus = 1
request_memory = 1GB
request_disk = 2GB

# Tell HTCondor to run our job once:
queue sample from list_of_samples.txt
```

Notice that using a **single submit file**, we now have **multiple jobs
in the queue**.

We can check on the status of our multiple jobs in HTCondor's queue by
using:

``` shell
condor_q
```

When ready, we can check our results in our `results/` directory:

``` shell
ls -lh results/
```

## Return the output to your local computer

Once you are done with your computational analysis, you will want to move the results to your local computer or to a long term storage location.

> [!IMPORTANT]  
> Remember to always backup your data off of the OSPool (and ideally in multiple places).

Let's practice copying our `.html` files to our local laptop. 

First, open a new terminal. Do not log into your OSPool account. Instead, navigate to where you want the files to go on your computer. We will store them in our `Downloads` folder. 

```shell
cd ~/Downloads
```
Then use `scp` ("secure copy") command to copy our results folder and it's contents:

```shell
scp -r username@hostname:/home/username/tutorial-fastqc/results ./
```
For many files, it will be easiest to create a compressed tarball (.tar.gz file) of your files and transfer that instead of each file individually.

An example of this could be `scp -r username@ap40.uw.osg-htc.org:/home/username/results ./`

Now, open the `.html` files using your internet browser on your local computer. 

### **Congratulations on finishing the first step of a sequencing analysis pipeline!**

