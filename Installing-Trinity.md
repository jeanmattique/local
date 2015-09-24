# Installing Trinity

After [downloading](https://github.com/trinityrnaseq/trinityrnaseq/releases) the software to a Linux server, simply type 
   
    make 

in the base installation directory.  This should build Inchworm and Chrysalis, both written in C++.  Butterfly should not require any special compilation, as its written in Java and already provided as portable precompiled software, but *Java-1.7* (or higher) is required.

Afterwards, you may want to build the additional plugin components that provide support for downstream analyses (such as abundance estimation using RSEM), in which case you would then type:

    make plugins


If you encounter any errors in building the RSEM software, simply

    cd trinity-plugins/tmp.rsem
     
    make

and assuming that succeeds, you can then cd back to the main Trinity installation directory and retype 'make plugins' to continue the remaining build.

Additional tools required for running Trinity include:

- [bowtie-1](http://bowtie-bio.sourceforge.net/index.shtml)

Trinity has been tested and is supported on Linux.


To test your installation of Trinity, try assembling the small sample data set provided with Trinity like so:

    cd sample_data/test_Trinity_Assembly/
    
    ./runMe.sh

## [OPTIONAL] Adapting Trinity to a computing grid for parallel processing of naively parallel steps

Trinity has many parallel-components, all of which can benefit from having multiple CPUs on a single server, but there are also cases such as in Chrysalis and Butterfly where tens of thousands to hundreds of thousands of commands can be executed in parallel, each having independent inputs and outputs.  These naively-parallel commands can be most efficiently computed in the context of a compute farm, submitting each of the commands (or batches of them) to individual nodes on the computing grid.  

There are several different computing grid job management systems that are in common use. Trinity supports LSF, SGE, SLURM, and PBS.  To leverage one, simply run 'Trinity --grid_conf your_conf_file.txt', where your_conf_file.txt is a very simple configuration file that indicates parameters for the grid job submission. For example, at the Broad and using LSF, a configuration file might contain the following:


    #-------------------------------------------------------------------------------------------
    # grid type: 
    grid=LSF
       
    # template for a grid submission
    cmd=bsub -q regevlab -R "rusage[mem=10]"
    # note -e error.file -o out.file are set internally, so dont set them in the above cmd. 
     
    # uses the LSF feature to pre-exec and check that the file system is mounted before executing.
    # this helps when you have some misbehaving grid nodes that lost certain file mounts.
    mount_test=T
     
    ##########################################################################################
    # settings below configure the Trinity job submission system, not tied to the grid itself.
    ##########################################################################################
     
    # number of grid submissions to be maintained at steady state by the Trinity submission system 
    max_nodes=500
     
    # number of commands that are batched into a single grid submission job.
    cmds_per_node=100
    
    #--------------------------------------------------------------------------------------------


where the above indicates that LSF is the grid type (either LSF or SGE are supported), the queue to submit to is our 'regevlab' named queue, and memory is set to 10 gigabytes. Up to 500 jobs will be submitted at any given time (throttled by the Trinity-included job management system), and the jobs are batched at 10 commands per submission (so, for example, 10 butterfly jobs will be submitted as a single grid job, each being executed serially).

For SGE, at the Broad Institute, we might specify a configuration:

    #--------------------------------------------------------------------------------------------
    # grid type: 
    grid=SGE
    # template for a grid submission
    cmd=qsub -V -cwd
    # number of grid submissions to be maintained at steady state by the Trinity submission system 
    max_nodes=500
    # number of commands that are batched into a single grid submission job.
    cmds_per_node=1
    #--------------------------------------------------------------------------------------------

where, SGE is indicated as the grid type.  We don't need to specify a queue name, apparently, as it gets submitted to the default queue, and the default memory allocation is sufficient. The project_code can also be left blank unless your SGE configuration requires it.  The maximum number of nodes to throttle the jobs at (500) and the number of commands executed in a single grid job (10) is the same as what we show above for our LSF configuration.

Likewise, for SLURM, we have:

    #---------------------------------------------------------------------------------------------
    # grid type: 
    grid=SLURM
    # template for a grid submission
    cmd=sbatch -p queue_name --mem=10000 --time=02:00:00 
    # number of grid submissions to be maintained at steady state by the Trinity submission system 
    max_nodes=4000
    # number of commands that are batched into a single grid submission job.
    cmds_per_node=20
    #----------------------------------------------------------------------------------------------


Example configuration files are provided under $TRINITY_HOME/htc_conf

