#!/bin/bash

# The name of the job
#$ -N JOB_ID

# The name of the execution queue
#$ -q short.q

# Specifies the interpreting shell for the job
#$ -S /bin/bash

# Specifies that all environment variables active within the qsub utility be exported to the context of the job.
#$ -V

# Declare the job as located in the current working directory.
#$ -cwd

# Defines  or  redefines  the  path  used  for the standard out and error stream of the job. Either absolute path or relative to the -cwd path
# #$ -o ./stdout/out.out
# #$ -e ./stderr/err.err

#$ -pe smp 2

# Defines  the size of the job array (numbers of identical jobs only differentiated by the parameter passed to the job)
#$ -t 1-NUM_FILES
# extract the content of the Nth number line from an argument file containing the parameters for the Nth number job

 TOPARSEFILE=VAR_FILE_PATH
 CURENTPARAM=$(sed -n -e "$SGE_TASK_ID p" $TOPARSEFILE)


# Create scratch 
 export MCR_CACHE_ROOT="/groups/dmf/gr_Pelet/YeastQuantX/scratch/spelet1_"$JOB_NAME"_"$SGE_TASK_ID"/mcrCache"
 mkdir -p /groups/dmf/gr_Pelet/YeastQuantX/scratch/spelet1_"$JOB_NAME"_"$SGE_TASK_ID"/mcrCache

/groups/dmf/gr_Pelet/YeastQuantX/YQ_CX/src/YQ_VX $CURENTPARAM


# Remove scratch space
rm -rf  /groups/dmf/gr_Pelet/YeastQuantX/scratch/spelet1_"$JOB_NAME"_"$SGE_TASK_ID"/mcrCache