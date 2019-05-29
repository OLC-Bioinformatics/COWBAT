#!/bin/bash

#Number of nodes to request. Leave at 1 unless you're running an MPI job.
#SBATCH -N 1
# Set number of processors you want for your job here.
#SBATCH --ntasks=55
# Set the amount of memory you want (in megabytes, so this is 191 gigs)
#SBATCH --mem=191000
# Set amount of time to give your job here. (D-HH:MM)
#SBATCH --time=1-00:00
# Standard out will go here (%j will be replaced by job ID when script is run).
#SBATCH -o /mnt/nas2/redmine/bio_requests/12430/slurm_logs/job_%j.out
# Standard error will go here (percent j replaced by job ID when script is run).
#SBATCH -e /mnt/nas2/redmine/bio_requests/12430/slurm_logs/job_%j.err

# Your code goes here. 
docker run -u ubuntu -i -v /mnt/nas2:/mnt/nas2 --name cowbat --rm cowbat:latest /bin/bash -c "source activate cowbat && python3 assembly_pipeline.py -s /mnt/nas2/redmine/bio_requests/12430/fastqs -r /mnt/nas2/databases/assemblydatabases/0.3.4"
