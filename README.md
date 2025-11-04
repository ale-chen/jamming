**Creating jobs for slurm:**
1. module load dSQ
2. dsq --job-file ___.txt --mem 50g -t 20\:00\:00 --partition pi\_ohern --cpus-per-task 1


For parfor:

dsq --job-file dive_mapping_jobs.txt -t 20\:00\:00 --account=pi_bp599 --cpus-per-task 16 --mem-per-cpu 2g

**In dsQ .txt batch file:**

"MATLAB COMMAND; exit; (linebreak for new worker)"
