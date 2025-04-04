**Creating jobs for slurm:**
1. module load dsQ
2. dsq --job-file ___.txt --mem 50g -t 20\:00\:00 --partition pi\_ohern --cpus-per-task 1

**In dsQ .txt batch file:**

"MATLAB COMMAND; exit; (linebreak for new worker)"
