bsub -o output_file -W 72:00 -R "rusage[mem=50000]" python Main_script_fitGPs_2.py /cluster/home/pmonika/Fitting/ ../../../scratch/pmonika/Fitting/ Betaine_2 FT
bsub -o output_file -W 72:00 -R "rusage[mem=50000]" python Main_script_fitGPs_2.py /cluster/home/pmonika/Fitting/ ../../../scratch/pmonika/Fitting/ Betaine_2 HT
bsub -o output_file -W 72:00 -R "rusage[mem=50000]" python Main_script_fitGPs_2.py /cluster/home/pmonika/Fitting/ ../../../scratch/pmonika/Fitting/ Proline_2 FT
bsub -o output_file -W 72:00 -R "rusage[mem=50000]" python Main_script_fitGPs_2.py /cluster/home/pmonika/Fitting/ ../../../scratch/pmonika/Fitting/ Proline_2 HT
