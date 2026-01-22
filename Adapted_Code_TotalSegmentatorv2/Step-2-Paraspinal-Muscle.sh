#!/bin/bash -l
#SBATCH -J TotalSeg_Para
#SBATCH --time=72:0:0
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH --mem-per-cpu=4GB
#SBATCH --quiet


module load anaconda
conda activate TotalSeg

for FILE in /path/to/NIFTI/repository/*.nii;
do filename=$(basename "$FILE");
filename="${filename%.*}";
TotalSegmentator -i $FILE -o /path/to/segmentations/folder/$filename -ta total --roi_subset autochthon_left autochthon_right;
done
