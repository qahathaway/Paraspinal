#!/bin/bash -l
#SBATCH -J TotalSeg_Para
#SBATCH --time=72:0:0
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH --mem-per-cpu=4GB
#SBATCH --quiet


module load anaconda
conda activate TotalSeg

for FILE in /home/sdemehr1/data_sdemehr1/MESA_Exam6_NIFTI/*.nii;
do filename=$(basename "$FILE");
filename="${filename%.*}";
TotalSegmentator -i $FILE -o /home/sdemehr1/data_sdemehr1/TotalSegmentator_Exam6/Segmentations2/$filename -ta total --roi_subset autochthon_left autochthon_right;
done
