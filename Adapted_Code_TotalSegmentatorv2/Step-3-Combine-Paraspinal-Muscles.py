#Import necessary libraries
import nibabel as nib
import numpy as np
import cv2
import glob
import os

for scanFilePath in sorted(glob.glob("/home/sdemehr1/data_sdemehr1/TotalSegmentator_Exam6/Segmentations2/**/*_right.nii.gz")):
    for scanFilePath2 in sorted(glob.glob("/home/sdemehr1/data_sdemehr1/TotalSegmentator_Exam6/Segmentations2/**/*_left.nii.gz")):
        if scanFilePath[67:78] == scanFilePath2[67:78]:
            path, filenames = os.path.split(scanFilePath)
            path2 = os.path.basename(path)
            scanParaR = nib.load(scanFilePath)
            scanParaR1 = scanParaR.get_fdata()

            scanParaL = nib.load(scanFilePath2)
            scanParaL1 = scanParaL.get_fdata()

            Para = cv2.add(scanParaR1, scanParaL1)
            Para_img = nib.Nifti1Image(Para, scanParaR.affine)
            nib.save(Para_img, "/home/sdemehr1/data_sdemehr1/TotalSegmentator_Exam6/Segmentations2Comb/" + path2 + "_Para.nii.gz")
