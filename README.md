# Paraspinal Muscle Biomarkers on Chest CT Predicts Vertebral Bone Mineral Density Loss and Fractures: An Opportunistic Deep Learning Study from MESA

## TotalSegmentator 2.0 for 3D Segmentation
#### This code was adapted from https://github.com/wasserth/TotalSegmentator with modifcations as detailed within this repository.
#### Please cite the manuscript by Wasserthal et al. at: https://pubs.rsna.org/doi/10.1148/ryai.230024 or https://pmc.ncbi.nlm.nih.gov/articles/PMC10546353/


### Key Results:
1.   Participants in the lowest quartile of paraspinal muscle attenuation had significantly lower baseline vertebral bone mineral density (vBMD) (176 vs 205 mg/cm³, P≤0.001), were older, White, and more often female.
2.   Adding muscle attenuation to clinical risk factors improved prediction of follow-up vBMD explained an additional ~7% of variance.
3.   Muscle attenuation provided incremental value for incident vertebral fracture prediction beyond vBMD alone (AUC 0.82 vs 0.77, P = 0.02).


### Summary Statement:
Deep learning–derived paraspinal muscle attenuation on noncontrast chest CT independently and incrementally predicts longitudinal vertebral bone mineral density loss and incident vertebral fractures beyond clinical risk factors.

### Background:
Muscle and bone loss commonly coexist with aging and share common risk factors, but paraspinal muscles could play a critical local biomechanical role in spinal loading.

### Purpose:
We evaluated whether CT-derived paraspinal muscle measures independently predict longitudinal thoracic vertebral bone mineral density loss.

### Materials and Methods:
Using the Multi-Ethnic Study of Atherosclerosis (MESA) Lung Study, 1316 participants’ chest CT scans in MESA Exam 5 (baseline) and 6 (median follow-up 6.2 years) were included after excluding those with incomplete clinical Fracture Risk Assessment Tool without BMD (FRAXnb) information or failed segmentations. Multilevel (T1-T10) deep learning (DL)–based 3D segmentation of thoracic paraspinal muscles was performed using TotalSegmentator. 3D paraspinal muscle measures (volume and attenuation) were dichotomized (1st vs. 2-4th quartiles) to identify the lowest-quality muscle group and predicted vertebral bone mineral density (vBMD) and incident vertebral fractures at Exam 6, adjusted for clinical risk factors and compared against baseline vBMD (Exam 5).

### Results:
DL-derived 3D segmentations were correlated with manual segmentations (Dice score: 0.94). At Exam 5, participants with low vs. high (i.e., 1st quartile (HU: 4) vs. 2-4th quartiles (27 HU)) paraspinal muscle attenuation shared common osteoporotic risk factors, such as being older (70 vs. 66 yrs.), White (49% vs. 30%), and female (77% vs. 45%). Low paraspinal muscle attenuation (176 vs. 205 mg/cm3, P≤0.001) and volume (180 vs. 203 mg/cm3, P≤0.001) were also associated with lower vBMD when compared to the higher 2-4th quartiles. Paraspinal muscle attenuation at Exam 5 was the strongest predictor of vBMD at Exam 6, beyond FRAXnb alone (Adjusted R²: 0.27 vs. 0.20, P<0.001), and provided incremental value in predicting incident vertebral fractures in Exam-6, above vBMD alone (AUC: 0.82 vs. 0.77, P=0.02).

### Conclusion:
DL-derived paraspinal muscle attenuation provides independent and incremental value for predicting vertebral bone mineral density and vertebral fractures.


#### ClinicalTrials.gov: NCT00005487


## Overview of Study Design and 3D Algorithm Development
![alt text](https://github.com/qahathaway/Paraspinal/blob/main/Overview.jpg)
#### 3D Segmentation Overview.
(A) Volumetric segmentation of paraspinal musculature was extracted in 1316 participants. (B) Manual segmentation was performed, compared to volumetric segmentation, and Dice Similarity Coefficient (DSC) and the Intersection over Union (IoU) were calculated. (C) The 3D segmentations from T1-T10 were individually extracted and combined into a single volume for analyses. (D) Muscle volume and median muscle attenuation were derived from the segmentation volumes. (E) Analyses were performed to assess the correlation of paraspinal muscle attenuation and volume with bone mineral density and the ability to predict vertebral fractures (VFx).



### Code is made freely available for academic research and teaching. The code within this repository cannot be freely used for commercial applications.
