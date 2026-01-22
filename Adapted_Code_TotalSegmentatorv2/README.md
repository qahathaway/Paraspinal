### DL-Derived 3D Paraspinal Muscle Segmentation

#### Step-1-NIFTI.py
Convert DICOM files into a single NIFTI file

#### Step-2-Paraspinal-Muscle.sh
Use the task "autochthon_left" and "autochthon_right" to segment the left and right paraspinal muscle, respectively.

#### Step-3-Combine-Paraspinal-Muscles.py
We combined the left and right paraspinal muscle volumes.

#### Step-4
We have previously published how we utilized TotalSegmentator to extract verterbal body information here:
https://pubmed.ncbi.nlm.nih.gov/40067103/

Please refer to https://github.com/qahathaway/vBMD/tree/main/Adapted_Code_TotalSegmentatorv2 for further information.

With the spatial information from the prior vertebral body extraction, we proceeded to extract level-by-level paraspinal muscle information.

#### Step-5-T1-T10-Paraspinal-Muscle.py

#### Overview of Methods
The original algorithm, TotalSegmentator v2 (1), has been recently developed using the nnU-net framework (2) on eight different clinical sites and 16 different scanners. Both 1.5- and 3.0-mm isotropic resolution CT scans were trained through batches of iterative annotation. Data augmentation was applied, and 1.5 mm slices were trained for 4,000 epochs, whereas 3.0 mm slices were trained for 8,000 epochs. The algorithm is capable of volumetric segmentation of 117 anatomic structures from a wide variety of CT acquisitions, demonstrating its generalizability (1). The volumetrically segmented structures include vertebral bodies and paraspinal muscles with DICE coefficients of 0.94 on internal testing and 0.93 on external testing for all anatomical regions measured.

1. Wasserthal J, Breit HC, Meyer MT, Pradella M, Hinck D, Sauter AW, Heye T, Boll DT, Cyriac J, Yang S, Bach M, Segeroth M. TotalSegmentator: Robust Segmentation of 104 Anatomic Structures in CT Images. Radiol Artif Intell 2023;5(5):e230024. doi: 10.1148/ryai.230024
