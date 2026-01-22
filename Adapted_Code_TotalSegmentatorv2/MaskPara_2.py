#####SERVER PARA#####

#Import necessary libraries
import nibabel as nib
import matplotlib.pyplot as plt
import numpy as np
import cv2
import glob
import os
import csv

import SimpleITK as sitk
import six
import radiomics
from radiomics import featureextractor

dataDir = '/home/sdemehr1/data_sdemehr1/TotalSegmentator_Exam6/Radiomics'
radiomicsCSV = os.path.join(dataDir, 'RadiomicsPara4-Exam6.csv')
params = os.path.join(dataDir, "Params.yaml")
headers = None

for scanFilePath in sorted(glob.glob("/home/sdemehr1/data_sdemehr1/MESA_Exam6_NIFTI/MESA4*.nii")):
    for scanFilePath2 in sorted(glob.glob("/home/sdemehr1/data_sdemehr1/TotalSegmentator_Exam6/Segmentations2Comb/MESA4*.nii.gz")):
        if scanFilePath[46:57] == scanFilePath2[71:82]:
            for scanFilePath3 in sorted(glob.glob("/home/sdemehr1/data_sdemehr1/TotalSegmentator_Exam6/Mask_Vert/MESA4*.nii.gz"), key=lambda name: (name[62:73], int(os.path.basename(name)[13:-12]))):
                if scanFilePath2[71:82] == scanFilePath3[62:73]:
                    path, filenames = os.path.split(scanFilePath3)
                    path2 = os.path.basename(path)

                    #Load the scan and extract data using nibabel 
                    scan = nib.load(scanFilePath)
                    scanArray = scan.get_fdata()
                    scan2 = nib.load(scanFilePath2)
                    scanArray2 = scan2.get_fdata()
                    scanArray2[scanArray2 >= 1] = 1

                    try:

                        if filenames[12:-12] == 'T1':
                            scan3 = nib.load(scanFilePath3)
                            scanArray3 = scan3.get_fdata()
                            
                            arrT1 = (scanArray3.mean(axis=(0,1)))
                            start = np.argmax(np.array(arrT1)>0)
                            number = np.where(arrT1>0)[0]
                            numberlist = number.tolist()
                            stop = np.count_nonzero(number) + start
                            T1 = np.zeros(scan3.shape)
                            for i in numberlist:
                                if i <= stop:
                                    T1[:, :, start] = scanArray2[:, :, start]
                                    start+=1
                            
                            T1_img = nib.Nifti1Image(T1, scan2.affine)
                            nib.save(T1_img, "/home/sdemehr1/data_sdemehr1/TotalSegmentator_Exam6/Mask_Para/" + filenames[:11] + "_T1_Para_Mask.nii.gz")
                            T1FilePath = "/home/sdemehr1/data_sdemehr1/TotalSegmentator_Exam6/Mask_Para/" + filenames[:11] + "_T1_Para_Mask.nii.gz"
                            extractor = featureextractor.RadiomicsFeatureExtractor(params)
                            resultT1 = extractor.execute(scanFilePath, T1FilePath, label=1)
                            with open(radiomicsCSV, 'a') as outputFile:
                                writer = csv.writer(outputFile, lineterminator='\n')
                                if headers is None:
                                    headers = list(resultT1.keys())
                                    headers.insert(0, 'Slice')
                                    headers.insert(0, 'Name')
                                    writer.writerow(headers)
                                types1 = [type(k) for k in resultT1.values()]
                                string_list = [str(element) for element in types1]
                                row = []
                                m = 0
                                for i in string_list:
                                    if i == "<class 'numpy.ndarray'>":
                                        row.append(float(list(resultT1.values())[m]))
                                        m+=1
                                    else:
                                        row.append(list(resultT1.values())[m])
                                        m+=1
                                row.insert(0, filenames[12:-12])
                                row.insert(0, filenames[:11])
                                writer.writerow(row)

                        elif filenames[12:-12] == 'T2':
                            scan3 = nib.load(scanFilePath3)
                            scanArray3 = scan3.get_fdata()
                           
                            arrT2 = (scanArray3.mean(axis=(0,1)))
                            start = np.argmax(np.array(arrT2)>0)
                            number = np.where(arrT2>0)[0]
                            numberlist = number.tolist()
                            stop = np.count_nonzero(number) + start
                            T2 = np.zeros(scan3.shape)
                            for i in numberlist:
                                if i <= stop:
                                    T2[:, :, start] = scanArray2[:, :, start]
                                    start+=1
                            
                            T2_img = nib.Nifti1Image(T2, scan2.affine)
                            nib.save(T2_img, "/home/sdemehr1/data_sdemehr1/TotalSegmentator_Exam6/Mask_Para/" + filenames[:11] + "_T2_Para_Mask.nii.gz")
                            T2FilePath = "/home/sdemehr1/data_sdemehr1/TotalSegmentator_Exam6/Mask_Para/" + filenames[:11] + "_T2_Para_Mask.nii.gz"
                            extractor = featureextractor.RadiomicsFeatureExtractor(params)
                            resultT2 = extractor.execute(scanFilePath, T2FilePath, label=1)
                            with open(radiomicsCSV, 'a') as outputFile:
                                writer = csv.writer(outputFile, lineterminator='\n')
                                if headers is None:
                                    headers = list(resultT2.keys())
                                    headers.insert(0, 'Slice')
                                    headers.insert(0, 'Name')
                                    writer.writerow(headers)
                                types1 = [type(k) for k in resultT2.values()]
                                string_list = [str(element) for element in types1]
                                row = []
                                m = 0
                                for i in string_list:
                                    if i == "<class 'numpy.ndarray'>":
                                        row.append(float(list(resultT2.values())[m]))
                                        m+=1
                                    else:
                                        row.append(list(resultT2.values())[m])
                                        m+=1
                                row.insert(0, filenames[12:-12])
                                row.insert(0, filenames[:11])
                                writer.writerow(row)        

                        elif filenames[12:-12] == 'T3':
                            scan3 = nib.load(scanFilePath3)
                            scanArray3 = scan3.get_fdata()
                           
                            arrT3 = (scanArray3.mean(axis=(0,1)))
                            start = np.argmax(np.array(arrT3)>0)
                            number = np.where(arrT3>0)[0]
                            numberlist = number.tolist()
                            stop = np.count_nonzero(number) + start
                            T3 = np.zeros(scan3.shape)
                            for i in numberlist:
                                if i <= stop:
                                    T3[:, :, start] = scanArray2[:, :, start]
                                    start+=1
                            
                            T3_img = nib.Nifti1Image(T3, scan2.affine)
                            nib.save(T3_img, "/home/sdemehr1/data_sdemehr1/TotalSegmentator_Exam6/Mask_Para/" + filenames[:11] + "_T3_Para_Mask.nii.gz")
                            T3FilePath = "/home/sdemehr1/data_sdemehr1/TotalSegmentator_Exam6/Mask_Para/" + filenames[:11] + "_T3_Para_Mask.nii.gz"
                            extractor = featureextractor.RadiomicsFeatureExtractor(params)
                            resultT3 = extractor.execute(scanFilePath, T3FilePath, label=1)
                            with open(radiomicsCSV, 'a') as outputFile:
                                writer = csv.writer(outputFile, lineterminator='\n')
                                if headers is None:
                                    headers = list(resultT3.keys())
                                    headers.insert(0, 'Slice')
                                    headers.insert(0, 'Name')
                                    writer.writerow(headers)
                                types1 = [type(k) for k in resultT3.values()]
                                string_list = [str(element) for element in types1]
                                row = []
                                m = 0
                                for i in string_list:
                                    if i == "<class 'numpy.ndarray'>":
                                        row.append(float(list(resultT3.values())[m]))
                                        m+=1
                                    else:
                                        row.append(list(resultT3.values())[m])
                                        m+=1
                                row.insert(0, filenames[12:-12])
                                row.insert(0, filenames[:11])
                                writer.writerow(row)

                        elif filenames[12:-12] == 'T4':
                            scan3 = nib.load(scanFilePath3)
                            scanArray3 = scan3.get_fdata()
                           
                            arrT4 = (scanArray3.mean(axis=(0,1)))
                            start = np.argmax(np.array(arrT4)>0)
                            number = np.where(arrT4>0)[0]
                            numberlist = number.tolist()
                            stop = np.count_nonzero(number) + start
                            T4 = np.zeros(scan3.shape)
                            for i in numberlist:
                                if i <= stop:
                                    T4[:, :, start] = scanArray2[:, :, start]
                                    start+=1
                            
                            T4_img = nib.Nifti1Image(T4, scan2.affine)
                            nib.save(T4_img, "/home/sdemehr1/data_sdemehr1/TotalSegmentator_Exam6/Mask_Para/" + filenames[:11] + "_T4_Para_Mask.nii.gz")
                            T4FilePath = "/home/sdemehr1/data_sdemehr1/TotalSegmentator_Exam6/Mask_Para/" + filenames[:11] + "_T4_Para_Mask.nii.gz"
                            extractor = featureextractor.RadiomicsFeatureExtractor(params)
                            resultT4 = extractor.execute(scanFilePath, T4FilePath, label=1)
                            with open(radiomicsCSV, 'a') as outputFile:
                                writer = csv.writer(outputFile, lineterminator='\n')
                                if headers is None:
                                    headers = list(resultT4.keys())
                                    headers.insert(0, 'Slice')
                                    headers.insert(0, 'Name')
                                    writer.writerow(headers)
                                types1 = [type(k) for k in resultT4.values()]
                                string_list = [str(element) for element in types1]
                                row = []
                                m = 0
                                for i in string_list:
                                    if i == "<class 'numpy.ndarray'>":
                                        row.append(float(list(resultT4.values())[m]))
                                        m+=1
                                    else:
                                        row.append(list(resultT4.values())[m])
                                        m+=1
                                row.insert(0, filenames[12:-12])
                                row.insert(0, filenames[:11])
                                writer.writerow(row)

                        elif filenames[12:-12] == 'T5':
                            scan3 = nib.load(scanFilePath3)
                            scanArray3 = scan3.get_fdata()
                           
                            arrT5 = (scanArray3.mean(axis=(0,1)))
                            start = np.argmax(np.array(arrT5)>0)
                            number = np.where(arrT5>0)[0]
                            numberlist = number.tolist()
                            stop = np.count_nonzero(number) + start
                            T5 = np.zeros(scan3.shape)
                            for i in numberlist:
                                if i <= stop:
                                    T5[:, :, start] = scanArray2[:, :, start]
                                    start+=1
                            
                            T5_img = nib.Nifti1Image(T5, scan2.affine)
                            nib.save(T5_img, "/home/sdemehr1/data_sdemehr1/TotalSegmentator_Exam6/Mask_Para/" + filenames[:11] + "_T5_Para_Mask.nii.gz")
                            T5FilePath = "/home/sdemehr1/data_sdemehr1/TotalSegmentator_Exam6/Mask_Para/" + filenames[:11] + "_T5_Para_Mask.nii.gz"
                            extractor = featureextractor.RadiomicsFeatureExtractor(params)
                            resultT5 = extractor.execute(scanFilePath, T5FilePath, label=1)
                            with open(radiomicsCSV, 'a') as outputFile:
                                writer = csv.writer(outputFile, lineterminator='\n')
                                if headers is None:
                                    headers = list(resultT5.keys())
                                    headers.insert(0, 'Slice')
                                    headers.insert(0, 'Name')
                                    writer.writerow(headers)
                                types1 = [type(k) for k in resultT5.values()]
                                string_list = [str(element) for element in types1]
                                row = []
                                m = 0
                                for i in string_list:
                                    if i == "<class 'numpy.ndarray'>":
                                        row.append(float(list(resultT5.values())[m]))
                                        m+=1
                                    else:
                                        row.append(list(resultT5.values())[m])
                                        m+=1
                                row.insert(0, filenames[12:-12])
                                row.insert(0, filenames[:11])
                                writer.writerow(row)

                        elif filenames[12:-12] == 'T6':
                            scan3 = nib.load(scanFilePath3)
                            scanArray3 = scan3.get_fdata()
                           
                            arrT6 = (scanArray3.mean(axis=(0,1)))
                            start = np.argmax(np.array(arrT6)>0)
                            number = np.where(arrT6>0)[0]
                            numberlist = number.tolist()
                            stop = np.count_nonzero(number) + start
                            T6 = np.zeros(scan3.shape)
                            for i in numberlist:
                                if i <= stop:
                                    T6[:, :, start] = scanArray2[:, :, start]
                                    start+=1
                            
                            T6_img = nib.Nifti1Image(T6, scan2.affine)
                            nib.save(T6_img, "/home/sdemehr1/data_sdemehr1/TotalSegmentator_Exam6/Mask_Para/" + filenames[:11] + "_T6_Para_Mask.nii.gz")
                            T6FilePath = "/home/sdemehr1/data_sdemehr1/TotalSegmentator_Exam6/Mask_Para/" + filenames[:11] + "_T6_Para_Mask.nii.gz"
                            extractor = featureextractor.RadiomicsFeatureExtractor(params)
                            resultT6 = extractor.execute(scanFilePath, T6FilePath, label=1)
                            with open(radiomicsCSV, 'a') as outputFile:
                                writer = csv.writer(outputFile, lineterminator='\n')
                                if headers is None:
                                    headers = list(resultT6.keys())
                                    headers.insert(0, 'Slice')
                                    headers.insert(0, 'Name')
                                    writer.writerow(headers)
                                types1 = [type(k) for k in resultT6.values()]
                                string_list = [str(element) for element in types1]
                                row = []
                                m = 0
                                for i in string_list:
                                    if i == "<class 'numpy.ndarray'>":
                                        row.append(float(list(resultT6.values())[m]))
                                        m+=1
                                    else:
                                        row.append(list(resultT6.values())[m])
                                        m+=1
                                row.insert(0, filenames[12:-12])
                                row.insert(0, filenames[:11])
                                writer.writerow(row)

                        elif filenames[12:-12] == 'T7':
                            scan3 = nib.load(scanFilePath3)
                            scanArray3 = scan3.get_fdata()
                           
                            arrT7 = (scanArray3.mean(axis=(0,1)))
                            start = np.argmax(np.array(arrT7)>0)
                            number = np.where(arrT7>0)[0]
                            numberlist = number.tolist()
                            stop = np.count_nonzero(number) + start
                            T7 = np.zeros(scan3.shape)
                            for i in numberlist:
                                if i <= stop:
                                    T7[:, :, start] = scanArray2[:, :, start]
                                    start+=1
                            
                            T7_img = nib.Nifti1Image(T7, scan2.affine)
                            nib.save(T7_img, "/home/sdemehr1/data_sdemehr1/TotalSegmentator_Exam6/Mask_Para/" + filenames[:11] + "_T7_Para_Mask.nii.gz")
                            T7FilePath = "/home/sdemehr1/data_sdemehr1/TotalSegmentator_Exam6/Mask_Para/" + filenames[:11] + "_T7_Para_Mask.nii.gz"
                            extractor = featureextractor.RadiomicsFeatureExtractor(params)
                            resultT7 = extractor.execute(scanFilePath, T7FilePath, label=1)
                            with open(radiomicsCSV, 'a') as outputFile:
                                writer = csv.writer(outputFile, lineterminator='\n')
                                if headers is None:
                                    headers = list(resultT7.keys())
                                    headers.insert(0, 'Slice')
                                    headers.insert(0, 'Name')
                                    writer.writerow(headers)
                                types1 = [type(k) for k in resultT7.values()]
                                string_list = [str(element) for element in types1]
                                row = []
                                m = 0
                                for i in string_list:
                                    if i == "<class 'numpy.ndarray'>":
                                        row.append(float(list(resultT7.values())[m]))
                                        m+=1
                                    else:
                                        row.append(list(resultT7.values())[m])
                                        m+=1
                                row.insert(0, filenames[12:-12])
                                row.insert(0, filenames[:11])
                                writer.writerow(row)

                        elif filenames[12:-12] == 'T8':
                            scan3 = nib.load(scanFilePath3)
                            scanArray3 = scan3.get_fdata()
                           
                            arrT8 = (scanArray3.mean(axis=(0,1)))
                            start = np.argmax(np.array(arrT8)>0)
                            number = np.where(arrT8>0)[0]
                            numberlist = number.tolist()
                            stop = np.count_nonzero(number) + start
                            T8 = np.zeros(scan3.shape)
                            for i in numberlist:
                                if i <= stop:
                                    T8[:, :, start] = scanArray2[:, :, start]
                                    start+=1
                            
                            T8_img = nib.Nifti1Image(T8, scan2.affine)
                            nib.save(T8_img, "/home/sdemehr1/data_sdemehr1/TotalSegmentator_Exam6/Mask_Para/" + filenames[:11] + "_T8_Para_Mask.nii.gz")
                            T8FilePath = "/home/sdemehr1/data_sdemehr1/TotalSegmentator_Exam6/Mask_Para/" + filenames[:11] + "_T8_Para_Mask.nii.gz"
                            extractor = featureextractor.RadiomicsFeatureExtractor(params)
                            resultT8 = extractor.execute(scanFilePath, T8FilePath, label=1)
                            with open(radiomicsCSV, 'a') as outputFile:
                                writer = csv.writer(outputFile, lineterminator='\n')
                                if headers is None:
                                    headers = list(resultT8.keys())
                                    headers.insert(0, 'Slice')
                                    headers.insert(0, 'Name')
                                    writer.writerow(headers)
                                types1 = [type(k) for k in resultT8.values()]
                                string_list = [str(element) for element in types1]
                                row = []
                                m = 0
                                for i in string_list:
                                    if i == "<class 'numpy.ndarray'>":
                                        row.append(float(list(resultT8.values())[m]))
                                        m+=1
                                    else:
                                        row.append(list(resultT8.values())[m])
                                        m+=1
                                row.insert(0, filenames[12:-12])
                                row.insert(0, filenames[:11])
                                writer.writerow(row)

                        elif filenames[12:-12] == 'T9':
                            scan3 = nib.load(scanFilePath3)
                            scanArray3 = scan3.get_fdata()
                           
                            arrT9 = (scanArray3.mean(axis=(0,1)))
                            start = np.argmax(np.array(arrT9)>0)
                            number = np.where(arrT9>0)[0]
                            numberlist = number.tolist()
                            stop = np.count_nonzero(number) + start
                            T9 = np.zeros(scan3.shape)
                            for i in numberlist:
                                if i <= stop:
                                    T9[:, :, start] = scanArray2[:, :, start]
                                    start+=1
                            
                            T9_img = nib.Nifti1Image(T9, scan2.affine)
                            nib.save(T9_img, "/home/sdemehr1/data_sdemehr1/TotalSegmentator_Exam6/Mask_Para/" + filenames[:11] + "_T9_Para_Mask.nii.gz")
                            T9FilePath = "/home/sdemehr1/data_sdemehr1/TotalSegmentator_Exam6/Mask_Para/" + filenames[:11] + "_T9_Para_Mask.nii.gz"
                            extractor = featureextractor.RadiomicsFeatureExtractor(params)
                            resultT9 = extractor.execute(scanFilePath, T9FilePath, label=1)
                            with open(radiomicsCSV, 'a') as outputFile:
                                writer = csv.writer(outputFile, lineterminator='\n')
                                if headers is None:
                                    headers = list(resultT9.keys())
                                    headers.insert(0, 'Slice')
                                    headers.insert(0, 'Name')
                                    writer.writerow(headers)
                                types1 = [type(k) for k in resultT9.values()]
                                string_list = [str(element) for element in types1]
                                row = []
                                m = 0
                                for i in string_list:
                                    if i == "<class 'numpy.ndarray'>":
                                        row.append(float(list(resultT9.values())[m]))
                                        m+=1
                                    else:
                                        row.append(list(resultT9.values())[m])
                                        m+=1
                                row.insert(0, filenames[12:-12])
                                row.insert(0, filenames[:11])
                                writer.writerow(row)

                        elif filenames[12:-12] == 'T10':
                            scan3 = nib.load(scanFilePath3)
                            scanArray3 = scan3.get_fdata()
                           
                            arrT10 = (scanArray3.mean(axis=(0,1)))
                            start = np.argmax(np.array(arrT10)>0)
                            number = np.where(arrT10>0)[0]
                            numberlist = number.tolist()
                            stop = np.count_nonzero(number) + start
                            T10 = np.zeros(scan3.shape)
                            for i in numberlist:
                                if i <= stop:
                                    T10[:, :, start] = scanArray2[:, :, start]
                                    start+=1
                            
                            T10_img = nib.Nifti1Image(T10, scan2.affine)
                            nib.save(T10_img, "/home/sdemehr1/data_sdemehr1/TotalSegmentator_Exam6/Mask_Para/" + filenames[:11] + "_T10_Para_Mask.nii.gz")
                            T10FilePath = "/home/sdemehr1/data_sdemehr1/TotalSegmentator_Exam6/Mask_Para/" + filenames[:11] + "_T10_Para_Mask.nii.gz"
                            extractor = featureextractor.RadiomicsFeatureExtractor(params)
                            resultT10 = extractor.execute(scanFilePath, T10FilePath, label=1)
                            with open(radiomicsCSV, 'a') as outputFile:
                                writer = csv.writer(outputFile, lineterminator='\n')
                                if headers is None:
                                    headers = list(resultT10.keys())
                                    headers.insert(0, 'Slice')
                                    headers.insert(0, 'Name')
                                    writer.writerow(headers)
                                types1 = [type(k) for k in resultT10.values()]
                                string_list = [str(element) for element in types1]
                                row = []
                                m = 0
                                for i in string_list:
                                    if i == "<class 'numpy.ndarray'>":
                                        row.append(float(list(resultT10.values())[m]))
                                        m+=1
                                    else:
                                        row.append(list(resultT10.values())[m])
                                        m+=1
                                row.insert(0, filenames[12:-12])
                                row.insert(0, filenames[:11])
                                writer.writerow(row)

                        else:
                            continue

                    except TypeError:
                        pass

                else:
                    continue
