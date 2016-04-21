# run_fsSVM.py
# Author: Zack Greenberg
# Department of Neurology
# Last Editd May 1st, 2014

# This script shows an example of slicing MR data (freesurfer reconstructions)
# into python numpy arrays and using them as inputs to a non-linear support-vector
# machine to perform pattern classification. 
#######################################################################################


#!/usr/bin/python

import os
import csv
from collections import Counter
import matplotlib
from matplotlib import pyplot as plt
from numpy import arange
import sklearn
from sklearn import svm
import numpy as np
from scipy.stats import norm


# Set data dir to MAC folder. The folder is 
# /mnt/fs4/Image_Processing_General/PIDN_ADID_ScanTracking/ on the cloud
MAC_dir = '/Users/zackgreenberg/Desktop/MAC/'

# Use CSV to open the master excel sheet in universal mode
f = open(MAC_dir + 'Spreadsheets/FS_volumes_SVM_mat.csv', 'rU')

# Read the open csv into a python list data structure
master_list =[];
mydata = csv.reader(f)
for row in mydata:
    master_list.append(row)
    
PIDN_list = [];
PIDN_count =0;

for item in master_list:
    PIDN_list.append(item[0])
    PIDN_count +=1;

# list of unique PIDNS
unique_PIDN_list = list(set(PIDN_list))

# Use CSV to open the master excel sheet in universal mode
d = open(MAC_dir + 'Spreadsheets/FS_Diagnosis.csv', 'rU')

# Read the t1 short list to check for any other missed PPG scans 
master_list2 =[];
mydata2 = csv.reader(d)
for row in mydata2:
    master_list2.append(row)

# containers for each subjects diagnosis        
Diagnosis_types=[]   
Diagnosis_list = [];
Diagnosis_count =0;

# Create Containers to count specific diagnoses  
AD_list =[];
AD_count =0;

# All types of Aphasia
Aphasia_list=[];
Aphasia_count = 0;

# PROGRESSIVE SUPRANUCLEAR PALSY
PSP_list=[];
PSP_count =0;

# Mild Cognitiv
MCI_list=[];
MCI_count=0;

# Corticobasal degeneration
CBD_list=[];
CBD_count =0;

# Fronto-temporal Dementia    
FTD_list=[];
FTD_count=0;

# Clinically Normal controls
Control_list=[];
Control_count=0;

# Semantic Dementia
SD_list=[];
SD_count=0;

# Any type of dementia
AllDementias_list=[];
AllDementias_count=0;

# Parkinson's Desease
PD_list=[];
PD_count=0;

# Dementia with lewy bodies
Lewy_list=[]
Lewy_count =0;

# Seizures list
Seizure_list =[];
Seizure_count =0;

# PARANEOPLASTIC SYNDROME
PNPS_list =[];
PNPS_count =0;

# POSTERIOR CORTICAL ATROPHY SYNDROME
PCAS_list =[];
PCAS_count =0;

# BRAIN TUMOR
Tumor_list =[];
Tumor_count =0;

# SPINOCEREBELLAR ATAXIA
SPCA_list =[];
SPCA_count =0;

New_PIDN_list=[];
New_PIDN_count =0;


# get all patient diagnoses
for item in master_list2:
    if item[0] in PIDN_list:
        Diagnosis_list.append([item[0],item[6]]);
        Diagnosis_types.append(item[6]);
        Diagnosis_count += 1;
        
for item in Diagnosis_list:
        if item[1].startswith('AD') or item[1].startswith("ALZHEIMER'S DISEASE"):
            AD_list.append(item[0]);
            AD_count +=1;
            AllDementias_list.append(item[0]);
            AllDementias_count +=1;
        elif item[1].startswith('MCI'):
            MCI_list.append(item[0]);
            MCI_count +=1;
        elif item[1].startswith('PROGRESSIVE SUPRANUCLEAR PALSY') or item[1].startswith('PSP'):
            PSP_list.append(item[0]);
            PSP_count +=1;
        elif item[1].startswith('CBD') or item[1].startswith('CORTICOBASAL DEGENERATION'):
            CBD_list.append(item[0]);
            CBD_count +=1;
        elif item[1].startswith('FTLD') or item[1].startswith('FTD') or item[1].startswith('FRONTOTEMPORAL DEMENTIA'):
            FTD_list.append(item[0]);
            FTD_count +=1;
            AllDementias_list.append(item[0]);
            AllDementias_count +=1;
        elif item[1].startswith('CLINICALLY NORMAL'):
            Control_list.append(item[0]);
            Control_count +=1;
        elif item[1].startswith('PROGRESSIVE NONFLUENT APHASIA SYNDROME') or item[1].startswith('PROGRESSIVE APHASIA - UNSPECIFIED') or item[1].startswith('PROGRESSIVE LOGOPENIC APHASIA SYNDROME'):
            Aphasia_list.append(item[0]);
            Aphasia_count +=1;
        elif item[1].startswith('SEMANTIC DEMENTIA'):
            SD_list.append(item[0]);
            SD_count +=1;
            AllDementias_list.append(item[0]);
            AllDementias_count +=1;
        elif item[1].startswith('POSTERIOR CORTICAL ATROPHY SYNDROME'):
            PCAS_list.append(item[0]);
            PCAS_count +=1;
        elif item[1].startswith('PARANEOPLASTIC SYNDROME'):
            PNPS_list.append(item[0]);
            PNPS_count +=1;
        elif item[1].startswith('SEIZURES'):
            Seizure_list.append(item[0]);
            Seizure_count +=1;
        elif item[1].startswith("PARKINSON'S DISEASE"):
            PD_list.append(item[0]);
            PD_count +=1;
        elif item[1].startswith("BRAIN TUMOR"):
            Tumor_list.append(item[0]);
            Tumor_count +=1;
        elif item[1].startswith("SPINOCEREBELLAR ATAXIA"):
            SPCA_list.append(item[0]);
            SPCA_count +=1;
            

# get rid of repeats
Diagnosis_types = list(set(Diagnosis_types));

# Add the # of unique pidns from each list
Total_Unique_PIDNS = len(list(set(Control_list))) + len(list(set(CBD_list))) + len(list(set(PSP_list))) + len(list(set(MCI_list))) + len(list(set(SD_list))) + len(list(set(PD_list))) + len(list(set(Tumor_list))) + len(list(set(FTD_list))) + len(list(set(Aphasia_list))) + len(list(set(AD_list))) + len(list(set(PCAS_list))) + len(list(set(PNPS_list))) + len(list(set(Seizure_list))) + len(list(set(SPCA_list)))

# Append diagnosis to end of data matrix for each patient
for item in master_list:
    for patient in Diagnosis_list:
        if item[0] == patient[0]:
            item.append(patient[1])


# to write out the labels (diagnosis) on the data sheet
#with open('fsMat_SVM_labeled.csv', 'wb') as file:
#    writer = csv.writer(file)
#    for rows in master_list:
#       writer.writerow(rows)

# Use CSV to open the master excel sheet in universal mode
g = open(MAC_dir + 'Spreadsheets/fsMat_SVM_labeled.csv', 'rU')

# Read in the new labeled csv with diagnosis to add column vector of 1's and 0s
master_list3 =[];
mydata3 = csv.reader(g)
for row in mydata3:
    master_list3.append(row)            
            
# add labeles based 
for item in master_list3:
    if item[-18].startswith('AD') or item[-18].startswith("ALZHEIMER'S DISEASE"):
        item[-17] = '1';
        item[-15] = '1';
        item[-1] = '1';
    else:
        item[-17] = '0';
        item[-15] = '0';
        

    if item[-18].startswith('FTLD') or item[-18].startswith('FTD') or item[-18].startswith('FRONTOTEMPORAL DEMENTIA'):
        item[-16] = '1';
        item[-15] = '1';
        item[-1] = '2';
    else:
        item[-16] = '0';
        item[-15] = '0';
        
        
    if item[-18].startswith('SEMANTIC DEMENTIA'):
        item[-14] = '1';
        item[-15] = '1';
        item[-1] = '3';
    else:
        item[-14] = '0';
        item[-15] = '0';
           
    
    if item[-18].startswith('MCI'):
        item[-13] = '1';
        item[-1] = '4';
    else:
        item[-13] = '0';
        
        
    if item[-18].startswith('CBD') or item[-18].startswith('CORTICOBASAL DEGENERATION'):
        item[-12] = '1';
        item[-1] = '5';
    else:
        item[-12] = '0';
        
        
    if item[-18].startswith('PROGRESSIVE SUPRANUCLEAR PALSY') or item[-18].startswith('PSP'):
        item[-11] = '1';
        item[-1] = '6';
    else:
        item[-11] = '0';
        
    
    if item[-18].startswith('PROGRESSIVE NONFLUENT APHASIA SYNDROME') or item[-18].startswith('PROGRESSIVE APHASIA - UNSPECIFIED') or item[-18].startswith('PROGRESSIVE LOGOPENIC APHASIA SYNDROME'):
        item[-10] = '1';
        item[-1] = '7';
    else:
        item[-10] = '0';
        
        
    if item[-18].startswith('CLINICALLY NORMAL') or item[-18].startswith('Clinically Normal'):
        item[-9] = '1';
        item[-1] ='8';
    else:
        item[-9] = '0';
        
        
    if item[-18].startswith("PARKINSON'S DISEASE"):
        item[-8] = '1';
        item[-1] ='9';
    else:
        item[-8] = '0';
        
        
    if item[-18].startswith("LEWY"):
        item[-7] = '1';
        item[-1] ='10';
    else:
        item[-7] = '0';
        
        
    if item[-18].startswith("SEIZURES"):
        item[-6] = '1';
        item[-1] ='11';
    else:
        item[-6] = '0';
        
        
    if item[-18].startswith('PARANEOPLASTIC SYNDROME'):
        item[-5] = '1';
        item[-1] ='12';
    else:
        item[-5] = '0';
        
    if item[-18].startswith('POSTERIOR CORTICAL ATROPHY SYNDROME'):
        item[-4] = '1';
        item[-1] ='13';
    else:
        item[-4] = '0';
        
    if item[-17].startswith('BRAIN TUMOR'):
        item[-3] = '1';
        item[-1] ='14';
    else:
        item[-3] = '0';
        
    if item[-18].startswith('SPINOCEREBELLAR ATAXIA'):
        item[-2] = '1';
        item[-1] ='15';
    else:
        item[-2] = '0';

# Containers for all   
X = [];
Y_all =[];
Y_AD =[];
Y_MCI =[];
Y_Control =[]
Y_SD =[];
Y_FTD=[];
Y_PSP=[];


# fill y labels
for item in X_data:
     Y_all.append(item[-15])

# Do SVM classification for all dementia types vs all non-dementia
All_clf = svm.SVC()
All_clf.fit(X,Y_all)
Y_predict_all = clf.predict(X)
All_acc=0
for item in Y_all:
     for prediction in Y_predict_all:
        if item == prediction:
           all_acc +=1
     break;


# convert to numpy array of floats
X= X_data.astype(float)

# look at the mean and std of one feature vector (column)
mean = np.mean(X[:0])
std = np.std(X[:0])
X_col = sorted(X[:,0])

# fit a gausian to plot the distribution
fit = norm.pdf(X_col, mean, std)
pl.plot(X_col, fit, '-o')
pl.hold();
pl.hist(X_col,normed=True)
pl.show()
