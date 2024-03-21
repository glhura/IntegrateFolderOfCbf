#!/usr/bin/env python

import argparse
import os
import sys
import multiprocessing
import math
from ctypes import *
import numpy as np
from numpy.ctypeslib import ndpointer

import ID_SeqFileType_In_Dir
import ProcessIzeroFile
import ChangeFileName
import GetIzeroValueForFile
import SubtractTwoDats

#Where Iend knows exposure begins and ends
ThresholdForIzero = 0.01
#Number of Processors to use
NumProc = 8

#The .so c code for integrating
so_file = "./pil_2m_integ_cbf.so"


#This function takes a file list and splits the list for each processor.
#Sometimes the split is not even and there are some orphan files
def SplitFileListN(FileList):
	NumFiles = len(FileList)
	SingleProcessorFlag = 0
	NumFilesPerProcess = int(math.trunc(NumFiles/(NumProc-1)))
	LeftOver = NumFiles - (NumFilesPerProcess*(NumProc-1))
	#print(NumFilesPerProcess, LeftOver)
	if NumFilesPerProcess < 1:
		Process = []
		Process.append([])
		Process[0] = FileList
	if NumFilesPerProcess >= 1:
		Process = []
		for i in range(0,NumProc):
			Process.append([])
		for i in range(0,(NumProc-1)):
			StartList = int(i*NumFilesPerProcess)
			FinishList = int((i+1)*NumFilesPerProcess)
			Process[i] = FileList[StartList:FinishList]
			
		Process[(NumProc-1)]=FileList[FinishList:(FinishList+LeftOver)]
	#print(Process)
	return(Process)

#Process the Izero file once here and get its normalization
def OperationsOnOneFile(File,SourceFlag):
	NormValues = GetIzeroValueForFile.GetNormValue(File,ThresholdForIzero)
	#print(SourceFlag)
	if int(SourceFlag) == 1:
		UserNorm = NormValues[0]
	if int(SourceFlag) == 2:
		UserNorm = NormValues[1]
	return(UserNorm)


"""This Function Identifies the cbf files in an input folder, it processes the Izero file
It then splits the folder into groups of 8 and returns 8 lists of files"""
def OperationsOnOneFold(Fold):
    if not os.path.exists(Fold):
        sys.exit("%s Does Not Exist" %  Fold)
    FoldFiles = ID_SeqFileType_In_Dir.ID_files(Fold, ".cbf")
    NumFrames = len(FoldFiles)
    directory_content_list = os.listdir(Fold)
    IzeroTxt = 0
    for list_item in directory_content_list:
        if (str(list_item).endswith("_Izero.txt")):
            IzeroTxt = 1
            FullListItem = Fold + "/" + list_item
            RunCommand = ProcessIzeroFile.MakeNormList(FullListItem,NumFrames,ThresholdForIzero)   
            
    if IzeroTxt == 0:
        sys.exit("No Izero.txt File in %s" % SampFold)
	
    SplitListN = SplitFileListN(FoldFiles)

    #print(NumFrames,SplitListN)
    return(NumFrames,SplitListN)

"""This function checks to see if Samp and Buff Folder exist, creates the out folders if they dont exist
Calls a function that creates NunProc (an int specified at top of code) lists of files from both samp and buff folders
it takes the number of files in each of the NumProc lists, it then reads the mask file and Exptparams and stores it into memory"""

def try_to_multi(SampFold,BuffFold,BeforeOrAfterFlag,SourceFlag,FullFlag):
    
    if not os.path.exists(SampFold):
        sys.exit("%s Does Not Exist" %  SampFold)
    if not os.path.exists(BuffFold):
        sys.exit("%s Does Not Exist" %  BuffFold)
			
    MakeSubDir = ("Results/Subtracted/" + SampFold + "_results")
    if BeforeOrAfterFlag == "b":
        MakeSubDir = (MakeSubDir + "/Buffer1")
    if BeforeOrAfterFlag == "a":
        MakeSubDir = (MakeSubDir + "/Buffer2")
    MakeUnSubSampDir = ("Results/Unsubtracted/" + SampFold + "_results")
    MakeUnSubBuffDir = ("Results/Unsubtracted/" + BuffFold + "_results")

    if not os.path.exists(MakeSubDir):
        os.makedirs(MakeSubDir)
    if not os.path.exists(MakeUnSubSampDir):
        os.makedirs(MakeUnSubSampDir)
    if not os.path.exists(MakeUnSubBuffDir):
        os.makedirs(MakeUnSubBuffDir)	
	
    StatsSampFiles = OperationsOnOneFold(SampFold)
    StatsBuffFiles = OperationsOnOneFold(BuffFold)
    NumFileSamp = StatsSampFiles[0]
    NumFileBuff = StatsBuffFiles[0]
    SampFileList = StatsSampFiles[1]
    BuffFileList = StatsBuffFiles[1]
    HowBigSamp = len(SampFileList)
    HowBigBuff = len(BuffFileList)
    
    
 
    if HowBigSamp != HowBigBuff:
        sys.exit("Sample Folder And Buffer Folder Dont Have Same Number of Files #1")
    group_samp = np.zeros(HowBigSamp)
    group_buff = np.zeros(HowBigSamp)
    for i in range(0,HowBigSamp):
        group_samp[i] = len(SampFileList[i])
        group_buff[i] = len(BuffFileList[i])
        if group_samp[i] != group_buff[i]: 
            sys.exit("Sample Folder And Buffer Folder Dont Have Same Number of Files #2")
        print(group_samp[i])
    
    mask_file = open("X_Y_mask.asc", 'r')
    mask_data = np.loadtxt(mask_file, dtype=np.int32)
    mask_file.close()
    mask_lines = int (mask_data.shape[0])
    
    expt_params = open("ExptParams", 'r')
    line = expt_params.readline()
    line = line.rstrip()
    split_line = line.split(' ')
    exp_arr = np.array(split_line, dtype = c_double)    

    for i in range(0,3):
        line = expt_params.readline()
        line = line.rstrip()
        split_line = line.split(' ')
        med_arr = np.array(split_line, dtype = c_double)
        exp_arr = np.append(exp_arr,med_arr)
    
    expt_params.close()    
    num_params = int (exp_arr.shape[0])    
    
    _params = ndpointer(c_double, flags='C')
    _mask = ndpointer(dtype=np.uintp, ndim=1, flags='C')
    C_prog = CDLL(so_file)
    C_prog.intensity_matrix_gen.argtypes = [c_char_p, c_size_t, _params, c_size_t, _mask]
    C_prog.intensity_matrix_gen.restype = None
    
    assert(mask_data.ndim == 2)
    mask_pp = (mask_data.__array_interface__['data'][0] +  np.arange(mask_data.shape[0])*mask_data.strides[0]).astype(np.uintp)
       
    """file_name1 = SampFileList[0][0].encode("UTF-8")
    C_prog.intensity_matrix_gen(file_name1, num_params, exp_arr, mask_lines, mask_pp)"""
    
    CurrentDir = os.getcwd()
    if HowBigSamp == 1:
        print(group_samp[0])

        for i in range(0, int (group_samp[0])):
            samp_name = SampFileList[0][i].encode("UTF-8")
            buff_name = BuffFileList[0][i].encode("UTF-8")
            print(samp_name)
            print(buff_name)
            C_prog.intensity_matrix_gen(samp_name, num_params, exp_arr, mask_lines, mask_pp)
            C_prog.intensity_matrix_gen(buff_name, num_params, exp_arr, mask_lines, mask_pp)
        for i in range(0, int (group_samp[0])):
            WrittenDatSamp = ChangeFileName.ChangeName(SampFileList[0][i],"dat",0)
            SplitTarget = WrittenDatSamp.split("/")
            WorkingFile = SplitTarget[-1]
            CreateUnsubFileName = os.path.join(MakeUnSubSampDir,WorkingFile)
            WrittenDatSamp = os.path.join(CurrentDir,WrittenDatSamp)
            CreateUnsubFileName = os.path.join(CurrentDir,CreateUnsubFileName)
            print(CreateUnsubFileName)
            os.rename(WrittenDatSamp, CreateUnsubFileName)
            WrittenDatSamp = ChangeFileName.ChangeName(BuffFileList[0][i],"dat",0)
            SplitTarget = WrittenDatSamp.split("/")
            WorkingFile = SplitTarget[-1]
            CreateUnsubFileName = os.path.join(MakeUnSubBuffDir,WorkingFile)
            WrittenDatSamp = os.path.join(CurrentDir,WrittenDatSamp)
            CreateUnsubFileName = os.path.join(CurrentDir,CreateUnsubFileName)
            print(CreateUnsubFileName)
            os.rename(WrittenDatSamp, CreateUnsubFileName)
            
    if HowBigSamp == NumProc:
        plist = []
        for i in range(0,HowBigSamp):
              proc_val = "p" + str (i)
              plist.append(proc_val)
        print(plist)
        for i in range(0,HowBigSamp):
            for j in range(0, int (group_samp[i])):
                samp_name = SampFileList[i][j].encode("UTF-8")
                buff_name = BuffFileList[i][j].encode("UTF-8")
                print(samp_name)
                print(buff_name)
            for k in range(0, int (group_samp[i])):
                plist[k] = multiprocessing.Process(target=C_prog.intensity_matrix_gen, args=(SampFileList[i][k].encode("UTF-8"), num_params, exp_arr, mask_lines, mask_pp))
                plist[k].start()
            for l in range(0, int (group_samp[i])):
                plist[l].join()
            for k in range(0, int (group_samp[i])):
                plist[k] = multiprocessing.Process(target=C_prog.intensity_matrix_gen, args=(BuffFileList[i][k].encode("UTF-8"), num_params, exp_arr, mask_lines, mask_pp))
                plist[k].start()
            for l in range(0, int (group_samp[i])):
                plist[l].join()
            for m in range(0, int (group_samp[i])):
                WrittenDatSamp = ChangeFileName.ChangeName(SampFileList[i][m],"dat",0)
                SplitTarget = WrittenDatSamp.split("/")
                WorkingFile = SplitTarget[-1]
                CreateUnsubFileName = os.path.join(MakeUnSubSampDir,WorkingFile)
                WrittenDatSamp = os.path.join(CurrentDir,WrittenDatSamp)
                CreateUnsubSampFileName = os.path.join(CurrentDir,CreateUnsubFileName)
                #print("This Spot", CreateUnsubFileName)
                os.rename(WrittenDatSamp, CreateUnsubFileName)
                WrittenDatBuf = ChangeFileName.ChangeName(BuffFileList[i][m],"dat",0)
                SplitTarget = WrittenDatBuf.split("/")
                WorkingBufFile = SplitTarget[-1]
                CreateUnsubBufFileName = os.path.join(MakeUnSubBuffDir,WorkingBufFile)
                WrittenDatBuf = os.path.join(CurrentDir,WrittenDatBuf)
                CreateUnsubBufFileName = os.path.join(CurrentDir,CreateUnsubBufFileName)
                print(CreateUnsubBufFileName)
                os.rename(WrittenDatBuf, CreateUnsubBufFileName)
                SampScale = OperationsOnOneFile(SampFileList[i][m],SourceFlag)
                BuffScale =  OperationsOnOneFile(BuffFileList[i][m],SourceFlag)
                #print(CreateUnsubFileName,SampScale,BuffScale)
                TheSubtractedFile = SubtractTwoDats.SubtractDats(CreateUnsubSampFileName,SampScale,CreateUnsubBufFileName,BuffScale,MakeSubDir)
            #print("Got Through a Round")
            
    
    return MakeSubDir
    


if __name__=="__main__":
    """Takes two folders as input and a few flags - (which normalizing intensity and folder output option), masks, integrates and subtracts data"""
    parser = argparse.ArgumentParser(description='This takes two folders as input an option --Imono allows for using monochrometer values for normalization')
    parser.add_argument('SampFold', help = 'Need Sample Folder')
    parser.add_argument('BuffFold', help = 'Need Buffer Folder')
    parser.add_argument("--iend", action = "store_true", help = "Define Source iend or imono")
    parser.add_argument("--imono", action = "store_true", help = "Define Source iend or imono")
    parser.add_argument("--full", action = "store_true", help = "Keeps temporary files open like .sub and .msk")
    parser.add_argument("--BeforeOrAfter", type= str, required=False, help = "Pick if Buffer1 or Buffer2")
    args = parser.parse_args()
    SourceFlag = 1
    SourceWord = "iend"
    if args.iend == True:
        SourceFlag = 1
        SourceWord = "iend"
    if args.imono == True:
        SourceFlag = 2
        SourceWord = "imono"
    FullFlag = 0
    if args.full == True:
        FullFlag = 1
    BeforeOrAfterFlag = "b"
    if args.BeforeOrAfter:
        BeforeOrAfterFlag = args.BeforeOrAfter			
    MAINACTION = try_to_multi(args.SampFold,args.BuffFold,BeforeOrAfterFlag,SourceFlag,FullFlag)