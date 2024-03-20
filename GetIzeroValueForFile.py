#!/usr/bin/env python
import os
import sys
import argparse
import numpy as np

import ProcessIzeroFile
import ID_SeqFileType_In_Dir

#def GetNormValue(FileName, Threshold, SourceFlag):
def GetNormValue(FileName, Threshold):

	if not os.path.isfile(FileName):
		sys.exit("%s Does Not Exist" %  FileName)
	WhereAmI = os.path.dirname(FileName)
	if not WhereAmI:
		WhereAmI =  "./"
	directory_content_list = os.listdir(WhereAmI)
	CbfFilesInDir = ID_SeqFileType_In_Dir.ID_files(WhereAmI, ".cbf")
	NumFiles = len(CbfFilesInDir)
	AlreadyIzeroAve = 0
	for list_item in directory_content_list:
		if (str(list_item).endswith("_Izero.ave")):
			#print(list_item)
			WorkWithThisFile = list_item
			#open(list_item,'r') as OpenWork
			FullNameWork = os.path.join(WhereAmI,list_item)
			WorkWithMe = np.loadtxt(FullNameWork)
			AlreadyIzeroAve = 1

			
	if (AlreadyIzeroAve == 0):
		IzeroCount = 0
		for list_item in directory_content_list:
			if (str(list_item).endswith("_Izero.txt")):
				IzeroFile = os.path.join(WhereAmI, list_item)
				IzeroCount += 1
		if IzeroCount == 0:
			sys.exit("No Izero File in %s" % WhereAmI)
		if IzeroCount > 1:
			sys.exit("Multiple Izero Files in %s" % WhereAmI)		
		RunCommand = ProcessIzeroFile.MakeNormList(IzeroFile,NumFiles,Threshold)
		#print(RunCommand)
		WorkWithMe = np.loadtxt(RunCommand)
#		print(RunCommand)
#		for list_item in directory_content_list:
#			if (str(list_item).endswith("_Izero.ave")):
#				WorkWithThisFile = list_item
#				#open(list_item,'r') as OpenWork
#				FullNameWork = os.path.join(WhereAmI,list_item)
#				WorkWithMe = np.loadtxt(FullNameWork)
#				AlreadyIzeroAve = 1
		
	#print(WorkWithMe)
	#print("LifeIsGood")
	
		
	FileNumber = int((FileName[-9:-4])) - 1
	#print("America")
#	WorkingArray = np.array(WorkWithMe)
	#print(np.shape(WorkWithMe))
	IendNormValue = WorkWithMe[FileNumber][1]
	
	ImonoNormValue = WorkWithMe[FileNumber][2]
	#print(FileName,IendNormValue,ImonoNormValue)	
	return (IendNormValue,ImonoNormValue)		
			
if __name__=="__main__":
	parser = argparse.ArgumentParser(description='This takes a file name, a threshold, and a flag and finds NormConstant')
	parser.add_argument('FileName', help = 'Need Sample Folder')
	parser.add_argument('Threshold', help = 'Need A Threshold')
	args = parser.parse_args()
	MAINACTION = GetNormValue(args.FileName, args.Threshold)
