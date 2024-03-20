#!/usr/bin/env python

import sys
import argparse
import os

#Purpose of this program is to manipulate an incoming file name and change it so that it 
#has a new extension or will keep original one.
#Written with the intended purpose of changing .cbf files into .sub.cbf files 

def ChangeName(FileName,Mod,KeepOrigFlag):
	#print(FileName,Mod,KeepOrigFlag)
	KeepOrigFlag = int(KeepOrigFlag)
	NameLength = len(FileName)
	str1 = "."
	LastDot= FileName.rfind(str1)
	#print(LastDot)
	WhereToCut=NameLength-LastDot
	Trouble = FileName[(NameLength-WhereToCut):(NameLength)]
	#print(Trouble)
	FirstFileRoot = FileName[:-4]
	#CurrentFileExtension
	if KeepOrigFlag == 0:
		NewFileName = (FirstFileRoot+"."+Mod)
	if KeepOrigFlag == 1:
		NewFileName = (FirstFileRoot+"."+Mod+Trouble)
	#print(NewFileName)
	return NewFileName

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='This takes the name of a file as input')
	parser.add_argument("FileName", help = 'Need File to Work On')
	parser.add_argument("Mod", help = 'Need Extension You Want to Add')
	parser.add_argument("--KeepOrig", action = "store_true", help = "KeepsLast4")
	args = parser.parse_args()
	KeepOrigFlag = 0
	if args.KeepOrig:
		KeepOrigFlag = 1
	MAINACTION = ChangeName(args.FileName, args.Mod, KeepOrigFlag)


