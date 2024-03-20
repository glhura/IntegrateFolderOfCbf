#!/usr/bin/env python

#This Code Takes an Izero File as input and an interger which describes how many 
#Dectris Frames were taken. It outputs a new file containing a list of index,intensity
#for each interval a frame was taken 

import sys
import argparse
import numpy as np
import os

import ChangeFileName


class CtimeImonoIend:
	def __init__(self, Ctime, Imono, Iend):
		self.Ctime = Ctime
		self.Imono = Imono
		self.Iend = Iend
	def __repr__(self):
		return repr((self.Ctime, self.Imono, self.Iend))
	 

def NumberOfLines(IzeroFileName):
	IzeroData = open(IzeroFileName, 'r')
	numlines = 0
	for line in IzeroData:
		numlines = numlines+1
	#print(numlines)
	IzeroData.close()
	return numlines

def DefineEdges(Data,Threshold):
	AmountData=len(Data)
	Threshold = float(Threshold)
	FrontEdge = 1
	BackEdge = AmountData
	for n in range(0,AmountData):
		if (float(Data[n]) > Threshold):
			FrontEdge = n
			break
	for n in range(FrontEdge,AmountData):
		if (float(Data[n]) < Threshold):
			BackEdge = n
			break
	#print(FrontEdge,BackEdge)
	return (FrontEdge,BackEdge)
	

		

def ReportResults(Data,NumberOfFrames,FrontEdge,BackEdge):
	TotalLength = BackEdge - FrontEdge
	#print(TotalLength)
	MakeInt = int(NumberOfFrames)
	#print(MakeInt)
	NumPointsPerFrame = int((TotalLength/MakeInt))
	#print(NumPointsPerFrame)
	StuffIwant = np.zeros((MakeInt,2))		
	k=0
	#print(np.average((Data[4948:4958].astype(np.float))))
	for j in range(0,(MakeInt)):
		StuffIwant[j,0] = int(j+1)
		FrontOfPacket = int(FrontEdge + (j*NumPointsPerFrame))
		BackOfPacket = int(FrontEdge + ((j+1)*NumPointsPerFrame))
		StuffIwant[j,1]= np.average(Data[FrontOfPacket:BackOfPacket].astype(float))
		#print(StuffIwant[j,0],StuffIwant[j,1])	
	#print(StuffIwant)
	return StuffIwant
		
			
	
#def MakeNormList(IzeroFileName,NumberOfFrames,ThresholdValue,SourceFlag):
def MakeNormList(IzeroFileName,NumberOfFrames,ThresholdValue):	
	HowManyLines = NumberOfLines(IzeroFileName)
	Column4 = []
	Column5 = []
	Column6 = []
	#print(HowManyLines)
	IzeroData = open(IzeroFileName, 'r')
	for n in range(1,HowManyLines):
		line_contents = IzeroData.readline()
		line_contents = line_contents.rstrip("\n") 
		line_contents= line_contents.split(" ")
		#print((line_contents[6]))
		if (len(line_contents) >= 7 and len(line_contents) <= 8):
			Column4.append(line_contents[4])
			Column5.append(line_contents[5])
			Column6.append(line_contents[6])
	Ctime_column = np.asarray(Column4)
	Imono_column = np.asarray(Column5)
	Iend_column = np.asarray(Column6)
	TheData = CtimeImonoIend(Ctime_column, Imono_column, Iend_column)
	#if SourceFlag == 1:
	IendSignal = TheData.Iend
	#print(len(IendSignal))
	#if SourceFlag == 2:
	ImonoSignal = TheData.Imono
	
	#print(type(TheData.Ctime))
	#if SourceFlag == 1:
	IendWhereFrontAndBack = DefineEdges(IendSignal,ThresholdValue)
	#print(IendWhereFrontAndBack)
	#if SourceFlag == 2:
	ImonoWhereFrontAndBack = (1,HowManyLines)
	IendFrontEdge = IendWhereFrontAndBack[0]
	IendBackEdge = IendWhereFrontAndBack[1]
	ImonoFrontEdge = ImonoWhereFrontAndBack[0]
	ImonoBackEdge = ImonoWhereFrontAndBack[1]	

	
	IendArrayToPrint = ReportResults(IendSignal,NumberOfFrames,IendFrontEdge,IendBackEdge)
	ImonoArrayToPrint = ReportResults(ImonoSignal,NumberOfFrames,ImonoFrontEdge,ImonoBackEdge)
	#print(IendArrayToPrint)
	MakeFirstInt = list((IendArrayToPrint[:,0].astype(int)))
	
	MakeSecondFloat = list(IendArrayToPrint[:,1].astype(float))
	#print(MakeSecondFloat)
	MakeThirdFloat = list(ImonoArrayToPrint[:,1].astype(float))
	NewFile1 = ChangeFileName.ChangeName(IzeroFileName,"ave",0)
	NewFile2 = open(NewFile1,"w")
	for n in range(0,int(NumberOfFrames)):
		NewFile2.write("%i %lf %lf\n" % (MakeFirstInt[n], MakeSecondFloat[n], MakeThirdFloat[n]))
		
	IzeroData.close()
	NewFile2.close()
	return(NewFile1)
	
if __name__=="__main__":
	parser = argparse.ArgumentParser(description='This takes a cbf file as input')
	parser.add_argument('IzeroFileName', help = 'Need Izero File')
	parser.add_argument('NumberOfFrames', help = 'Need An Integer Number of Frames')
	parser.add_argument('ThresholdValue', help = 'Threshold Value')
	#parser.add_argument("--iend", action = "store_true", help = "Define Source iend or imono")
	#parser.add_argument("--imono", action = "store_true", help = "Define Source iend or imono")
	args = parser.parse_args()
	#SourceFlag = 1
	#if args.iend == True:
		#SourceFlag = 1
	#if args.imono == True:
		#SourceFlag = 2	
	#print(SourceFlag)
	MAINACTION = MakeNormList(args.IzeroFileName,args.NumberOfFrames, args.ThresholdValue)
	
	#MAINACTION = MakeNormList(args.IzeroFileName,args.NumberOfFrames, args.ThresholdValue, SourceFlag)


