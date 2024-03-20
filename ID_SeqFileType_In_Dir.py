#!/usr/bin/env python
import os
import sys
import argparse


#This function looks through the directory for anything with a .Extention extension and assumes this is the input
def ID_files(DirName,Extention): 

	directory_content_list = os.listdir(DirName)
	The_first_list = []
	for list_item in directory_content_list:
		if (str(list_item).endswith(Extention)):
			The_first_list.append(os.path.join(DirName, list_item))
	
	NumList = len(The_first_list)
	The_second_list = []
	for n in range(0,NumList):
		FileNameSplit = The_first_list[n].split(Extention)
		EndString = len(FileNameSplit[0])
		ThingToCheck = (FileNameSplit[0][-5:EndString])
		if (ThingToCheck.isdigit()):
			The_second_list.append(The_first_list[n])
	
	The_second_list.sort()
	#print (The_second_list)
	return (The_second_list)

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='Need A Directory Name And Extention to Look For')
	parser.add_argument('DirName', help = 'Need Sample Folder')
	parser.add_argument('Extention', help = 'Need Sample Folder')
	args = parser.parse_args()
	MAINACTION = ID_files(args.DirName,args.Extention)
