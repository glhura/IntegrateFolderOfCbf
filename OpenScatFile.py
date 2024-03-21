#!/usr/bin/env python
import numpy as np

class Q_I_Error:
	def __init__(self, q_values, I_values, E_values):
		self.q_values = q_values
		self.I_values = I_values
		self.E_values = E_values
	def __repr__(self):
		return repr((self.q_values, self.I_values, self.E_values))

#This function reads a file and outputs a class Q_I_Error
#If Error is missing it adds Error of 3% of I

def OpenScatFile(file_name):
	# Read the data in from a file to a list
	Column1 = []
	Column2 = []
	Column3 = []
	num_lines = 0
	
	with open(file_name, 'r') as f:
		for line in f:
			if line.find('#') > 0:
				line = f.readline()
			line = line.rstrip("\n")
			line = line.rstrip("\t")
			line = line.rstrip(" ")		
			line_contents= line.split(" ")
		
			if (len(line_contents) == 3 and float(line_contents[0])):

				Column1.append(line_contents[0])
				Column2.append(line_contents[1])
				Column3.append(line_contents[2])
		
			if (len(line_contents) == 2 and float(line_contents[0])):	
				Column1.append(line_contents[0])
				Column2.append(line_contents[1])
				Column3.append(float(line_contents[1])*0.03) 
    	
	Q_column = np.asarray(Column1).astype(float)
	I_column = np.asarray(Column2).astype(float)
	E_column = np.asarray(Column3).astype(float)
    
	FinalOutput = Q_I_Error(Q_column, I_column, E_column)
	f.close()
	return FinalOutput
