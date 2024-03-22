# IntegrateFolderOfCbf

pil_2m_integ_cbf.c integrates cbf images

Stands for pilatus 2m integrated cbf
To compile the c code:

gcc -o pil_2m_integ_cbf pil_2m_integ_cbf.c -lm

To run this code two other files need to be available in the "current working" directory 1) ExptParams and 2) X_Y_mask.asc

Meant to be included is a TestFiles directory.
An example of running this command is

./pil_2m_integ_cbf TestFiles/C2_C2_snap_00001.cbf

Should create TestFiles/C2_C2_snap_00001.dat

For this code to be used by python scripts a pil_2m_integ_cbf.so file needs to be created as follows

gcc -shared -o pil_2m_integ_cbf.so -fPIC pil_2m_integ_cbf.c

Once available, the python script called MIS_CBF_Foder.py will act on two folders
In the TestFiles directory the folders "C2" and "C1b" are examples

MIS_CBF_Foder.py TestFiles/C2 TestFiles/C1b

MIS stands for Mask, Integrate and Subtract

Should create folders:

Results/Subtracted/TestFiles/C2_results/Buffer1

Results/Unsubtracted/TestFiles/C1b_results

Results/Unsubtracted/TestFiles/C2_results

To run, these other python scripts are necessary:

ChangeFileName.py
GetIzeroValueForFile.py
ID_SeqFileType_In_Dir.py
OpenScatFile.py
ProcessIzeroFile.py
SubtractTwoDats.py



