import sys
import os

APP_PATH  =  "./../src/bscite.exe" # path to B-SCITE executable
OUTPUT_FOLDER = "./example"
SUMMARY_FILE_FOLDER = "./"

SCFile   = "./SCFile-n_50-m_100.txt" # path input SC File
bulkFile = "./bulkFile-n_50.txt"     # path input bulk file
fp = 0.00001 # esimated false positive rate of SCS experiment
fn = 0.15    # estimated false negative rate of SCS experiment
n  = 50      # number of mutations
m  = 100     # number of cells
r  = 1       # number of repeats
l  = 2000    # number of loops

run_command =""
run_command += APP_PATH + " "
run_command += "-SCFileLocation "   + SCFile   + " "
run_command += "-bulkFileLocation " + bulkFile + " "
run_command += "-n "  + str(n)  + " "
run_command += "-m "  + str(m)  + " "
run_command += "-fd " + str(fp) + " "
run_command += "-ad " + str(fn) + " "
run_command += "-r "  + str(r)  + " "
run_command += "-l "  + str(l)  + " "
run_command += "-o "  + OUTPUT_FOLDER + " "
#run_command += "-summaryFolder " + SUMMARY_FILE_FOLDER.rstrip("/")
os.system(run_command)
