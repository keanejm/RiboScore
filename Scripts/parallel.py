from os import listdir, fork, waitpid, _exit, EX_OK, WNOHANG
import os
from sys import exit, stdout
from time import sleep
from rust_bam import writerustfiles
from gbd_bam import gini, gbd_gini

def removeIfDone(pidList, startedCount, finishedCount):
    for pid in pidList:
        sleep(1)                        # Reduce processor load 
        returnedPid, status = waitpid(pid, WNOHANG)
        if returnedPid > 0:
            pidList.remove(returnedPid)
            finishedCount += 1

            # Carriage return to write over previous line
            stdout.write("Processes: started %i, finished %i      \r" % (startedCount, finishedCount))
            stdout.flush()

    return pidList, finishedCount

def rustforker(seqDict, pathToRefgene, pathToStudy, offset, readlen_range, outputDir, maxProcesses):    
    pidList = []
    startedCount = 0
    finishedCount = 0

    # Ensure one end slash
    if not pathToStudy[-1] == "/":
        pathToStudy += "/"
           
    for root, dirs, filenames in os.walk(pathToStudy):
	for treatmentDirContent in filenames:
	    if treatmentDirContent[-4:] == ".bam":
                filepath = os.path.abspath(os.path.join(root,treatmentDirContent))
                startedCount += 1

                # Carriage return to write over previous line
                stdout.write("Processes: started %i, finished %i      \r" % (startedCount, finishedCount))
                stdout.flush()
                childPid = fork()                               # Create new process and store id
                if childPid < 0:
                    exit("Error: failed to create child process")

                elif childPid == 0:                             # To be executed by a child process
                    writerustfiles(seqDict, pathToRefgene, filepath, offset[treatmentDirContent[:-4]], readlen_range, outputDir)
                    _exit(EX_OK)

                else:                                           # To be executed by the parent process
                    pidList.append(childPid)    
                    while len(pidList) >= maxProcesses:         # Limit number of processes
                        pidList, finishedCount = removeIfDone(pidList, startedCount, finishedCount)

    while len(pidList) > 0:
        pidList, finishedCount = removeIfDone(pidList, startedCount, finishedCount)

    print "Processes: started %i, finished %i      " % (startedCount, finishedCount)    

    
def gbdforker(pathToStudy, offset, tree_dict, maxProcesses):    
    pidList = []
    startedCount = 0
    finishedCount = 0

    # Ensure one end slash
    if not pathToStudy[-1] == "/":
        pathToStudy += "/"
           
    for root, dirs, filenames in os.walk(pathToStudy):
	for treatmentDirContent in filenames:
	    if treatmentDirContent[-4:] == ".bam":
                filepath = os.path.abspath(os.path.join(root,treatmentDirContent))
                startedCount += 1

                # Carriage return to write over previous line
                stdout.write("Processes: started %i, finished %i      \r" % (startedCount, finishedCount))
                stdout.flush()
                childPid = fork()                               # Create new process and store id
                if childPid < 0:
                    exit("Error: failed to create child process")

                elif childPid == 0:                             # To be executed by a child process
                    gbd_gini(filepath, offset[treatmentDirContent[:-4]], tree_dict)
                    _exit(EX_OK)

                else:                                           # To be executed by the parent process
                    pidList.append(childPid)    
                    while len(pidList) >= maxProcesses:         # Limit number of processes
                        pidList, finishedCount = removeIfDone(pidList, startedCount, finishedCount)

    while len(pidList) > 0:
        pidList, finishedCount = removeIfDone(pidList, startedCount, finishedCount)

    print "Processes: started %i, finished %i      " % (startedCount, finishedCount)      
    
    
    
    
    
    
    
    
    
    
