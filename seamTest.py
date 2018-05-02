#/bin/env python
"""
This python script runs test for the new seam implementation in GeoFold, pdb2seams2.
"""
import subprocess
import sys
import argparse
import glob
import os.path



def cleanup(basename,log):
    '''Clean up temporary files generated in a run of seamtest'''
    # Log can either be to sys.stdout or a file
    openedLog = False
    if log != sys.stdout and type(log) == str:
        log = open(log,'ab')
        openedLog = True
    # delete the files
    subprocess.run(["rm","-v","%s_1.pdb"%(basename),"%s.hb"%(basename)],stdout=log)
    if openedLog: log.close()

def runTest(pdb,log):
    '''This is the meat of the program.  Given a pdb file, performs preprocessing
    necessary for pdb2seams2 and performs pdb2seams2 barrel detection.'''
    # Deal with the logfile
    openedLog = False
    if log != sys.stdout and type(log) == str: 
        log = open(log,'ab')
        openedLog = True
    # get PDB code
    basename = os.path.basename(pdb)[:-4]
    # renumber, number the pdb sequentially
    subprocess.run(["GeoFold/xrenumber_one",pdb,"%s_1.pdb"%(basename)],stdout=log)
    # pdb2hb, find hydrogen bonds
    hbfile = open("%s.hb"%(basename),"w+")
    subprocess.run(["GeoFold/xpdb2hb","default.par","%s_1.pdb"%(basename)],stdout=hbfile,stderr=log)
    hbfile.close()
    # seams, perform pdb2seams2 barrel detection
    seamfile = open("seamfiles/%s.seams"%(basename),"w+")
    subprocess.run(["GeoFold/seams/xpdb2seams2","%s.hb"%(basename)],stdout=seamfile,stderr=log)
    seamfile.close()
    if openedLog: log.close()
    return basename

def makepar(pdb,parfile):
    '''Generates parameters file, parfile, for a full run of GeoFold on pdb'''
    fout = open(parfile,'w+')
    fout.write("LNAME %s\n"%(os.path.basename(pdb)))
    fout.write("EMAIL walcob@rpi.edu\n")
    fout.write("PDBCODE %s\n"%(os.path.basename(pdb)[:-4]))
    fout.write("OMEGA 1.\nINTERMEDIATE 0\nBARRELMOVES 1\nORANGE 0. 0.5 1.0 1.5\nRUNGEOFOLD 1\nMOLSCRIPT 0\nBREAKCUT 0.05\nPIVOTCUT 0.25\nHINGECUT 0.5\nSEAMCUT 10\nBREAKPOINTENTROPY 90.\nHINGEPOINTENTROPY 30.\nTEMPERATURE 300.\nCONCENTRATION 1.\nVOIDENTROPY 0.\nSOLIDITY 1000.\nHBONDENERGY 100.\nHAMMONDSCALE 1000.\nSIDECHAINENTROPY 1.\nHINGEBARRIER 0.\nPIVOTBARRIER 0.\nWATER 1.\nMAXSPLIT 4\nMAXTIME 10.\nMINSEG 4\nCAVITATION 0.000001\nFLORY 0\nSUBMIT submit\nREDUCING 0\nMOLSCRIPT 1\nHLFE 0\nFOLDING 0\nCHAIN +\nwc 1.\nwn 8.\n")
    fout.close()

def fullTest(pdb,parfile,log):
    '''Performs full run of GeoFold given parfile parameters file and pdb'''
    # logfile handling
    openedLog = False
    if log != sys.stdout and type(log) == str: 
        log = open(log,'ab')
        openedLog = True
    # create the parameters file
    makepar(pdb,parfile)
    # run GeoFold
    subprocess.run(["python","GeoFold/rungeofold.py",parfile,"/Users/walcob/seamTest/seamTest.conf"],stdout=log,stderr=log)
    if openedLog: log.close()
    
def getSCOPe(pdb):
    '''Extracts the SCOPe structural class from a given pdb'''
    with open(pdb) as pdbIn:
        lines = "\n".join([line for line in pdbIn if "REMARK  99 ASTRAL SCOPe-sccs:" in line])
    return lines.split(":")[1].strip()
    
def findBarrel(basename):
    '''Check pdb2seams2 output for barrels'''
    # Find the line in the output showing the number of barrels detected
    with open("seamfiles/%s.seams"%(basename)) as seamIn:
        line = "".join([line for line in seamIn if line[:8] == "NBARRELS"])
        print(line);sys.stdout.flush()
        try:
            return int(line.split()[1]) > 0
        except ValueError or IndexError:
            print("Error: ",basename,line);sys.stdout.flush()
            return False
    
def SCOPeTest():
    '''Detect barrels in entire SCOPe database'''
    # MPI setup to make this faster
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    if rank == 0: print("comm initialized");sys.stdout.flush()
    # initialize result output
    barrels = {}
    nonbarrels = {}
    # Read in SCOPe PDBs
    pdbs = sorted(glob.glob("/Users/walcob/SCOPe/*/*.ent"))
    if rank == 0: print("PDBs read in");sys.stdout.flush()
    # get SCOPe identifiers
    SCOPe = {}
    for pdb in pdbs:
        SCOPeStr = getSCOPe(pdb)
        SCOPe[pdb] = SCOPeStr
        if rank == 0: print(pdb,SCOPeStr,"       ",end='\r');sys.stdout.flush()
    if rank == 0: print("SCOPe read in");sys.stdout.flush()
    # run tests
    for i in range(rank,len(pdbs),size):
        print(pdbs[i]);sys.stdout.flush()
        basename = runTest(pdbs[i],log=sys.stdout)
        foundBarrel = findBarrel(basename)
        print(i,pdbs[i],foundBarrel)    ;sys.stdout.flush()
        # add results to correct dictionary
        if foundBarrel:
            try:
                barrels[SCOPe[pdb]].append(pdb)
            except KeyError:
                barrels[SCOPe[pdb]] = [pdb]
        else:
            try:
                nonbarrels[SCOPe[pdb]].append(pdb)
            except KeyError:
                nonbarrels[SCOPe[pdb]] = [pdb]
        # cleanup unneeded files
        cleanup(basename,log=sys.stdout)
    # write output
    writeSCOPe(barrels,"SCOPeBarrels_%i.txt"%(rank),nonbarrels,"SCOPeNonbarrels_%i.txt"%(rank))

def writeSCOPe(barrels,barrelsout,nonbarrels,nonbarrelsout):
    '''Write out SCOPe domains detected as barrels and nonbarrels'''
    bout = open(barrelsout,"w+")
    nout = open(nonbarrelsout,"w+")
    # for every SCOPe ID, write the id and all entries of that class identified as barrels
    for SCOPe in sorted(barrels.keys()):
        bout.write("%s\n"%(SCOPe))
        bout.write("\n".join(barrels[SCOPe]))
        bout.write("\n")
    bout.close()
    # like above, but for nonbarrels
    for SCOPe in sorted(nonbarrels.keys()):
        nout.write("%s\n"%(SCOPe))
        nout.write("\n".join(nonbarrels[SCOPe]))
        nout.write("\n")
    nout.close()
    
def doFullTests(log,pdbs=['database/barrels/2B3P.pdb','database/barrels/2AWJ.pdb','database/barrels/1I1B.pdb','database/barrels/1R2T.pdb','database/barrels/1JOO.pdb']):
    '''Perform full runs of GeoFold so pathway prediction accuracy may be assessed.
    Only perform this for entries with experimental folding data.'''
    # MPI setup to make this faster
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    # logfile parsing
    if log != sys.stdout: log = open("%s_%i.%s"%(log.split('.')[0],rank,log.split('.')[1]),'ab')
    # load pdbs and run GeoFold
    for i in range(rank,len(pdbs),size):
        print(i,rank,pdbs[i]);sys.stdout.flush()
        fullTest(pdbs[i],"%s.par"%(os.path.basename(pdbs[i])[:-4]),log)
    # close log
    if log != sys.stdout: log.close()
    
def main():
    # parser setup
    parser = argparse.ArgumentParser(description="Test seams")
    # run on all pdbs in database directory
    parser.add_argument("--all",action="store_true",default=False)
    # use debug output
    parser.add_argument("--debug",action="store_true")
    # files to test
    parser.add_argument("-f",nargs='+')
    # run on the entire SCOPe database
    parser.add_argument("--SCOPe",action="store_true",default=False)
    # do a full run of GeoFold
    parser.add_argument("--full",action="store_true",default=False)
    # logfile
    parser.add_argument("-log",default=None,help="log output to logfile. set to logfile path")
    
    # parse arguments
    args = parser.parse_args()
    pdbs = []
    # logfile setup
    if args.log:
        log = args.log
    else:
        log = sys.stdout
    # test full SCOPe
    if args.SCOPe:
        SCOPeTest()
    # perform full GeoFold test
    elif args.full:
        if args.f is None:
            doFullTests(log)
        else:
            doFullTests(log,args.f)
    else:
        # check full database
        if(args.all):
            pdbs += glob.glob("database/barrels/*") + glob.glob("database/nonbarrels/*")
        # check specified files
        if(args.f is not None): pdbs += args.f
        # run the tests
        for pdb in pdbs:
            basename = runTest(pdb,log)
            # cleanup
            if(not args.debug):
                cleanup(basename,log)
            
if __name__ == "__main__": main()