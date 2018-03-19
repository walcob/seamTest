#/bin/env python
import subprocess
import sys
import argparse
import glob
import os.path

def cleanup(basename):
    subprocess.run(["rm","-v","%s_1.pdb"%(basename),"%s.hb"%(basename)],stdout=sys.stdout)

def runTest(pdb):
    basename = os.path.basename(pdb)[:-4]
    print("="*25,basename,"="*25)
    # renumber
    print("="*25,"Renumbering","="*25)
    subprocess.run(["GeoFold/xrenumber_one",pdb,"%s_1.pdb"%(basename)],stdout=sys.stdout)
    # pdb2hb
    print("="*25,"Extracting H-bonds","="*25)
    hbfile = open("%s.hb"%(basename),"w+")
    subprocess.run(["GeoFold/xpdb2hb","default.par","%s_1.pdb"%(basename)],stdout=hbfile,stderr=sys.stdout)
    hbfile.close()
    # seams
    print("="*25,"Finding Seams","="*25)
    seamfile = open("seamfiles/%s.seams"%(basename),"w+")
    subprocess.run(["GeoFold/seams/xpdb2seams2","%s.hb"%(basename)],stdout=seamfile,stderr=sys.stdout)
    seamfile.close()
    return basename

def makepar(pdb,parfile):
    fout = open(parfile,'w+')
    fout.write("LNAME %s\n"%(os.path.basename(pdb)))
    fout.write("EMAIL walcob@rpi.edu\n")
    fout.write("PDBCODE %s\n"%(os.path.basename(pdb)[:-4]))
    fout.write("OMEGA 1.\nINTERMEDIATE 0\nBARRELMOVES 1\nORANGE 0. 0.2 0.4 0.6 0.8 1.0 1.2 1.4\nRUNGEOFOLD 1\nMOLSCRIPT 0\nBREAKCUT 0.05\nPIVOTCUT 0.01\nHINGECUT 0.5\nSEAMCUT 10\nBREAKPOINTENTROPY 90.\nHINGEPOINTENTROPY 30.\nTEMPERATURE 300.\nCONCENTRATION 1.\nVOIDENTROPY 0.\nSOLIDITY 1000.\nHBONDENERGY 100.\nHAMMONDSCALE 1000.\nSIDECHAINENTROPY 1.\nHINGEBARRIER 0.\nPIVOTBARRIER 0.\nWATER 1.\nMAXSPLIT 4\nMAXTIME 10.\nMINSEG 4\nCAVITATION 0.000001\nFLORY 0\nSUBMIT submit\nREDUCING 0\nMOLSCRIPT 1\nHLFE 0\nFING 0\nCHAIN +\nwc 1.\nwn 8.\n")
    fout.close()

def fullTest(pdb,parfile):
    makepar(pdb,parfile)
    subprocess.run(["python","GeoFold/rungeofold.py",parfile,"/Users/walcob/seamTest/seamTest.conf"],stdout=sys.stdout,stderr=sys.stdout)
    
    
def main():
    # Get list of pdbs
    parser = argparse.ArgumentParser(description="Test seams")
    parser.add_argument("--all",action="store_true",default=False)
    parser.add_argument("--debug",action="store_true")
    parser.add_argument("-f",nargs='+')
    args = parser.parse_args()
    pdbs = []
    if(args.all):
        pdbs += glob.glob("database/barrels/*") + glob.glob("database/nonbarrels/*")
    if(args.f is not None): pdbs += args.f
    for pdb in pdbs:
        basename = runTest(pdb)
        # cleanup
        if(not args.debug):
            cleanup(basename)
    
if __name__ == "__main__": main()