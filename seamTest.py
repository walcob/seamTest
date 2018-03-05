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
    print(basename)
    print("="*10,pdb,"="*10)
    # renumber
    subprocess.run(["GeoFold/xrenumber_one",pdb,"%s_1.pdb"%(basename)],stdout=sys.stdout)
    # pdb2hb
    hbfile = open("%s.hb"%(basename),"w+")
    subprocess.run(["GeoFold/xpdb2hb","default.par","%s_1.pdb"%(basename)],stdout=hbfile,stderr=sys.stdout)
    hbfile.close()
    # seams
    seamfile = open("seamfiles/%s.seams"%(basename),"w+")
    subprocess.run(["GeoFold/seams/xpdb2seams2","%s.hb"%(basename)],stdout=seamfile,stderr=sys.stdout)
    seamfile.close()
    
def main():
    # Get list of pdbs
    barrels = glob.glob("database/barrels/*")
    nonBarrels = glob.glob("database/nonbarrels/*")
    pdbs = barrels+nonBarrels
    for pdb in pdbs:
        runTest(pdb)
        # cleanup
        cleanup(basename)
    
if __name__ == "__main__": main()