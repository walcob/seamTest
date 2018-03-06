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
    return basename
    
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
        if(args.debug):
            cleanup(basename)
    
if __name__ == "__main__": main()