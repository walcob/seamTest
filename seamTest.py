#/bin/env python
import subprocess
import sys
import argparse
import glob

def cleanup(pdb):
    subprocess.run(["rm","-v","%s_1.pdb"%(pdb),"%s.hb"%(pdb)],stdout=sys.stdout)

def main():
    # Get list of pdbs
    barrels = glob.glob("database/barrels/*")
    nonBarrels = glob.glob("database/nonbarrels/*")
    pdbs = barrels+nonBarrels
    for pdb in pdbs:
        # parameters
        
        # renumber
        subprocess.run(["GeoFold/xrenumber_one",pdb,"%s_1.pdb"%(pdb)],stdout=sys.stdout)
        # pdb2hb
        hbfile = open("%s.hb"%(pdb),"w+")
        subprocess.run(["GeoFold/xpdb2hb","default.par","%s_1.pdb"%(pdb)],stdout=hbfile,stderr=sys.stdout)
        hbfile.close()
        # seams
        seamfile = open("%s.seams"%(pdb),"w+")
        subprocess.run(["GeoFold/seams/xpdb2seams2","%s.hb"%(pdb)],stdout=seamfile,stderr=sys.stdout)
        seamfile.close()
        # cleanup
        cleanup(pdb)