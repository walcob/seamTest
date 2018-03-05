#/bin/env python
import subprocess
import sys
import argparse
import glob

def cleanup():
    None
    
def writeParamaters(parfile):
    None

def main():
    # Get list of pdbs
    barrels = glob.glob("database/barrels/*")
    nonBarrels = glob.glob("database/nonbarrels/*")
    # getchain
    # renumber
    # 3to1
    # pdb2cij
    # pdb2hb
    # contactmask
    # seams
    # cleanup
    None