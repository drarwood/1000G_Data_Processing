'''
This python script takes in a list of variants and extracts allele freq differences
across the 5 broad genetic ancestries, and across all 26 sub-ancestries
Written by Andrew wood
Question to A.R.Wood@exeter.ac.uk
'''
import argparse
import struct
import numpy as np

def LoadPopulationPairs(a):
    '''
    This function creates a dictionary of index (int) -> population pair label
    Parameters:
      a (string) : file prefix associated with 1000G variant frequency difference data
    Returns:
      x (dict.) : dictionary of variant_id -> KG_Variant class objects
    '''
    x = {}
    with open(a+".pp") as f:
        for l in f:
            d=l.rstrip("\n").split("\t")
            x[d[0]] = d[1]
    return x


class KGVariant:
    def __init__(self, a1, a2, start_byte, non_zero_diffs):
        self.a1 = a1
        self.a2 = a2
        self.start_byte = int(start_byte)
        self.non_zero_diffs = int(non_zero_diffs)

def LoadIndexFile(a):
    '''
    This function creates a dictionary of class objects for 1000G variants
    Parameters:
      a (string) : file prefix associated with 1000G variant frequency difference data
    Returns:
      x (dict.) : dictionary of variant_id -> KGVariant class objects
    '''
    x = {}
    with open(a+".va") as f:
        for l in f:
            d=l.rstrip("\n").split("\t")
            x[d[0]] = KGVariant(d[1], d[2], d[3], d[4])
    return x


class UserVariant:
    def __init__(self, primary_allele, other_allele, primary_allele_match):
        self.primary_allele = primary_allele
        self.other_allele   = other_allele
        self.primary_allele_match = primary_allele_match

def LoadUserVars(a,b):
    '''
    This function creates a dictionary of class objects for user-defined variants for lookup
    Parameters:
      a (string): filename with user-defined variants for lookup
      b (dict.) : kg_variant_id -> KGVariant class objects
    Returns:
     x (dict.): dictionary of variant_id -> UserVariant class objects
    '''
    x = {}
    with open(a) as f:
        for l in f:
            d=l.rstrip("\n").split("\t")
            # set default id to be used for lookup in 1000G data:
            id = ":".join(d)
            if id in b or ":".join([d[0],d[1],d[3],d[2]]) in b:
                # update ID to be used for lookup in 1000G data if required:
                if ":".join([d[0],d[1],d[3],d[2]]) in b:
                    id = ":".join([d[0],d[1],d[3],d[2]])
                # Determine if the allele that freq diffs is based on == primary user-defined allele (column 3 in user list)
                primary_allele_match = False
                if d[2] == b[id].a1:
                   primary_allele_match = True
                # add variant to dictionary of UserVariant class objects
                x[id] = UserVariant(d[2], d[3], primary_allele_match)
    return x



def Extract1000GAlleleFreqDiffs(a, b, c, d, e):
    '''
    Parameters:
      a (string) : file prefix associated with 1000G variant frequency difference data
      b (string) : filename of output file
      c (dict.)  : index -> population pair label
      d (dict.)  : varid -> KGVariant class object
      e (dict.)  : varid -> UserVariant class object
    Returns:
      NA
    '''
    with open(a+".fd", "rb") as f, open(b, "w") as o:

        # Form header for output file: a1 = primary allele by user diffs refer to
        o.write("\t".join(["variant\ta1\ta2\t" + "\t".join(list(c.values()))]) +"\n")

        # cycle through e (dictionary containing user variant data for lookup)
        for k,v in e.items():
            print(k)
            # create numpy array of size len(c) initialised to zero
            diffs = np.zeros(len(c))
            # Get the start byte of the variant, and the number of non-zero diffs
            start_byte = d[k].start_byte
            non_zero_diffs = d[k].non_zero_diffs

            # only attempt to read block > 0
            if non_zero_diffs > 0:
                # Go to starting byte associated with variant
                f.seek(start_byte)
                buffer = f.read(6*non_zero_diffs)
                for fields in struct.iter_unpack("<Hf", buffer):
                    pop_index, diff = fields
                    diffs[pop_index] = round(diff, 4)

            # determine whether we need to multiply diffs by -1 if 
            # primary allele used for freq diff != primary allele for lookup defined by user
            # Action if required:
            if not v.primary_allele_match:
                diffs = (diffs * -1)+0.0

           # Write to file
                o.write("\t".join([k, v.primary_allele, v.other_allele]) + "\t" +  "\t".join(str(x) for x in diffs) + "\n")

if __name__ == "__main__":

    # Get command line arguments for parameters
    parser = argparse.ArgumentParser()
    parser.add_argument("--freq-diffs", required=True)
    parser.add_argument("--vars", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    # Load population pairs
    pop_pairs = LoadPopulationPairs(args.freq_diffs)

    # Load index file into variant class and return dictionary of variant class objects
    kg_vars = LoadIndexFile(args.freq_diffs)

    # Read user defined variant list and extract relevant data
    # Tab delimited: chr bp primary_allele other_allele
    user_vars = LoadUserVars(args.vars, kg_vars)

    # Extract allele frequency differences across population pairs
    Extract1000GAlleleFreqDiffs(args.freq_diffs, args.out, pop_pairs, kg_vars, user_vars)
