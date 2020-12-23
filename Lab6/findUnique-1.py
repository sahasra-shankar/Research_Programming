#!/usr/bin/env python3
# Name: Sahasra Shankar (sshanka8)
# Group Members: Shweta Jones and Pranav Muthuraman (shsujones and pmuthura)
from collections import defaultdict
import sys
class FastAreader:

    def __init__(self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname

    def doOpen(self):
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):

        header = ''
        sequence = ''

        with self.doOpen() as fileH:

            header = ''
            sequence = ''

            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>'):
                if not line:  # we are at EOF
                    return header, sequence
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header, sequence
########################################################################
# Main
# Here is the main program
########################################################################
from functools import reduce
class findUnique:
    """Takes in input file and finds unique and non-essential substrings for each tRNA sequence to
    determine essential subsequences and list those under the corresponding section in the original sequence"""
    def __init__(self, filename):
        """Initializes essential and unique dicts to self and computes powersets of input after
        filtering out unnecessary characters from headers and sequences"""
        self.origDict = defaultdict(set) #initializes as a default dict of sets
        self.uniqueDict = defaultdict(set) #initializes as a default dict of sets
        self.finalDict = {}  #dict for essentials
        self.origSeqDict = {} #dict for original seqs from input
        for (header, sequence) in FastAreader(filename).readFasta(): #for header and seq in input
                header = header.replace(" ", "") #removes whitespace in header
                seqEdit = sequence.replace('-', '').replace('_', '').replace('.', '') #removes unneeded chars
                subSets = set(seqEdit[i:j + 1] for i in range(0, len(seqEdit)) for j in range(i, len(seqEdit) + 1)) #gets powerset
                self.origSeqDict[header] = sequence.replace('-', '').replace('_', '').replace('.', '') #removes unneeded chars from seq and stores to origSeqDict
                self.origDict[header] = subSets #stores powerset in origDict with key as header

    def checkDicts(self):
        """Uses origDict to create set of other tRNA subsequences and perform a union on those to create
        a dictionary of unique subsequences"""
        for tRnaMe in self.origDict.keys(): #compare key 1(header)
            others = set() #creates set for others
            for tRna in self.origDict.keys(): #compare key 2(header)
                if tRna != tRnaMe: #if not key not same
                    others = others.union(self.origDict[tRna]) #set union
            self.uniqueDict[tRnaMe] = self.origDict[tRnaMe] - others #unique dict has unique subseqs as value and

    def checkEssentials(self):
        """Uses uniqueDict from checkDicts method to filter out essential subsequences from unique and sort to
        find position of essential subs in original seq for formatting"""
        self.checkDicts()
        for tRna in self.uniqueDict.keys(): #for key(header) in uniqueDict
            essential_set = set() #create essential set
            trna_set = sorted(self.uniqueDict[tRna], key=len) #sort by key length
            non_essential = set()  # create nonessential set
            essential_dict = defaultdict() #sets essential_dict as default dict
            sorted_essential_dict = defaultdict() #sets sorted_essential_dict to default
            for i in range(0,len(trna_set)): #for the range of set length
                for j in range(i+1, len(trna_set)): #compare next in set
                    if trna_set[i] in trna_set[j]: #if subseq found in next subseq
                        non_essential.add(trna_set[j])#put next subseq in nonessential
            essential_set = set(trna_set) - non_essential #essential is set of unique- set of nonessential
            for essential_seq in essential_set: #for each entry in essential set
                orig_seq = self.origSeqDict[tRna] #set orig_seq to edited seq
                subseqPos_in_trnaSeq = orig_seq.find(essential_seq) #finds position of essential subseq within orig seq
                essential_dict[essential_seq] = subseqPos_in_trnaSeq #sets pos as value in essential_dict
            sorted_essential_dict = dict(sorted(essential_dict.items(), key=lambda x: x[1])) #sorts values in essential dict to put into sorted_essential_dict
            self.finalDict[tRna] = sorted_essential_dict #sets finaldict values to sorted_essential_dict

    def final_output(self):
        """Uses finalDict from checkEssentials method to format essential subsequences under each
        original tRNA sequence using . as alignment guides"""
        self.checkEssentials() #gets dict from checkEssentials method
        for tRna, ess_dict in self.finalDict.items(): #for each header and essential dict in finalDict
            print(tRna) #print header
            print(self.origSeqDict[tRna]) #print original seq
            for subseq in ess_dict.keys(): #for each subseq in the keys of the essential dict within finalDict
                output_str = '.' * ess_dict[subseq] + subseq #output essential with posvalue-1(.)
                print(output_str) #print essential with this format

def main(inCL=None):
    """Sends input file to findUnique class and then calls final_output method to print
    output with correct formatting"""
    find_unique = findUnique() #uses whatever input file is specified HERE
    find_unique.final_output() #calls final_output method from findUnique class

if __name__ == "__main__":
    main()
