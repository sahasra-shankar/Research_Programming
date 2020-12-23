import sequenceAnalysis
sequenceAnalysis.NucParams.rnaCodonTable
from collections import defaultdict, Counter #tools to help with condensing dictionaries
def main ():
    """Calls methods from sequenceAnalysis and uses returned dictionaries to create compiled dictionaries with stats describing total codon counts 
    and frequencies in alphabetical order by amino acids in alphabetical order as well as calculating sequence length and GC content
    and prints with proper formatting"""
    myReader = sequenceAnalysis.FastAreader() #make sure to change this to use stdin
    masterDict = defaultdict(int)
    for head, seq in myReader.readFasta() : #calls all methods from sequenceAnalysis
        myNuc = sequenceAnalysis.NucParams(seq)
        myNuc.addSequence()
        myNuc.aaComposition()
        myNuc.nucComposition()
        fastaDict = myNuc.codonComposition()
        if not bool(masterDict):  #checks if dictionary is empty and places fastaDict in it if yes
            masterDict = fastaDict
        else:
            masterDict = Counter(masterDict) + Counter(fastaDict) #combines dictionaries if masterDict not empty
        #masterDict is dictionary with key being codon and value being how many times codon occurs in input
        
    sequenceLength = 0 
    gcContentVal =0  #sets initial value to 0 in order to add on it
    for aa in masterDict.keys():
        sequenceLength += len(aa)*masterDict[aa] #length of codon is multiplied with number of times codon occurs
        p=0
        while p < len(aa):
            if aa[p] in ['G','C']: #searches for occurrences of G and C in codons in masterDict
                gcContentVal += masterDict[aa] #counts GC occurrences
            p += 1 #adds 1 to p at the end of loop to continue iteration through codon
            
    codonStats = dict() #new dictionary for codon stats
    for nuc in masterDict.keys():
        aa = sequenceAnalysis.NucParams.rnaCodonTable
        codonStats[nuc] = aa[nuc] #matches key in masterDict with value in rnaCodonTable

    codonStatsSorted = sorted(codonStats.items() ,  key=lambda x: (x[1],x[0]))  #sort codons in alpha order, by amino acid
    #https://www.geeksforgeeks.org/ways-sort-list-dictionaries-values-python-using-lambda-function/
            
    print("Sequence Length {:.2f} Mb".format(round((sequenceLength/(1000000)),2))) #prints sequence length rounded to 2 decimals
    print("") #prints empty line
    print("GC content = {:.1f}%".format(100*(gcContentVal/sequenceLength))) #calculates and prints GC content rounded to 1 decimal
    print("")
    
    for nuc in codonStatsSorted: #calculate relative codon usage for each codon and print
        filtCodon = [key for key, val in codonStats.items() if val ==codonStats[nuc[0]]] #list of codons coding for same amino acid
        codonVals = [val for key, val in masterDict.items() if key in filtCodon] #list of number of occurrences of that codon
        codonFreqCount = sum(codonVals) #sums total occurrences for codons coding for 1 amino acid
        print ('{:s} : {:s} {:5.1f} ({:6d})'.format(nuc[0], nuc[1], (100*masterDict[nuc[0]]/codonFreqCount), masterDict[nuc[0]]))
        #prints codon, amino acid, codon freq, and number of occurrences for each codon

if __name__ == "__main__":
    main()