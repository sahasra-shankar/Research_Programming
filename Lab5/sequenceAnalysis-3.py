class ProteinParam:
    """For a given protein sequence, calculate physical-chemical properties such as the number of amino acids in the string,
    the molecular weight, molar and mass extinction coefficients, theoretical PI and the amino acid composition"""
    # These tables are for calculating:
    #     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
    #     absorbance at 280 nm (aa2abs280)
    #     pKa of positively charged Amino Acids (aa2chargePos)
    #     pKa of negatively charged Amino acids (aa2chargeNeg)
    #     and the constants aaNterm and aaCterm for pKa of the respective termini
    #  Feel free to move these to appropriate methods as you like

    # As written, these are accessed as class attributes, for example:
    # ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    aa2mw = {
        'A': 89.093, 'G': 75.067, 'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225, 'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
    }

    mwH2O = 18.015
    aa2abs280 = {'Y': 1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R': 12.4, 'H': 6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__(self, protein):
        """initialize self from input protein sequence and compute aaComposition dictionary"""
        self.protein = protein.upper()  # make input uppercase
        self.aaCompositionDict = dict()
        for aa in self.aa2mw.keys():
            aaCnt = 0
            for seq in list(self.protein):
                if aa == seq:
                    aaCnt += 1
            self.aaCompositionDict[aa] = aaCnt

    def aaCount(self):
        """return single int count of aa characters found in input string"""
        for aa in self.aa2mw.keys():
            result = 0
            for seq in list(self.protein):
                if aa == seq:
                    result += 1
                elif aa != seq:
                    result += 0
        return result

    def pI(self):
        """find pH producing neutral net charge when put into charge method and return"""
        bestCharge = 100000
        bestPh = 0
        for value in range(0, 1400 + 1):
            value = value / 100
            charge = abs(self._charge_(value))
            if charge < bestCharge:
                bestCharge = charge
                bestPh = value
        return bestPh

    def aaComposition(self):
        """return aaComposition dictionary computed in init"""
        return self.aaCompositionDict

    def _charge_(self, bestPh):
        """calculates individual charge of amino acids having charge in input string and returns net charge """
        for aa in self.aa2mw.keys():
            for seq in list(self.protein):
                if aa == seq:
                    rCharge = self.aaCompositionDict['R'] * (
                                10 ** self.aa2chargePos['R'] / ((10 ** self.aa2chargePos['R']) + (10 ** bestPh)))
                    kCharge = self.aaCompositionDict['K'] * (
                                10 ** self.aa2chargePos['K'] / ((10 ** self.aa2chargePos['K']) + (10 ** bestPh)))
                    hCharge = self.aaCompositionDict['H'] * (
                                10 ** self.aa2chargePos['H'] / ((10 ** self.aa2chargePos['H']) + (10 ** bestPh)))
                    dCharge = self.aaCompositionDict['D'] * (
                                10 ** bestPh / ((10 ** self.aa2chargeNeg['D']) + (10 ** bestPh)))
                    eCharge = self.aaCompositionDict['E'] * (
                                10 ** bestPh / ((10 ** self.aa2chargeNeg['E']) + (10 ** bestPh)))
                    cCharge = self.aaCompositionDict['C'] * (
                                10 ** bestPh / ((10 ** self.aa2chargeNeg['C']) + (10 ** bestPh)))
                    yCharge = self.aaCompositionDict['Y'] * (
                                10 ** bestPh / ((10 ** self.aa2chargeNeg['Y']) + (10 ** bestPh)))

                    pChargeTotal = rCharge + kCharge + hCharge + (
                                (10 ** self.aaNterm) / ((10 ** self.aaNterm) + (10 ** bestPh)))
                    nChargeTotal = dCharge + eCharge + cCharge + yCharge + (
                                (10 ** bestPh) / ((10 ** self.aaCterm) + (10 ** bestPh)))
                    netCharge = pChargeTotal - nChargeTotal
                elif aa != seq:  # assume if non-amino acid exists then net charge value does not exist
                    netCharge = 0
        return netCharge

    def molarExtinction(self):
        """calculates protein light absorption at wavelength of 280nm and returns total value"""
        yAbs = ((self.aa2abs280['Y']) * self.aaCompositionDict['Y'])
        wAbs = ((self.aa2abs280['W']) * self.aaCompositionDict['W'])
        bAbs = ((self.aa2abs280['C']) * self.aaCompositionDict['C'])
        total = yAbs + wAbs + bAbs
        return total

    def massExtinction(self):
        """calculates by dividing molar extinction coefficient by molecular weight and returns value"""
        myMW = self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0

    def molecularWeight(self):
        """returns molecular weight of protein sequence"""
        newList = []
        for aa in self.aa2mw.keys():
            for seq in list(self.protein):
                if aa == seq:
                    myAA = self.aa2mw[aa]
                    newList.append(myAA)
                    sumAA = sum(newList)  # sums all weights in new list
                    # https://appdividend.com/2019/08/09/python-sum-example-sum-function-in-python-tutorial/
                    molWeight = self.mwH2O + (sumAA - (self.aaCount() * self.mwH2O))
                elif aa != seq:  # assume if non-amino acid exists then molecular weight cannot be calculated
                    molWeight = 0
        return molWeight


class NucParams:
    """For all input, filters and returns dictionaries displaying stats of amino acid count, nucleotide base (singular)counts, codon counts and total base count while handling invalid codons/bases appropriately--information is received by genomeAnalyzer which prints final stats with formatting
    ex.
    input: testGenome.fa file
    output:
    sequence length= 3.14 Mb

    GC content= 60.2%

    UAA: -32.6 (1041)
    UAG: -38.6 (1230)
    UGA: -28.8 (918)
    GCA: A 14.1 (10605)
    GCC: A 40.5 (30524)
    GCG: A 30.5 (22991)
    GCU: A 14.9 (11238)
    UGC: C 67.2 (4653)
    .
    .
    .
    """
    rnaCodonTable = {
        # RNA codon table
        # U
        'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
        'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
        'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
        'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
        # C
        'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
        'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
        'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
        'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
        # A
        'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
        'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
        'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
        'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
        # G
        'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
        'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
        'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
        'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    dnaCodonTable = {key.replace('U', 'T'): value for key, value in rnaCodonTable.items()}

    def __init__(self, inSeq, inString=''):
        """Splits input by every three characters and puts into new list while accounting for any additional input that is added"""
        self.inString = inString
        self.inSeq = inSeq
        finalInput = inString + inSeq  # sets input to initial input plus anything that was added to it after
        n = 3  # https://stackoverflow.com/questions/9475241/split-string-every-nth-character
        codonList = [finalInput[i:i + n] for i in range(0, len(finalInput), n)]  # list of codons
        self.codonList = codonList  # initialized codonList to self

    def addSequence(self):
        """Returns codon list computed in init object because init cannot return lists"""
        return self.codonList

    def aaComposition(self):
        """Return dictionary of all 20 amino acids with count of amino acid characters after decoding codons"""
        aaList = [self.rnaCodonTable[k] for k in self.codonList if
                  k in self.rnaCodonTable]  # puts amino acid in list if codon from list is a key in rnaCodonTable
        newDict = {}  # create dict for amino acids with count
        aaDict = {'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0, 'M': 0, 'N': 0,
                  'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0}
        for aa in aaDict.keys():
            aaCnt = 0
            for c in aaList:
                if aa == c:
                    aaCnt += 1  # adds 1 to aaCnt if amino acid is found in list and in aaDict
            newDict[aa] = aaCnt  # store count in newDict
        aaDict.update(newDict)  # update aaDict with counts of amino acids stored in newDict
        return aaDict

    def nucComposition(self):
        """Returns dictionary of all valid nucleotide bases in input with counts"""
        validNucsDict = {'A': 0, 'G': 0, 'C': 0, 'T': 0, 'U': 0, 'N': 0}  # dictionary with only valid bases
        codonStr = ''
        codonStr = codonStr.join(self.codonList)  # turn list back to string to iterate through bases
        newNucDict = {}  # create new dict for nucleotides with count
        for aa in validNucsDict.keys():
            aaCnt = 0
            for c in codonStr:  # compare to each element in codonStr
                if aa == c:
                    aaCnt += 1
            newNucDict[aa] = aaCnt  # store count in newNucDict
        validNucsDict.update(newNucDict)
        return validNucsDict

    def codonComposition(self):
        """Returns dictionary of RNA codons and each of their counts after filtering out any with invalid bases"""
        from collections import defaultdict
        newCodons = []  # new codon list only with valid base codons
        rnaCodonList = [codon.replace('T', 'U') for codon in
                        self.codonList]  # replaces T in codon for U to use rnaCodonTable
        for codon in rnaCodonList:
            p = 0
            while p < len(codon):
                if codon[p] not in ['A', 'C', 'G', 'U']:
                    break  # removes codons with invalid bases including N
                p += 1  # to continue iteration through codon
            newCodons.append(codon)
        codonDict = defaultdict(int)  # create new dict to put valid codons with count
        for codon in newCodons:
            codonDict[codon] += 1  # for each occurrence of codon, 1 is added to total count of codon
        return codonDict

    def nucCount(self):
        """Returns sum value from dictionary returned in nucComposition"""
        return sum(self.validNucsDict)  # sum valid bases dictionary to get total number of bases

import sys
class FastAreader:
    '''
    Define objects to read FastA files.

    instantiation:
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''

    def __init__(self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname

    def doOpen(self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''

        with self.doOpen() as fileH:

            header = ''
            sequence = ''

            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>'):
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

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

class OrfFinder:
    """Analyzes input given through a FASTA-formatted data sequence file and finds ORFs
    that are at leats 100 bases long and defined by start and stop codons"""

    def __init__(self, seq):
        """Initializes input sequence as self"""
        self.seq = seq #initializes as self
        my_dna = Seq(self.seq, generic_dna) #uses Bio.Seq to take input with RNA alphabet
        self.reverseComp = my_dna.reverse_complement() #gets reverse complement from input

    def p1Frame(self):
        """Identifies ORFs within plus 1 frame and returns to plus 1 frame list"""
        inputSeq = Seq(self.seq)
        i = 0 #set all iteration values to 0
        t = 0
        start = 0 #initial start value to reset if loop entered
        n = 3 #splice factor
        p1List = []  # list to store P1 frames
        startCodons = ["ATG", "TTG", "GTG"]  # list of start codons
        stopCodons = ["TAG", "TGA", "TAA"]  # list of stop codons

        for i in range(0, len(inputSeq), n):  #check for start codon after splitting input into codons
            p1Frame = [inputSeq[i:i + n]] #create p1Frame list of all codons in sequence
            thisFrame = "+1" #name of frame
            if p1Frame[0] in startCodons: #check first codon for start and enter loop if condition met
                start = 1 #set start value to 1 to confirm loop entered
                begFrame = i #beginning frame set to position of first nuc base of start codon
                endFrame = 0 #end frame set to 0 because not found yet
                t = i + 3 #next iteration value set to beginning of next codon
                for t in range(t, len(inputSeq), 3):  #check for end if start
                    if inputSeq[t:t + 3] in stopCodons: #iterates through rest of spliced sequence
                        endFrame = t + 3 #if found end frame set to last nuc base position of end codon
                        break; #break out of for loop once end is found

                    if endFrame == 0:  #if no end
                        endFrame = len(inputSeq) #end frame set to length of input
                frameLen = endFrame - begFrame #calculates frame length
                newOrfP1 = (thisFrame + " " + str(begFrame + 1) + " " + str(endFrame) + " " + str(frameLen)) #for print statement
                p1List.append(newOrfP1) #appends frame into frame list

            if start == 0:  #no start
                begFrame = 0
                endFrame = 0
                start = 1 #reset start val to 1 to confim loop entry
                z = 0 #new iteration value set to 0 if start not found
                for z in range(z, len(inputSeq), 3):  #check for end if no start
                    if inputSeq[z:z + 3] in stopCodons: #
                        endFrame = z + 3

                    if endFrame == 0:  # if no end
                        endFrame = len(inputSeq)
                frameLen = endFrame - begFrame
                newFrameP1 = (thisFrame + " " + str(begFrame + 1) + " " + str(endFrame) + " " + str(frameLen))
                p1List.append(newFrameP1) #appends frame to frame list

        return p1List #returns frame list

    def p2Frame(self):
        """Identifies ORFs within plus 2 frame and returns to plus 2 frame list"""
        inputSeq = Seq(self.seq)
        i = 0 #iteration values set to 0
        t = 0
        start = 0 #start values set to 0
        n = 3 #splice factor
        p2List = []  # list to store P2 frames
        startCodons = ["ATG", "TTG", "GTG"]  # ist of start codons
        stopCodons = ["TAG", "TGA", "TAA"]  #list of stop codons

        for i in range(1, len(inputSeq), n):  #check for start
            p2Frame = [inputSeq[i:i + n]] #list of input spliced into codons
            thisFrame = "+2" #frame name
            if p2Frame[0] in startCodons: #iterates through codons checking for start
                start = 1 #if condition met loop entered and start value set to 1
                begFrame = i #beg frame val set to position of first nuc in start codon
                endFrame = 0 #end set to 0
                t = i + 3 #next iteration val set to beg of next codon
                for t in range(t, len(inputSeq), 3):  #check for end if start
                    if inputSeq[t:t + 3] in stopCodons: #checks from next codon for stop codon
                        endFrame = t + 3 #if found sets end val to last nuc of stop codon
                        break; #breaks out of loop once stop found

                    if endFrame == 0:  #if no end
                        endFrame = len(inputSeq) #end set to length of seq
                frameLen = endFrame - begFrame #calculates frame length
                newOrfP2 = (thisFrame + " " + str(begFrame + 1) + " " + str(endFrame) + " " + str(frameLen)) #for print statement
                p2List.append(newOrfP2) #appends frame to frame list

            if start == 0:  #if no start
                begFrame = 0 #start and end position set to 0
                endFrame = 0
                start = 1 #if condition met loop entered and start val set to 1
                z = 1 #checks from first full codon
                for z in range(z, len(inputSeq), 3):  #check for end if no start
                    if inputSeq[z:z + 3] in stopCodons: #checks next occurring codons for stop
                        endFrame = z + 3 #if found end frame set to last nuc in stop codon

                    if endFrame == 0:  #if no end
                        endFrame = len(inputSeq) #set end frame to length of input seq
                frameLen = endFrame - begFrame #calculates frame length
                newFrameP2 = (thisFrame + " " + str(begFrame + 1) + " " + str(endFrame) + " " + str(frameLen)) #for input seq
                p2List.append(newFrameP2) #appends frame to frame list

        return p2List #returns frame list


    def p3Frame(self):
        """Identifies ORFs within plus 3 frame and returns to plus 3 frame list"""
        inputSeq = Seq(self.seq)
        i = 0 #iteration values set to 0
        t = 0
        start = 0 #start set to 0
        n = 3 #splice factor
        p3List = []  # list to store P3 frames
        startCodons = ["ATG", "TTG", "GTG"]  # list of start codons
        stopCodons = ["TAG", "TGA", "TAA"]  # list of stop codons

        for i in range(2, len(inputSeq), n):  #check for start from after second base in input
            p3Frame = [inputSeq[i:i + n]] #list of codons after spliced
            thisFrame = "+3" #frame name
            if p3Frame[0] in startCodons: #iterates through codons looking for start
                start = 1 #if start found start value set to 1
                begFrame = i #beg val set to first nuc in start codon
                endFrame = 0 #end frame set to 0 because not found yet
                t = i + 3 #next iteration value set to beg of next codon
                for t in range(t, len(inputSeq), 3):  #check for end if start
                    if inputSeq[t:t + 3] in stopCodons: #cehcks rest of codons for stop
                        endFrame = t + 3 #if found end val set to last nuc of stop codon
                        break; #breaks out of loop once stop found

                    if endFrame == 0:  #if no end
                        endFrame = len(inputSeq) #if condition met sets end val to length of sequence
                frameLen = endFrame - begFrame #calculates frame length
                newOrfP3 = (thisFrame + " " + str(begFrame + 1) + " " + str(endFrame) + " " + str(frameLen)) #for print statement
                p3List.append(newOrfP3) #appends frame to frame list

            if start == 0:  #no start
                begFrame = 0 #if met sets beg and end values to 0
                endFrame = 0
                start = 1 #start value set to 1 once loop entered
                z = 2 #iteration value starts from first valid codonn
                for z in range(z, len(inputSeq), 3):  #check for end if no start
                    if inputSeq[z:z + 3] in stopCodons: #iterates through codons checking for stop codon
                        endFrame = z + 3 #if found stop set to position of last nuc base in codon

                    if endFrame == 0:  #if no end
                        endFrame = len(inputSeq) #end is length of input
                frameLen = endFrame - begFrame #calculates frame length
                newFrameP3 = (thisFrame + " " + str(begFrame + 1) + " " + str(endFrame) + " " + str(frameLen)) #for print statement
                p3List.append(newFrameP3) #appends frame to frame list

        return p3List #returns frame list


    def m1Frame(self):
        """Identifies ORFs within minus 1 frame and returns to minus 1 frame list"""
        i = 0 #iteration vals set to 0
        t = 0
        start = 0 #start value set to 0
        n = 3 #splice factor
        m1List = []  # list to store M1 frames
        startCodons = ["ATG", "TTG", "GTG"]  # ist of start codons
        stopCodons = ["TAG", "TGA", "TAA"]  #list of stop codons

        for i in range(0, len(self.reverseComp), n):  # check for start
            m1Frame = [self.reverseComp[i:i + n]] #list of codons in reverse seq after seq spliced
            thisFrame = "-1" #frame name
            if m1Frame[0] in startCodons: #iterates through all codons in reverse seq
                start = 1 #if loop entered start sets to 1
                begFrame = i #if start codon then beg frame value is first nuc of start codon
                endFrame = 0 #end value set to 0
                t = i #next iteration starts from start codon
                for t in range(t, len(self.reverseComp), 3):  #check for end if start
                    if self.reverseComp[t:t + 3] in stopCodons:
                        endFrame = len(self.reverseComp) - i #end position set to length of reverse seq minus first nuc of start
                        begFrame = len(self.reverseComp) - (t+3) #beg frame position set to length of reverse seq minus end of end codon
                        break; #breaks out of for loop

                if endFrame == 0:  #if no end
                    endFrame = len(self.reverseComp) - i #end set to length minus beg of start codon
                    begFrame = len(self.reverseComp) - (t + 3) #beg set to length minus end of end codon
                frameLen = endFrame - begFrame #calculates frame length
                newOrfM1 = (thisFrame + " " + str(begFrame + 1) + " " + str(endFrame) + " " + str(frameLen)) #for print
                m1List.append(newOrfM1) #appends frame to frame list

            if start == 0:  #if no start
                begFrame = 0 #start and end vals set to 0
                endFrame = 0
                start = 1 #if loop entered then start val sets to 1
                z = 0 #iteration starts from first nuc
                for z in range(z, len(self.reverseComp), 3):  #check for end if no start
                    if self.reverseComp[z:z + 3] in stopCodons: #checks each codon for stop
                        endFrame = len(self.reverseComp) #if found then end is length of reverse seq
                        begFrame = len(self.reverseComp) - (z + 3) #start set to length minus last nuc of stop codon

                if endFrame == 0:  # if no end
                    endFrame = len(self.reverseComp) #sets end to length of reverse seq
                frameLen = endFrame - begFrame #calculates frame length
                newFrameM1 = (thisFrame + " " + str(begFrame + 1) + " " + str(endFrame) + " " + str(frameLen)) #for print statement
                m1List.append(newFrameM1) #appends frame to frame list

        return m1List #returns frame list

    def m2Frame(self):
        i = 0 #iteration set to 0
        t = 0
        start = 0 #start set to 0 so that it can be reset if wanted codon is found
        n = 3 #splice factor
        m2List = []  # list to store M2 frames
        startCodons = ["ATG", "TTG", "GTG"]  #list of start codons
        stopCodons = ["TAG", "TGA", "TAA"]  #list of stop codons

        for i in range(1, len(self.reverseComp), n):  #check for start
            m2Frame = [self.reverseComp[i:i + n]] #frame list of spliced reverse input in codons
            thisFrame = "-2" #name of frame
            if m2Frame[0] in startCodons: #iterates through each codon looking for start
                start = 1 #if found loop entered and start value set to 1
                begFrame = i #begin frame val set to position of first nuc in start codon
                endFrame = 0 #end frame set to 0 because not found yet
                t = i + 3
                for t in range(t, len(self.reverseComp), 3):  # check for end if start
                    if self.reverseComp[t:t + 3] in stopCodons:
                        endFrame = len(self.reverseComp) - i #if stop found set end to length minus start codon position
                        begFrame = len(self.reverseComp) - (t + 3) #beg set to length minus last nuc position in stop codon
                        break; #break out of loop once stop found

                if endFrame == 0:  #if no end
                    endFrame = len(self.reverseComp) - i #end set to length of seq minus start codon position
                    begFrame = len(self.reverseComp) - (t + 3)
                frameLen = endFrame - begFrame #calculates frame length
                newOrfM2 = (thisFrame + " " + str(begFrame + 1) + " " + str(endFrame) + " " + str(frameLen)) #for print statement
                m2List.append(newOrfM2) #appends frame to frame list

            if start == 0:  #if no start
                begFrame = 0 #set beg and end frame vals to 0
                endFrame = 0
                start = 1 #if loop entered set start to 1
                z = 1 #iterates from first valid codon
                for z in range(z, len(self.reverseComp), 3):  #check for end if no start
                    if self.reverseComp[z:z + 3] in stopCodons: #checks following codons for stop
                        endFrame = len(self.reverseComp) #if found sets to length of seq
                        begFrame = len(self.reverseComp) - (z + 3) #beg frame set to length minus last nuc in codon

                if endFrame == 0:  #if no end
                    endFrame = len(self.reverseComp) #set end to length of reverse comp input
                frameLen = endFrame - begFrame #calculates frame length
                newFrameM2 = (thisFrame + " " + str(begFrame + 1) + " " + str(endFrame) + " " + str(frameLen)) #for print
                m2List.append(newFrameM2) #appends frame to frame list

        return m2List #returns frame list

    def m3Frame(self):
        i = 0 #iteration set to 0
        t = 0
        start = 0 #start set to 0 so that it can be reset if wanted codon is found
        n = 3
        m3List = []  # list to store M3 frames
        startCodons = ["ATG", "TTG", "GTG"]  # list of start codons
        stopCodons = ["TAG", "TGA", "TAA"]  # list of stop codons

        for i in range(2, len(self.reverseComp), n):  #checks in range starting from val 2 to end of sequence for frame -3
            m3Frame = [self.reverseComp[i:i + n]] #sets m3Frame to be list of codons in reverseComp strand
            thisFrame = "-3" #identifies which frame is being checked to put into print statement
            if m3Frame[0] in startCodons: #check for start
                start = 1 #start value is set to 1 if loop is entered
                begFrame = i #start position set to i
                endFrame = 0
                t = i + 3 #sets second iteration val to end of start codon
                for t in range(t, len(self.reverseComp), 3): #for rest of sequence after start is found
                    if self.reverseComp[t:t + 3] in stopCodons: #check for end if start
                        endFrame = len(self.reverseComp) - i #calculates end and start backwards
                        begFrame = len(self.reverseComp) - (t + 3)
                        break; #break out of loop once end codon is found so that if another end is found it is not taken

                if endFrame == 0:  #if no end
                    endFrame = len(self.reverseComp) - i #reverse end and beg frame start calculations as of plus frames
                    begFrame = len(self.reverseComp) - (t + 3)
                frameLen = endFrame - begFrame #calculates frame length
                newOrfM3 = (thisFrame + " " + str(begFrame + 1) + " " + str(endFrame) + " " + str(frameLen)) #for print statement
                m3List.append(newOrfM3) #appends to frame list

            if start == 0:  # no start
                begFrame = 0
                endFrame = 0
                start = 1 #resets start value to 1 once loop is entered for no start
                z = 2
                for z in range(z, len(self.reverseComp), 3):  #check for end if no start
                    if self.reverseComp[z:z + 3] in stopCodons: #splices rest of sequence into codons to check for stop codon
                        endFrame = len(self.reverseComp) #end set to length of reverse seq
                        begFrame = len(self.reverseComp) - (z + 3) #beg set to length minus last nuc in stop

                if endFrame == 0:  #if no end
                    endFrame = len(self.reverseComp) #if no end sets end position to reverse comp input length
                frameLen = endFrame - begFrame #calculates frame length
                newFrameM3 = (thisFrame + " " + str(begFrame + 1) + " " + str(endFrame) + " " + str(frameLen)) #for print statement
                m3List.append(newFrameM3) #appends to frame list

        return m3List #returns frame list with all ORFS for that frame

    def fullComp(self):
        """Appends all frame lists to one full list and returns full list"""
        fullList = [] #create a full list to put all frame lists into
        fullList.extend(self.p1Frame()) #extends full list with addition of each frame
        fullList.extend(self.p2Frame())
        fullList.extend(self.p3Frame())
        fullList.extend(self.m1Frame())
        fullList.extend(self.m2Frame())
        fullList.extend(self.m3Frame())
        return fullList #returns full list of all frame ORFs

    def sortFullComp(self):
        """Uses lambda function to sort fullList by decreasing ORF frame size and returns edited full list"""
        res = sorted(self.fullComp(), key=sort_by, reverse=True) #https://www.geeksforgeeks.org/python-sort-given-list-of-strings-by-part-of-string/
        return res #returns sorted fullComp

def sort_by(x):
    """Sorts full list by splitting string, sorting by last integer and
    in the case of same value, sorts by beginning frame position"""
    t,b,e,l = x.split() #algorithm of sorts showing how to sort
    return (int(l),-int(b)) #returns algorithm for use in sort method




