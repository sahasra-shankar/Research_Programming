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
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
    }

    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__ (self, protein): 
        """initialize self from input protein sequence and compute aaComposition dictionary"""
        self.protein = protein.upper() #make input uppercase
        self.aaCompositionDict = dict()
        for aa in self.aa2mw.keys():
            aaCnt=0
            for seq in list(self.protein):
                if aa == seq:
                    aaCnt+=1
            self.aaCompositionDict[aa]=aaCnt

    def aaCount (self): 
        """return single int count of aa characters found in input string"""
        for aa in self.aa2mw.keys():
            result = 0 
            for seq in list(self.protein):
                if aa == seq:
                    result += 1
                elif aa != seq:
                    result += 0
        return result

    def pI (self): 
        """find pH producing neutral net charge when put into charge method and return"""
        bestCharge = 100000
        bestPh = 0
        for value in range(0,1400+1):
            value = value/100
            charge = abs(self._charge_(value))
            if charge < bestCharge:
                bestCharge = charge
                bestPh = value
        return bestPh

    def aaComposition (self): 
        """return aaComposition dictionary computed in init"""
        return self.aaCompositionDict

    def _charge_ (self, bestPh): 
        """calculates individual charge of amino acids having charge in input string and returns net charge """
        for aa in self.aa2mw.keys():
            for seq in list(self.protein):
                if aa == seq:
                    rCharge = self.aaCompositionDict['R']*(10**self.aa2chargePos['R']/((10**self.aa2chargePos['R'])+(10**bestPh)))
                    kCharge = self.aaCompositionDict['K']*(10**self.aa2chargePos['K']/((10**self.aa2chargePos['K'])+(10**bestPh)))
                    hCharge = self.aaCompositionDict['H']*(10**self.aa2chargePos['H']/((10**self.aa2chargePos['H'])+(10**bestPh)))
                    dCharge = self.aaCompositionDict['D']*(10**bestPh/((10**self.aa2chargeNeg['D'])+(10**bestPh)))
                    eCharge = self.aaCompositionDict['E']*(10**bestPh/((10**self.aa2chargeNeg['E'])+(10**bestPh)))
                    cCharge = self.aaCompositionDict['C']*(10**bestPh/((10**self.aa2chargeNeg['C'])+(10**bestPh)))
                    yCharge = self.aaCompositionDict['Y']*(10**bestPh/((10**self.aa2chargeNeg['Y'])+(10**bestPh)))
        
                    pChargeTotal = rCharge+kCharge+hCharge+((10**self.aaNterm)/((10**self.aaNterm)+(10**bestPh)))
                    nChargeTotal = dCharge+eCharge+cCharge+yCharge+((10**bestPh)/((10**self.aaCterm)+(10**bestPh)))
                    netCharge = pChargeTotal - nChargeTotal
                elif aa != seq: #assume if non-amino acid exists then net charge value does not exist 
                    netCharge = 0
        return netCharge
    
    def molarExtinction (self): 
        """calculates protein light absorption at wavelength of 280nm and returns total value"""
        yAbs = ((self.aa2abs280['Y'])*self.aaCompositionDict['Y'])
        wAbs = ((self.aa2abs280['W'])*self.aaCompositionDict['W'])
        bAbs = ((self.aa2abs280['C'])*self.aaCompositionDict['C'])
        total = yAbs+wAbs+bAbs
        return total

    def massExtinction (self): 
        """calculates by dividing molar extinction coefficient by molecular weight and returns value"""
        myMW =  self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0

    def molecularWeight (self): 
        """returns molecular weight of protein sequence"""
        newList = []
        for aa in self.aa2mw.keys():
            for seq in list(self.protein):
                if aa == seq:
                    myAA = self.aa2mw[aa]
                    newList.append(myAA)
                    sumAA = sum(newList) #sums all weights in new list
                    #https://appdividend.com/2019/08/09/python-sum-example-sum-function-in-python-tutorial/
                    molWeight = self.mwH2O + (sumAA-(self.aaCount()*self.mwH2O)) 
                elif aa != seq: #assume if non-amino acid exists then molecular weight cannot be calculated
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
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}

    def __init__ (self, inSeq, inString=''):
        """Splits input by every three characters and puts into new list while accounting for any additional input that is added"""
        self.inString = inString
        self.inSeq = inSeq
        finalInput = inString + inSeq #sets input to initial input plus anything that was added to it after
        n = 3  #https://stackoverflow.com/questions/9475241/split-string-every-nth-character
        codonList = [finalInput[i:i+n] for i in range(0, len(finalInput), n)] #list of codons
        self.codonList = codonList #initialized codonList to self
        
    def addSequence (self):
        """Returns codon list computed in init object because init cannot return lists"""
        return self.codonList
    
    def aaComposition(self):
        """Return dictionary of all 20 amino acids with count of amino acid characters after decoding codons"""
        aaList =[self.rnaCodonTable[k] for k in self.codonList if k in self.rnaCodonTable] #puts amino acid in list if codon from list is a key in rnaCodonTable
        newDict = {}  #create dict for amino acids with count 
        aaDict = {'A':0,'C':0,'D':0,'E':0,'F':0,'G':0,'H':0,'I':0,'K':0,'L':0,'M':0,'N':0,'P':0,'Q':0,'R':0,'S':0,'T':0,'V':0,'W':0,'Y':0}
        for aa in aaDict.keys():
            aaCnt = 0
            for c in aaList:
                if aa == c:
                    aaCnt += 1  #adds 1 to aaCnt if amino acid is found in list and in aaDict
            newDict[aa] = aaCnt #store count in newDict
        aaDict.update(newDict) #update aaDict with counts of amino acids stored in newDict
        return aaDict
    
    def nucComposition(self):
        """Returns dictionary of all valid nucleotide bases in input with counts"""
        validNucsDict = {'A':0, 'G':0, 'C':0, 'T':0, 'U':0, 'N':0} #dictionary with only valid bases
        codonStr = ''
        codonStr = codonStr.join(self.codonList) #turn list back to string to iterate through bases
        newNucDict = {}  #create new dict for nucleotides with count 
        for aa in validNucsDict.keys():
            aaCnt = 0
            for c in codonStr: #compare to each element in codonStr
                if aa == c:
                    aaCnt += 1
            newNucDict[aa] = aaCnt #store count in newNucDict
        validNucsDict.update(newNucDict)
        return validNucsDict
    
    def codonComposition(self):
        """Returns dictionary of RNA codons and each of their counts after filtering out any with invalid bases"""
        from collections import defaultdict
        newCodons = []  #new codon list only with valid base codons
        rnaCodonList = [codon.replace('T', 'U') for codon in self.codonList] #replaces T in codon for U to use rnaCodonTable
        for codon in rnaCodonList:
            p = 0
            while p < len(codon):
                if codon[p] not in ['A','C','G','U']:
                    break  #removes codons with invalid bases including N
                p+=1 #to continue iteration through codon
            newCodons.append(codon)
        codonDict = defaultdict(int)   #create new dict to put valid codons with count  
        for codon in newCodons:
            codonDict[codon] += 1  #for each occurrence of codon, 1 is added to total count of codon
        return codonDict
    
    def nucCount(self):
        """Returns sum value from dictionary returned in nucComposition"""
        return sum(self.validNucsDict) #sum valid bases dictionary to get total number of bases

import sys    
class FastAreader :
    ''' 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence
        
from Bio.Seq import Seq
class OrfFinder:
    """Analyzes input given through a FASTA-formatted data sequence file and finds ORFs
    that are at leats 100 bases long and defined by start and stop codons"""

    def __init__(self, seq):
        """Initializes input sequence as self"""
        self.seq = seq
        transDict = {"A":"T","T":"A","C":"G", "G":"C"}
        reverseComp = seq[::-1].maketrans(transDict)  # reverse input
        self.reverseComp = reverseComp
        
    def p1Frame(self):
        inputSeq = Seq(self.seq)
        i = 0
        t = 0
        start = 0
        n = 3
        p1List = [] # list to store P1 frames
        startCodons = ["ATG", "TTG", "GTG"]  # list of start codons
        stopCodons = ["TAG", "TGA", "TAA"]  # list of stop codons

        for i in range(0, len(inputSeq), n): #check for start
            p1Frame = [inputSeq[i:i+n]]
            thisFrame = "+1"
            if p1Frame[0] in startCodons:
                start = 1
                begFrame = i
                endFrame = 0
                t = i+3
                for t in range(t, len(inputSeq), 3):  # check for end if start
                    if inputSeq[t:t+3] in stopCodons:
                        endFrame = t+3
                        break;
                        
                    if endFrame == 0: #if no end if start
                        endFrame = len(inputSeq)
                frameLen = endFrame - begFrame
                newOrfP1 = (thisFrame + " " + str(begFrame+1) + ".." + str(endFrame) + " " + str(frameLen))
                p1List.append(newOrfP1)
                        
            if start == 0:  # no start
                begFrame = 0
                endFrame = 0
                z = 0
                for z in range(z, len(inputSeq), 3):  # check for end if no start
                    if inputSeq[z:z+3] in stopCodons:
                        endFrame = z+3

                    if endFrame == 0:  # if no end
                        endFrame = len(inputSeq)
                frameLen = endFrame - begFrame
                newFrameP1 = (thisFrame + " " + str(begFrame+1) + ".." + str(endFrame) + " " + str(frameLen))
                p1List.append(newFrameP1)
            
        return p1List
    
    def p2Frame(self):
        inputSeq = Seq(self.seq)
        i = 0
        t = 0
        start = 0
        n = 3
        p2List = [] # list to store P2 frames
        startCodons = ["ATG", "TTG", "GTG"]  # list of start codons
        stopCodons = ["TAG", "TGA", "TAA"]  # list of stop codons

        for i in range(1, len(inputSeq), n): #check for start
            p2Frame = [inputSeq[i:i+n]]
            thisFrame = "+2"
            if p2Frame in startCodons:
                start = 1
                begFrame = i
                endFrame = 0
                t = i
                for t in range(t, len(inputSeq), 3):  # check for end if start
                    if inputSeq[t:t+3] in stopCodons:
                        endFrame = t+3
                        break;
                        
                if endFrame == 0:
                    endFrame = len(inputSeq)
                    frameLen = endFrame - begFrame
                    newOrfP2 = (thisFrame + " " + str(begFrame+1) + ".." + str(endFrame) + " " + str(frameLen))
                    p2List.append(newOrfP2)
                
            if start == 0:  # no start
                begFrame = 0
                endFrame = 0
                z = 0
                for z in range(z, len(inputSeq), 3):  # check for end if no start
                    if inputSeq[z:z+3] in stopCodons:
                        endFrame = z+3
                        frameLen = endFrame - begFrame

                    if endFrame == 0:  # if no end
                        endFrame = len(inputSeq)
                        frameLen = endFrame - begFrame
                        newFrameP2 = (thisFrame + " " + str(begFrame+1) + ".." + str(endFrame) + " " + str(frameLen))
                        p2List.append(newFrameP2)
            return p2List
        
    def p3Frame(self):
        pass
    def m1Frame(self):
        
        pass
    def m2Frame(self):
        pass
    def m3Frame(self):
        pass
    
    def fullComp(self):
        fullList = []
        fullList.append(self.p1Frame())
        #fullList.append(self.p2Frame())
        #fullList.append(p3List)
        #fullList.append(m1List)
        #fullList.append(m2List)
        #fullList.append(m3List)
        return fullList
                                
                

