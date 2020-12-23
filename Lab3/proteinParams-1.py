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
        """Initializes input protein sequence to self after making it uppercase"""
        protein = protein.upper()  #makes input uppercase and initializes to self
        self.protein = protein.replace(" ", "") #gets rid of whitespaces

    def aaCount(self):
        """Returns single int count of aa characters found in input string if amino acid from input are found in given aa2mw dictionary"""
        count = 0 #sets count to 0
        for c in self.protein: #for every c found in sequence
            for aa in self.aa2mw.keys(): #for keys of aa2mw dict
                if aa == c: #if amino acids from input and keys match
                    count +=1 #add 1 to count
        return count #returns count

    def aaComposition(self):
        """Computes and returns aaCompDict by searching for all occurrences of an amino acid and adding 1 to the count if repeats
        of an amino acid are found"""
        aaCompDict = {'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'L': 0, 'K': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0,
         'R': 0, 'S': 0, 'T': 0, 'V': 0, 'Y': 0, 'W': 0} #creates aa dict where all counts of an amino acid are initially set to 0
        editDict = {} #creates empty dict to use for updating aaComp
        for aa in self.protein: #for each aa in input
            if aa in editDict: #if the aa exists in editDict
                editDict[aa] += 1 #add 1 to count
            else: #if not
                editDict[aa] = 1 #set count to 1
        aaCompDict.update(editDict) #update counts of aaComp with editDict
        return aaCompDict #returns updated aaComp dict

    def _charge_(self, bestPh):
        """Calculates individual charge of charge-holding amino acids in input string by positive and negative groups and
        returned best pH value from pI methodand returns net charge """
        chargeDict = self.aaComposition()  #uses aaComp dict
        kCharge = chargeDict['K'] * ((10 ** ProteinParam.aa2chargePos['K']) / ((10 ** ProteinParam.aa2chargePos['K']) + 10 ** (bestPh)))  #bestpH value returned by pI
        rCharge = chargeDict['R'] * ((10 ** ProteinParam.aa2chargePos['R']) / ((10 ** ProteinParam.aa2chargePos['R']) + 10 ** (bestPh))) #calculation for positive charge
        hCharge = chargeDict['H'] * ((10 ** ProteinParam.aa2chargePos['H']) / ((10 ** ProteinParam.aa2chargePos['H']) + 10 ** (bestPh))) #** means to the power of
        nTermCharge = ((10 ** ProteinParam.aaNterm) / ((10 ** ProteinParam.aaNterm) + 10 ** (bestPh))) #calculates positive N term charge
        posTotal = kCharge + rCharge + hCharge + nTermCharge  #sum of all positive charges

        dCharge = chargeDict['D'] * ((10 ** (bestPh)) / ((10 ** ProteinParam.aa2chargeNeg['D']) + 10 ** (bestPh))) #calculation for negative charge
        eCharge = chargeDict['E'] * ((10 ** (bestPh)) / ((10 ** ProteinParam.aa2chargeNeg['E']) + 10 ** (bestPh)))
        cCharge = chargeDict['C'] * ((10 ** (bestPh)) / ((10 ** ProteinParam.aa2chargeNeg['C']) + 10 ** (bestPh)))
        yCharge = chargeDict['Y'] * ((10 ** (bestPh)) / ((10 ** ProteinParam.aa2chargeNeg['Y']) + 10 ** (bestPh)))
        cTermCharge = ((10 ** bestPh) / ((10 ** ProteinParam.aaCterm) + 10 ** (bestPh))) #calculates negative C term charge
        negTotal = dCharge + eCharge + cCharge + yCharge + cTermCharge  #sum of all negative charges

        total = posTotal - negTotal #subtract negCharge total from posCharge total
        return total #returns total charge of input

    def pI(self):
        """Finds pH producing as close to a neutral net charge as possible when put into charge method and returns that value"""
        bestCharge = 100000 #best charge is initially set to some arbitrarily large value to ensure condition wont be met
        bestPh = 0 #bestpH initially set to 0
        for value in range(0, 1400 + 1): #for value in the standard pH range
            value = value / 100 #to get standard pH values (1,2,etc)
            charge = abs(self._charge_(value)) #take absolute value of charge so negative values are not obtained
            if charge < bestCharge: #as long as abs charge is less than bestCharge; definitely will be bc bestCharge is very high
                bestCharge = charge #set abs charge value to be new bestCharge value
                bestPh = value #new bestpH value is the one that results in the neutral net charge
        return bestPh #returns bestpH value

    def molarExtinction(self):
        """Calculates protein light absorption at wavelength of 280nm and returns total value"""
        yAbs = ((self.aa2abs280['Y']) * self.aaComposition()['Y'])
        wAbs = ((self.aa2abs280['W']) * self.aaComposition()['W'])
        bAbs = ((self.aa2abs280['C']) * self.aaComposition()['C'])
        total = yAbs + wAbs + bAbs #sums all absorption values and sets to total
        return total #returns total absorption

    def massExtinction(self):
        """Divides molar extinction coefficient by molecular weight and returns mass extinction coefficient value"""
        myMW = self.molecularWeight() #obtains mol weight from molecularWeight method
        return self.molarExtinction() / myMW if myMW else 0.0 #divides  and returns molar extinction coeff by mol weight if possible
        #and if not returns mass extinction coeff as 0

    def molecularWeight(self):
        """Returns molecular weight of protein sequence by obtaining the weight of each amino acid
         from the dictionary, subtracting the weight of water and then multiplying by the aa count which is then
         added back to the weight of water to return the final molecular weight"""
        weightDict = self.aaComposition() #uses aaComp dict
        molWeight = 0  #set initial weight to 0
        for key in weightDict:  #the weight of every aa in the dict is calculated and added to final weight
            weight = ProteinParam.aa2mw.get(key)  #obtains aa weight from dict
            mw1 = weight - ProteinParam.mwH2O  #the weight of water subtracted from weight of aa
            mw2 = weightDict[key] * mw1  #last value is multiplied with count of  amino aicd
            molWeight += mw2  #this is added to sum of the molecular weights
        molWeight = molWeight + ProteinParam.mwH2O  #final sum of molecular weights is added to weight of water
        return molWeight #returns final molecular weight

import sys
def main():
    """Prompts the user to enter a protein sequence which upon receiving will be sent to the above methods to
     calculate different physical-chemical properties of the protein sequence and then print a formatted version of these values"""
    inString = input('protein sequence?') #prompts user for input seq
    while inString:
        myParamMaker = ProteinParam(inString) #sends to methods
        myAAnumber = myParamMaker.aaCount() #sends to count method to get count of aa in input seq
        print("Number of Amino Acids: {aaNum}".format(aaNum=myAAnumber)) #puts aaCount value into print statement
        print("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight())) #rounds molecular weight to 1 decimals
        print("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction())) #rounds molar coeff to 2 decimals
        print("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction())) #rounds mass coeff to 2 decimals
        print("Theoretical pI: {:.2f}".format(myParamMaker.pI())) #rounds pI to 2 decimals
        print("Amino acid composition:")
        myAAcomposition = myParamMaker.aaComposition()
        keys = list(myAAcomposition.keys()) #lists all aa in aaComp dict
        keys.sort() #to sort keys
        if myAAnumber == 0: myAAnumber = 1  #if no AAs are present, for loop executed
        for key in keys:
            print("\t{} = {:.2%}".format(key, myAAcomposition[key] / myAAnumber)) #formats by putting aa followed by its percent makeup of input seq

        inString = input('protein sequence?') #prompts new input

if __name__ == "__main__":
    main() #calls main