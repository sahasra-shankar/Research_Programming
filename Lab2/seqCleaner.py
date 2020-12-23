#!/usr/bin/env python3 
# Name: Sahasra Shankar (sshanka8) 
# Group Members: Shweta Jones and Pranav Muthuraman (shsujone and pmuthura)

'''
Read a DNA string from user input and return a collapsed substring of embedded Ns to: {count}.

Example: 
input: AaNNNNNNGTC
output: AA{6}GTC

Any lower case letters are converted to uppercase
'''
class DNAstring (str): 
    """Return length and get count of Ns in input stored in {}"""
    def length (self):       
        return(length(self))

    def purify(self):        
        ''' Return an upcased version of the string, collapsing a single run of Ns.'''        
        upperDNA = self.upper() #made the DNA all uppercase       
        nFirst = upperDNA.find('N')      
        nCount = upperDNA.count('N')     
        print (upperDNA[:nFirst] +"{"+ str (nCount) +"}"+ upperDNA[(nFirst+nCount):])  #returns count of N in {}       
    
def main():   
    ''' Get user DNA data and clean it up.'''   
    data = input('DNA data?')   
    thisDNA = DNAstring (data)   
    pureData = thisDNA.purify()   
    
main()
