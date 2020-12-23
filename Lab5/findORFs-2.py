#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/env python3
# Name: Sahasra Shankar (sshanka8)
# Group Members: Shweta Jones and Pranav Muthuraman (shsujone and pmuthura)


########################################################################
# CommandLine
########################################################################
import sequenceAnalysis
import sys
class CommandLine():
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond. 
    it implements a standard command line argument parser with various argument options,
    a standard usage and help.

    attributes:
    all arguments received from the commandline using .add_argument will be
    available within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
    '''
    #Pseudocode
    #For OrfFinder class
    #initialize input as self and set as first input strand
    #use Bio.Seq package to get input sequence with RNA alphabet
    #.reverseComplement gets the reverse compliment of the input
    #make empty list for total ORFs from all frames
    #make methods for each plus and minus frame
    #in each method, look for ORFs within its designated frame
    #find start and end codons to find ORF
    #set initial start value and iteration values to 0
        #for loop for range from 0 to length of sequence splices input sequence by 3 to get codons
            #nested if loop searches for start codons in list of codons
            #if start codon is found then start value is set to 1, i is set to start position of ORF, and second iteration value set to i+3
                #once start value is found enter second for loop for second iteration value that searches for stop codons
                #if stop codon is found then mark end position of frame and calculate frame length
                #append this frame into plus frame list
            #if stop codon is not found, mark end position as input sequence length and calculate this frame length instead
        #if start codon is not found then start value remains at 0 and enters if loop for the condition that start==0
            #once this loop is entered, start value gets set to 1 so that the same frames do not get appended to frame list
            #after if loop is entered, stop codon is searched for and if found, end frame position is set to value of last base of stop codon
            #if stop codon is not found in then end frame position is set to length of input sequence and frame length is calculated accordingly
            #frame gets appended to plus frame list
    #method then returns frame list
    #repeat for all plus frames with modified ranges based on how each frame needs to be spliced
    #for minus frames use reverse compliment as input
    #same logic as plus frames except append to minus frame lists
    #make method fullComp to create list to put frame lists in and return fullList
    #make method to sort list by decreasing ORF size
        #in the case of same frame size then sort by first occurring frame position
    #call sorted list from main
    #edit so that stdin and stdout are used
        #leave file name arg blank so that code reads from stdin
        #stdout will make use of print statements

    def __init__(self, inOpts=None) :
        '''
        Implement a parser to interpret the command line argv string using argparse.
        '''
        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog - a brief description of what this thing does', 
                                             epilog = 'Program epilog - some other stuff you feel compelled to say', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s [options] -option1[default] <input >output'
                                             )
        self.parser.add_argument('inFile', action = 'store', help='input file name') #comment out for submission
        self.parser.add_argument('outFile', action = 'store', help='output file name') #comment out
        self.parser.add_argument('-lG', '--longestGene', action = 'store', nargs='?', const=True, default=False, help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, choices= (100,200,300,500,1000), default=100, action = 'store', help='minimum Gene length')
        self.parser.add_argument('-s', '--start', action = 'append', default = ['ATG'],nargs='?', help='start Codon') #allows multiple list options
        self.parser.add_argument('-t', '--stop', action = 'append', default = ['TAG','TGA','TAA'],nargs='?', help='stop Codon') #allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')  
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)

########################################################################
# Main
# Here is the main program
# 
#
########################################################################
   

def main(inCL=None):
    '''
    Calls method from sequenceAnalysis that returns sorted ORFs for six frames from the input file and
    prints with formatting
    '''
    if inCL is None:
        myCommandLine = CommandLine() #if command line not given
    else :
        myCommandLine = CommandLine(inCL) #if command line given

    print (myCommandLine.args)

    myReader = sequenceAnalysis.FastAreader(myCommandLine.args.inFile) #sends test file to FastAreader and then to OrfFinder
    if myCommandLine.args.outFile:
        f = open(myCommandLine.args.outFile, "a") #opens and reads input file

    for head, seq in myReader.readFasta():
        if myCommandLine.args.outFile:
            f.write('{:s}\n'.format(head)) #prints file information and then new line character
        else:
            print(head)

        finalNuc = sequenceAnalysis.OrfFinder(seq) #sets finalNuc to use OrfFinder from sequenceAnalysis
        geneOut = finalNuc.sortFullComp() #sets geneOut to run sortFullComp for given input


        for x in geneOut:
            geneX = x.split() #splits output from sequenceAnalysis for easy sorting
            if int(geneX[3]) <= int(myCommandLine.args.minGene): #only print higher than minGene
                continue
            if myCommandLine.args.outFile:
                f.write('{:s} {:s}..{:s}  {:s}\n'.format(geneX[0], geneX[1], geneX[2], geneX[3])) #formats geneOut with all parameters given and new line character
            else:
                print('{:s} {:s}..{:s}  {:s}'.format(geneX[0],geneX[1],geneX[2],geneX[3])) #in case of no outFile given prints with this format
        #break

    
###### replace the code between comments.
    # print (myCommandLine.args)
    # print(myCommandLine.args.minGene)
        # myCommandLine.args.inFile has the input file name
        # myCommandLine.args.outFile has the output file name
        # myCommandLine.args.longestGene is True if only the longest Gene is desired
        # myCommandLine.args.start is a list of start codons
        # myCommandLine.args.minGene is the minimum Gene length to include
        #['tass2.fa', 'tass2ORFdata-ATG-100.txt', '-lG']
        #tass2.fa tass2ORFdata-ATG-100.txt -mG=300 -lG
        #lab5test.fa
#######


if __name__ == "__main__":
    main()  # delete the list when you want to run with STDIN


# In[ ]:




