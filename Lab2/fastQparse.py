#!/usr/bin/env python3 
# Name: Sahasra Shankar (sshanka8) 
# Group Members: Shweta Jones and Pranav Muthuraman (shsujone and pmuthura)

'''
Parse sequence name information from single line of a FASTQ formatted file
'''

class FastqString (str):
    ''' Put input from main through object to parse information'''
           
    def parse(self):
        ''' Index information and assign to variables '''
        inst = self.find(":")
        print ("Instrument = "+self[1:inst])
        run = self.find(":",inst+1)
        print ("Run ID = "+self[inst+1:run])
        flowId = self.find(":",run+1)
        print ("Flow Cell ID = "+self[run+1:flowId])
        flowLane = self.find(":",flowId+1)
        print ("Flow Cell Lane = "+self[flowId+1:flowLane])
        tileNo = self.find(":",flowLane+1)
        print ("Tile Number = "+self[flowLane+1:tileNo])
        xCoord = self.find(":",tileNo+1)
        print ("X-Coord = "+self[tileNo+1:xCoord])
        yCoord = self.find(":",xCoord+1)
        print ("Y-Coord = "+self[xCoord+1:yCoord])
    
def main():
    ''' Get input and run through class object for indexing'''
    #assume inputs are in the same format as given example
    data = input('FASTQ sequence name line: ')
    thisData = FastqString(data)
    finalData = thisData.parse()
    
main()