#!/usr/bin/env python3 
# Name: Sahasra Shankar (sshanka8) 
# Group Members: Shweta Jones, Pranav Muthuraman (shsujone and pmuthura)
'''
Calculate bond lengths and angles in degrees using dot product for inputted coordinates
'''

import math
class Triad : 
    """
    Calculate angles and distances among a triad of points.
 
    Author: David Bernick
    Date: March 21, 2013
    Points can be supplied in any dimensional space as long as they are consistent.
    Points are supplied as tupels in n-dimensions, and there should be three
    of those to make the triad. Each point is positionally named as p,q,r
    and the corresponding angles are then angleP, angleQ and angleR.
    Distances are given by dPQ(), dPR() and dQR()
 
    Required Modules: math
    initialized: 3 positional tuples representing Points in n-space
             p1 = Triad( p=(1,0,0), q=(0,0,0), r=(0,1,0) )
    attributes: p,q,r the 3 tuples representing points in N-space
    methods:  angleP(), angleR(), angleQ() angles measured in radians
          dPQ(), dPR(), dQR() distances in the same units of p,q,r
 
    """
 
    def __init__(self,p,q,r) :
        """ Initialize self from input of three coordinate sets to use in objects 
        
        Example construction:
            p1 = Triad( p=(1.,0.,0.), q=(0.,0.,0.), r=(0.,0.,0.) ). 
        """
        self.p = p
        self.q = q
        self.r = r
# private helper methods
    def d2 (self,a,b) : # calculate squared distance of point a to b
        return float(sum((ia-ib)*(ia-ib)  for  ia,ib in zip (a,b)))
    
    def dot (self,a,b) : # dotProd of standard vectors a,b
        return float(sum(ia*ib for ia,ib in zip(a,b)))
    
    def ndot (self,a,b,c) : # dotProd of vec. a,c standardized to b
        return float(sum((ia-ib)*(ic-ib) for ia,ib,ic in zip (a,b,c)))
    
# calculate lengths(distances) of segments PQ, PR and QR
    def dPQ (self):
        """ Provides the distance between point p and point q """
        return math.sqrt(self.d2(self.p,self.q))
    
    def dPR (self):
        """ Provides the distance between point p and point r """
        return math.sqrt(self.d2(self.p,self.r))
    
    def dQR (self):
        """ Provides the distance between point q and point r """
        return math.sqrt(self.d2(self.q,self.r))
    
    def angleP (self) :
        """ Provides the angle made at point p by segments pq and pr (radians). """
        return math.acos(self.ndot(self.q,self.p,self.r) /   math.sqrt(self.d2(self.q,self.p)*self.d2(self.r,self.p)))
    
    def angleQ (self) :
        """ Provides the angle made at point q by segments qp and qr (radians). """
        return math.acos(self.ndot(self.p,self.q,self.r) /  math.sqrt(self.d2(self.p,self.q)*self.d2(self.r,self.q)))
 
    def angleR (self) :
        """ Provides the angle made at point r by segments rp and rq (radians). """
        return math.acos(self.ndot(self.p,self.r,self.q) /  math.sqrt(self.d2(self.p,self.r)*self.d2(self.q,self.r)))

def main():
    ''' Split input into points'''
    coords = input("Enter three coordinate sets: ") #https://kite.com/python/answers/how-to-split-a-string-into-a-list-of-integers-in-python
    p1 = coords.find("(")
    p2 = coords.find(")")
    p = (coords[p1+1:p2])
    p = p.split(',')
    pMap = map(float, p) #use float to accept larger decimal values
    pList = list(pMap)

    q1 = coords.find("(",p2+1)
    q2 = coords.find(")",p2+1)
    q = (coords[q1+1:q2])
    q = q.split(',')
    qMap = map(float, q)
    qList = list(qMap)
   
    r1 = coords.find("(",q2+1)
    r2 = coords.find(")",q2+1)
    r = (coords[r1+1:r2])
    r = r.split(',')
    rMap = map(float, r)
    rList = list(rMap)
    
    
    newCoords = Triad(pList,qList,rList)
    pq = newCoords.dPQ()
    qr = newCoords.dQR()
    pqr = newCoords.angleQ()
    pqr = (pqr*180)/math.pi #conversion of radians to degrees
    
    print ("N-C bond length = " + ("%0.2f" % pq))
    print ("N-Ca bond length = " + ("%0.2f" % qr))
    print ("C-N-Ca bond length = " + ("%0.1f" % pqr))
    
    
main()