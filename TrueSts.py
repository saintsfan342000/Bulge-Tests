
def TrueSts(P,e1,e2,R1,R2,E,v,to):
    
    from numpy import exp, sqrt,min
    
    #Include differing radii for different directions
    tnew=10000 # Initialize to erroneously high value

    #Initial approximation of thickness
    ta=to*exp(-(e1+e2))

    k=1
    while min([ta/tnew,tnew/ta])<.999:
        if k!=1:
            ta=tnew  #Update ta to the previously calculated tnew if we're not on the first iteration
        tru1=P*R1/(2*ta) #Approximate true stress based on ta
        tru2=P*R2/(2*ta)
        e1p=e1-(1/E)*(tru1-v*tru2) #e_plastic strain is e1 minus e1_elastic
        e2p=e2-(1/E)*(tru2-v*tru1)
        e3p=-(e1p+e2p) #Assume plastic incompressibility
        e3=e3p-(v/E)*(tru1+tru2)   #Then add on the elastic part of e3 (plane stress) to get e3_total
        #epeq=sqrt((2/3)*(e1p**2+e2p**2+e3p**2))   #Calcu epeq 
        tnew=to*exp(e3)  #Get new thickness
        #trunew1=P*R1/(2*tnew) #New true stress
        #trunew2=P*R2/(2*tnew) #New true stress
        k+=1
        
        #return trunew1, trunew2,e1p,e2p,e3p,epeq,k,tnew
    
    return tnew
