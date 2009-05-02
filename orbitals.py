import math

class Orbital(object):
    
    def __init__(self,Z=1):
        
        self.Z = Z
        self.coef = 1
        
    def psi_sphr_r(self,r):
        pass
    
    def psi_sphr_a(self,p,t):
        pass
        
    def psi_sphr(self,r,p,t):
        return self.coef*self.psi_sphr_r(r)*self.psi_sphr_a(p,t)
    def psi_rect(self,x,y,z):
        ro2 = x*x+y*y
        r = (x*x + y*y + z*z)**0.5
        p = math.atan2(y,x)
        t = math.atan2(ro2**0.5,z)

        return self.psi_sphr(r,p,t)
    
    
    
class Orb_1s(Orbital):

    def __init__(self,*args,**kwargs):
        Orbital.__init__(self,*args,**kwargs)
        self.coef = (self.Z/math.pi)**0.5*self.Z**3
    
    def psi_sphr_r(self,r):
        return math.exp(-self.Z*r)

    def psi_sphr_a(self,p=0,t=0):
        return 1.

class Orb_2s(Orbital):
    
    def __init__(self,*args,**kwargs):
        Orbital.__init__(self,*args,**kwargs)
        self.coef = 0.25*(self.Z/(2*math.pi))**0.5*self.Z**3
         
    def psi_sphr_r(self,r):

        return (2.-self.Z*r)*math.exp(-self.Z*r/2.)

    def psi_sphr_a(self,p=0,t=0):
        return 1.
    
class Orb_2p(Orbital):
    
    def __init__(self,label,*args,**kwargs):

        Orbital.__init__(self,*args,**kwargs)
        self.coef = 0.25*(self.Z/(2*math.pi))**0.5*self.Z**5
        self.label = label.lower()
        
    def psi_sphr_r(self,r):
        return r*math.exp(-self.Z*r/2)
    def psi_sphr_a(self,p,t):
        if self.label == "x":
            return math.sin(t)*math.cos(p)
        elif self.label == "y":
            return math.sin(t)*math.sin(p)
        else:
            return math.cos(t)

class Orb_3s(Orbital):
    
    def __init__(self,*args,**kwargs):
        Orbital.__init__(self,*args,**kwargs)
        self.coef = 1/81.*(self.Z/(3*math.pi))**0.5*self.Z**3
         
    def psi_sphr_r(self,r):
        return (27.-18*self.Z*r+2*self.Z**2*r**2)*math.exp(-self.Z*r/3.)

    def psi_sphr_a(self,p=0,t=0):
        return 1.
        

class Orb_3p(Orbital):
    
    def __init__(self,label,*args,**kwargs):

        Orbital.__init__(self,*args,**kwargs)

        self.coef = 1/81.*(2*self.Z/math.pi)**0.5*self.Z**5
        if label not in ["x","y","z"]:
            raise TypeError("%s not in allowed types" % (label))
        else:
            self.label = label.lower()
        
    def psi_sphr_r(self,r):
        return r*(6-self.Z*r)*math.exp(-self.Z*r/3)
    
    def psi_sphr_a(self,p,t):
        if self.label == "x":
            return math.sin(t)*math.cos(p)
        elif self.label == "y":
            return math.sin(t)*math.sin(p)
        else:
            return math.cos(t)


class Orb_3d(Orbital):
    
    def __init__(self,label,*args,**kwargs):

        Orbital.__init__(self,*args,**kwargs)
        self.coef = 1/81.*(self.Z/(6*math.pi))**0.5*self.Z**7
        if label not in ["z2","x2y2","xy","xz","yz"]:
            raise TypeError("%s not in allowed types" % (label))
        else:
            self.label = label.lower()
        
    def psi_sphr_r(self,r):
        return r**2*math.exp(-self.Z*r/3.)
    def psi_sphr_a(self,p,t):
        if self.label == "z2":
            return 3*math.cos(t)**2-1
        elif self.label == "x2y2":
            return math.sin(t)**2*math.cos(2*p)
        elif self.label == "xy":
            return math.sin(t)**2*math.cos(2*p)
        elif self.label == "xz":
            return math.sin(2*t)*math.cos(p)
        elif self.label == "yz":
            return math.sin(2*t)*math.sin(p)

        


import pylab
import numpy
import parsexyz

atoms = parsexyz.ParseXYZ('./data/geom/xyz/benzene.xyz').parse()


#orb2 = Orb_2p("x",Z=8)

s = 3 
x = numpy.arange(-s,s,s/40.)
y = numpy.arange(-s,s,s/40.)
z = numpy.zeros((len(x),len(y)),float)
print z
atomic_nums = {'C':6,'H':1}

for a in atoms:
    sym,(xo,yo,zo) = a
    if sym == 'H':
        orb1 = Orb_1s(Z=atomic_nums[sym])
    else:
        orb1 = Orb_2p('Z',Z=atomic_nums[sym])

    zd = numpy.array([ orb1.psi_rect(a-xo,b-yo,zo)  for b in y for a in x]).reshape((len(x),len(y)))
    print sym,xo,yo,zo,zd.min(),zd.max()

    z+= zd

V = numpy.arange(z.min(),z.max(),(z.max()-z.min())/100.)
c = pylab.contourf(x,y,z,V)

pylab.show()

    