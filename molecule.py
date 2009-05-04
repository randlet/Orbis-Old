from __future__ import with_statement
import os,sys
import numpy
import math
import settings
import yaml

from scipy.optimize import fmin_l_bfgs_b
#===================================================================================================
class Bond():
    BOND_TYPES = {}
    FMT = '%s-%s'
    def __init__(self,a,b,k_xy=1.0):
        self.a = a
        self.b = b
        self.k_xy = k_xy

        if len(self.BOND_TYPES) == 0:
            self._readBondTypes()
            
        
    
    @property
    def atoms(self):
        return (self.a,self.b)
    
    def XCoords(self):
        return (self.a.x,self.b.x)
    def YCoords(self):
        return (self.a.y,self.b.y)
    def coords(self):
        return (self.XCoords,self.YCoords)

    @classmethod
    def AddBond(self,a,b,k):
        b1 = self.FMT % (a,b)
        b2 = self.FMT % (b,a)

        if b1 not in self.BOND_TYPES.keys() and b2 not in self.BOND_TYPES:
            Bond.BOND_TYPES[b1] = {settings.k_delim:k}
            
    @classmethod
    def GetBondType(self,a,b):
        b1 = self.FMT % (a,b)
        b2 = self.FMT % (b,a)

        if b1 in self.BOND_TYPES.keys():
            return Bond.BOND_TYPES[b1]  
        elif b2 in self.BOND_TYPES.keys():
            return Bond.BOND_TYPES[b2]  
        
    @classmethod
    def BondNames(self,types = {}):
        if len(types) == 0:
            types = self.BOND_TYPES
        return sorted(types.keys())
    
    @classmethod
    def WriteBondLib(self,bond_types):
        import yaml
        outf = open(settings.bond_data_file,'w')
        yaml.dump(bond_types,outf)
        outf.close()
        Bond._readBondTypes()
        
    @classmethod
    def _readBondTypes(self,fname=""):
        """load data for a given atom"""
        try:
            if fname =="":
                fname = settings.bond_data_file
            self.BOND_TYPES.clear()
            with file(fname,'r') as stream:
                for key,val in yaml.safe_load(stream).iteritems():
                    self.BOND_TYPES[key] = val
            
        except:
            print sys.exc_info()[0]
            #print "Unexpected error while trying to read from  " path,sys.exc_info()[0])
            raise        
    
    
        
        
#===================================================================================================
class BondDraw(Bond):
    def __init__(self,artist,*args,**kwargs):
        Bond.__init__(self,*args,**kwargs)
        self.artist = artist

    def refresh(self):
        xd = [self.a.x,self.b.x]
        yd = [self.a.y,self.b.y]
        self.artist.set_data(xd,yd)
        
#===================================================================================================        
class Atom():
    
    ATOM_TYPES = {}
    
    def __init__(self,x,y,hx=1.0,sym="C"):
        self.x = x
        self.y = y
        self.hx = hx
        self.sym = sym
        if len(self.ATOM_TYPES) == 0:
            self._readAtomTypes()
        self.coord_num = 0
    @classmethod
    def AtomNames(self,types ={}):
        if len(types) == 0:
            types = self.ATOM_TYPES
        return ["%s - %s" % (x,types[x]["description"]) for x in sorted(types.keys())]

    @classmethod
    def WriteAtomLib(self,atom_types):
        import yaml
        outf = open(settings.atomic_data_file,'w')
        yaml.dump(atom_types,outf)
        outf.close()
        Atom._readAtomTypes()
        
    
    @classmethod
    def _readAtomTypes(self,fname=""):
        """load data for a given atom"""
        try:
            if fname =="":
                fname = settings.atomic_data_file
                
            self.ATOM_TYPES.clear()
            with file(fname,'r') as stream:
                for key,val in yaml.safe_load(stream).iteritems():
                    self.ATOM_TYPES[key] = val
            
        except:
            print sys.exc_info()
            #print "Unexpected error while trying to read from  " path,sys.exc_info()[0])
            raise        
        
#===================================================================================================
class AtomDraw(Atom):
    
    def __init__(self,artist,*args,**kwargs):
        Atom.__init__(self,*args,**kwargs)
        self.artist = artist
        
    def setPos(self,x,y):
        self.x,self.y = x,y
        self.artist.set_data([x],[y])

    def setData(self,x,y,hx,sym):
        self.setPos(x,y)
        self.hx = hx
        self.sym = sym
        self.artist.set_color(settings.COLOUR_LIST[sorted(Atom.ATOM_TYPES.keys()).index(self.sym)].lower())
        
#===================================================================================================        
class Molecule():
    
    def __init__(self,atoms=[],bond_pairs=[]):
        """
        atoms is a list of tuples of the form (sym,hx)
        bond_pairs is list of tuples of the form (n_a,n_b,k_xy)
        """
        self.atom_stack = []
        self.bond_stack = []
        self.bond_pairs = []
        
    def reset(self):
        self.atom_stack = []
        self.bond_stack = []

    def addAtoms(self,atoms):
        for atom in atoms:
            
            if len(atom)>3:
                x,y,hx,sym = atom
            else:
                x,y = atom
                hx = 0.
                sym = "C"
                
            self.addAtom(x,y,hx,sym)
            
    def addAtom(self,artist,x,y,hx=0.,sym="C"):
        
        a =AtomDraw(artist,x,y,hx,sym)
        self.atom_stack.append(a)
        
        

    def addBonds(self,bonds):
        for bond in bonds:
            if len(bond)>2:
                a,b,k_xy = bond
            else:
                a,b = bond
                k_xy = 1.

            self.addBond(None,a,b,k_xy)
            
    def addBond(self,artist,a,b,k_xy=1.0):
        dup = False
        self.bond_pairs.append((a,b))
        a = self.atom_stack[a]
        b = self.atom_stack[b]
        
        for bond in self.bond_stack:
            if a in bond.atoms and b in bond.atoms:
                return
            
        a.coord_num += 1
        b.coord_num += 1
        bond = BondDraw(artist,a,b,k_xy)
        self.bond_stack.append(bond)
            
        
    def getBond(self,a,b):
        a = self.atom_stack[a]
        b = self.atom_stack[b]
        
        for bond in self.bond_stack:
            if a in bond.atoms and b in bond.atoms:
                return bond
        return None

        
    def removeAtom(self,idx):
        atom = self.atom_stack[idx]
        self.bond_stack = [bond for bond in self.bond_stack if atom not in bond.atoms]
        self.atom_stack.pop(idx)
    
    def removeBond(self,a,b):

        if a == b:
            return
        
        self.bond_pairs = [pair for pair in self.bond_pairs if not ((a,b) == pair or (b,a) == pair)]
        
        a = self.atom_stack[a]
        b = self.atom_stack[b]
        a.coord_num -= 1
        b.coord_num -= 1 
        #for idx,bond in enumerate(self.bond_stack):
            #if a in bond.atoms and b in bond.atoms:
                #bond = self.bond_stack.pop(idx)
                #break
        self.bond_stack = [bond for bond in self.bond_stack if not (a in bond.atoms and b in bond.atoms)]
        

    @property
    def numAtoms(self):
        return len(self.atom_stack)
    
    def huckelMatrix(self):
        
        huckel = numpy.zeros((self.numAtoms,self.numAtoms),float)
        
        for ii,atom in enumerate(self.atom_stack):
            #diagonal terms
            huckel[ii,ii] = atom.hx
            
            #off diagonals
            for jj in range(ii):
                bond  = self.getBond(ii,jj)
                if bond:
                    k_xy = bond.k_xy
                else:
                    k_xy = 0.0

                huckel[ii,jj] = k_xy
                huckel[jj,ii] = k_xy
                
        return huckel
    
    def updateFromHuckel(self,huckel,artist=None):
        
        for ii,atom in enumerate(self.atom_stack):
            #diagonal terms
            atom.hx = huckel[ii,ii]
            
            #off diagonals
            for jj in range(ii):
                bond  = self.getBond(ii,jj)
                if bond:
                    bond.k_xy = huckel[ii,jj] 

                    
    def _getRs(self,ps,pairs):
        #future note: converting to list comprehension was slower
        rs = []
        for ii,jj in pairs:
            diff1 = ps[ii][1]-ps[jj][1]
            diff2 = ps[ii][0]-ps[jj][0]
            r2 = math.sqrt(diff1*diff1+diff2*diff2)
            rs.append(r2)
        
        return rs
#------------------------------------------------------------------------------------------------          
    def _getR3s(self,ps,pairs):
        #future note: converting to list comprehension was slower
        rs = []
        for ii,jj in pairs:
            p1,p2 = ps[ii],ps[jj]
            diff1 = p2[1] - p1[1]
            diff2 = p2[0] - p1[0]
            r2 = math.sqrt(diff1*diff1+diff2*diff2)
            rs.append(-1.*r2*r2*r2)
        
        return rs
    
#------------------------------------------------------------------------------------------------          
    def _sumV(self,positions,A,bond,non_bond):
        xs = positions[0::2]
        ys = positions[1::2]
        ps = zip(xs,ys)
    
        rs = self._getRs(ps,bond)
        s = sum(map(lambda r: (r-A)*(r-A),rs))
        
        rs = self._getR3s(ps,non_bond)
        s += sum(map(math.exp,rs))
        
        return s
#------------------------------------------------------------------------------------------------          
    def _getInitPos(self):

        xs = [atom.x for atom in self.atom_stack]
        ys = [atom.y for atom in self.atom_stack]
        
        xys = zip(xs,ys)
        init_pos = []
        
        for pos in xys:
            init_pos.append(pos[0])
            init_pos.append(pos[1])
        
        init_pos = numpy.array(init_pos)
        
        return init_pos
#------------------------------------------------------------------------------------------------          
    def _getBondPairs(self,huck):
        bond_pairs = []
        non_bond_pairs = []
        for ii in range(huck.shape[0]):
            for jj in range(ii):
                if abs(huck[ii,jj])>0:
                    bond_pairs.append([jj,ii])
                else:
                    non_bond_pairs.append([ii,jj])
        return bond_pairs,non_bond_pairs
    
        
#------------------------------------------------------------------------------------------------          
    def minimizePositions(self):
        
        huck = self.huckelMatrix()
        self._bonding,self._non_bonding = self._getBondPairs(huck)

        max_r = 40
        A=1
        init_pos = self._getInitPos()
        
        ranges = [(-max_r,max_r)]*self.numAtoms*2
        
        self.positions = fmin_l_bfgs_b(self._sumV,init_pos,approx_grad=True,m=200,args=(A,self._bonding,self._non_bonding),bounds=ranges,maxfun=150000)[0]
        
        xs = self.positions[0::2]
        ys = self.positions[1::2]
        
        for ii,atom in enumerate(self.atom_stack):
            atom.setPos(xs[ii],ys[ii])
            
        for bond in self.bond_stack:
            bond.refresh()
            

        
#a = Atom(0,0)
#b = Atom(1,1)

#m = Molecule()
#atoms = [(0,0,1.),(1,0),(1,1),(0,1)]
#bonds = [(0,1),(1,2),(2,3),(3,0)]
#m.addAtoms(atoms)
#m.addBonds(bonds)

#for ii,atom in enumerate(m.atom_stack):
    #print ii,atom.x,atom.y
#m.minimizePositions()
#for ii,atom in enumerate(m.atom_stack):
    #print ii,atom.x,atom.y

#h = m.huckelMatrix() 

#h[0,0] = 0.5

#m.updateFromHuckel(h)
#print m.huckelMatrix()

