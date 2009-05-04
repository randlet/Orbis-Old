import os
import sys

class GeomParser():

    def __init__(self,path):
        self.atoms = []
        self.path = path
        self.name = ""
        try:
            self._stream = open(self.path,'r')
        except IOError:
            print "unable to open %s" % (self.path)
            
            
            
    def parse(self):
        pass
    
    def close(self):
        self._stream.close()
        
    def createBonds(self,bond_dist):
        bd2 = bond_dist**2
        bonds = []
        for ii,a1 in enumerate(self.atoms):
            for jj,a2 in enumerate(self.atoms):
                if ii != jj:
                    sym,x1,y1,z1 = a1
                    sym,x2,y2,z2 = a2
                    r2 = (x2-x1)**2+(y2-y1)**2+(z2-z1)**2
                    print r2**0.5,bond_dist
                    if r2 <= bd2:
                        bonds.append((ii,jj,1))
                        
        self.bonds = tuple(bonds)
        
        
    
class ParseCSV(GeomParser):
    def __init__(self,path):
        GeomParser.__init__(self,path)
        
    def parse(self,bond_dist=1.6):
        import csv
        row = self._stream.readline()
        dialect = csv.Sniffer().sniff(row)
        self._stream.seek(0)
        reader = csv.reader(self._stream,dialect)
        
        for row in reader:
            row = [x for x in row if x.strip() != '']
            if len(row) == 3:
                x,y,z = map(float,row+[0])
                sym = "A"
            elif len(row) == 4:
                x,y,z = map(float,row[1:])
                sym = row[0]
            self.atoms.append((sym,x,y,z))
        
        self.createBonds(bond_dist)
        self.close()
        
class ParseXYZ(GeomParser):
    
    def __init__(self,path):
        GeomParser.__init__(self,path)
            
    def parse(self,bond_dist=1.6):
        self.n_atoms = int(self._stream.readline().strip())
        self.name = self._stream.readline().strip()
        for atom in self._stream.readlines():
            atom = [x.strip() for x in atom.split()]

            sym = atom[0]
            x,y,z= map(float,atom[1:])
            self.atoms.append((sym,x,y,z))

        self.atoms = tuple(self.atoms)

        self.createBonds(bond_dist)
        assert(len(self.atoms)==self.n_atoms)

        


class ParseMOL(GeomParser):
    
    def __init__(self,path):
        GeomParser.__init__(self,path)
            
            
    def parse(self):
        self.name = self._stream.readline().strip()
        self.comment = self._stream.readline().strip()
        self.fmt = ""
        while self.fmt == "":
            line = self._stream.readline()
            idx = line.rfind('V')
            if idx>0:
                line = line.split()
                self.n_atoms = int(line[0].strip())
                self.n_bonds = int(line[1].strip())
                self.fmt = line[-1].strip()
        
        self.atoms = []
        for idx in range(self.n_atoms):
            line = self._stream.readline().split()

            sym = line[3]
            x,y,z= map(float,line[0:3])
            
            self.atoms.append((sym,x,y,z))

        
        self.bonds = []
        for idx in range(self.n_bonds):
            line = self._stream.readline().split()
            a,b,t = line[0:3]
            self.bonds.append((int(a)-1,int(b)-1,t))

        self.atoms = tuple(self.atoms)
        self.bonds = tuple(self.bonds)
        
        assert(len(self.atoms)==self.n_atoms)
        assert(len(self.bonds)== self.n_bonds)

        return self.atoms

        
