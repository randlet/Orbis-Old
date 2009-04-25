import numpy
import settings

class AbstractModel(object): 
    def __init__(self): 
        self.listeners = [] 
        self.update_required = True
        
    def addListener(self, listenerFunc): 
        if listenerFunc not in self.listeners:
            self.listeners.append(listenerFunc) 
        
    def removeListener(self, listenerFunc): 
        if listenerFunc in self.listeners:
            self.listeners.remove(listenerFunc) 
        
    def update(self,event): 
        
        if self.update_required:
            for eachFunc in self.listeners: 
                eachFunc(self)
            
            self.update_required = False

class HuckelSolver(AbstractModel):
    MIN_SHAPE = (0,0)
#------------------------------------------------------------------------------------------------    
    def __init__(self,data = numpy.mat(numpy.zeros(MIN_SHAPE),float)):
        
        AbstractModel.__init__(self)
        assert type(data==numpy.matrix)
        
        self.setData(data)
#------------------------------------------------------------------------------------------------        
    def getEigenValsVecs(self):
        """Get the list of 
        Inputs: None
        Outputs: ordered list of eigenvalue, eigenvector pairs"""
        
        return zip(self.eigen_vals,self.eigen_vecs)

   
#------------------------------------------------------------------------------------------------        
    def getSize(self):
        return self.data.shape[0]
#------------------------------------------------------------------------------------------------        
    def getTotalEnergy(self):
        """Get the total energy of the system by summing population of orbital*orbital energy level
        Inputs: None
        Outputs: system energy"""
        pi_energy = 0.
        for energy,vec,pop in self.populated_levels:
            pi_energy += energy*pop
            
        return pi_energy

#------------------------------------------------------------------------------------------------    
    def reset(self):
        shape =  self.MIN_SHAPE #self.data.shape
        size = shape[0]
       
        self.num_e = size
        self.data = numpy.mat(numpy.zeros(self.MIN_SHAPE),float)
        
        self.eigen_vals = numpy.array([0]*size,float)
        self.eigen_vecs = numpy.mat(numpy.zeros(shape,float))
        self.aa_polar = []
        self.ab_polar = []
        self.unique_levels = []
        self.bond_orders =  numpy.mat(numpy.zeros(shape),float)
        self.spin_densities = numpy.array([0]*size,float)
        self.net_charges = numpy.array([0]*size,float)
        
        self.update_required = True


    def setData(self,data):
        """Set all elements of the huckel matrix and then calculate everything and trigger updater
        Inputs: data -> nxn numpy array/matrix
                note data should all be less than zero!
        Outputs: None"""
        
        assert (type(data) in (numpy.matrix,numpy.array,numpy.ndarray))
        #assert (data.shape[0] == data.shape[1] or data.shape == (1,0))
        self.reset()
        if type(data)!= numpy.matrix:
            data = numpy.mat(data)
        for ii in range(data.shape[0]):
            for jj in range(ii):
                data[ii,jj] = -1.*abs(data[ii,jj])
                data[jj,ii] = data[ii,jj]
        
        self.data = data

        self._doCalcs()
        
        
    
#------------------------------------------------------------------------------------------------    

    def setNumElectrons(self,num):
        """Set the number of electrons in the system and then update pi-bond orders, spin densities and net charges
        Inputs: num -> int, Number of electrons. Automatically set to be > 0
        Outputs: None"""
        
        
        self.num_e = max(settings.MIN_NUM_E,int(num))

        self._populateLevels()
        
        self._calcBondOrders()

        self._calcNetCharges()
        
        #self._calcAAPolarizability()
        
        #self._calcABPolarizability()
        
        self.update_required = True        
        #self.update()
        
    
#------------------------------------------------------------------------------------------------    

    def setMatrixElement(self,row,col,val):
        """Set the matrix element at data[row,col] and then calculate all and trigger updater
        Inputs: row,col -> int, int >= 0
                val -> val is converted to being a negative number automatically
        Outputs: None"""
        if row != col:
            val = -1.*abs(val)
        self.data[max(0,row),max(0,col)] = val
        self._doCalcs()
        
    
#------------------------------------------------------------------------------------------------
    def _calcAAPolarizability(self):
        """Atom-Atom polarizabilities fom Computing methods in quantum organic chemistry - Greenwood: pg 54"""
        
        N = self.getSize()
        self.aa_polar = numpy.mat(numpy.zeros((N,N),float))
        levels = self.populated_levels 
        
        dbl_orbitals = filter(lambda orb: abs(orb[2] - 2.)<settings.eps,levels)
        M = len(dbl_orbitals) #number of doubly occupied orbitals
        
        for rr in range(N):
            for uu in range(rr+1):
                aap = 0.
                for jj in range(M):
                    c_rj = self.eigen_vecs[jj][rr]
                    c_uj = self.eigen_vecs[jj][uu]
                    e_j = self.eigen_vals[jj]
                    tmp = 0
                    for kk in range(M,N):
                        c_rk = self.eigen_vecs[kk][rr]
                        c_uk = self.eigen_vecs[kk][uu]
                        e_k = self.eigen_vals[kk]
                        
                        tmp = tmp+c_rk*c_uk/(e_j-e_k)
                    aap = aap+c_rj*c_uj*tmp
                self.aa_polar[rr,uu] = aap
                self.aa_polar[uu,rr] = aap
        return 4*self.aa_polar
    
#------------------------------------------------------------------------------------------------
    def _calcABPolarizability(self):
        """Atom-Bond polarizabilities fom Computing methods in quantum organic chemistry - Greenwood: pg 46 eq 3-16"""
        
        N = self.getSize()
        self.ab_polar = []
        for uu in range(N):
            self.ab_polar.append(self._calcSingleABPolarizability(uu))
            
    def _calcSingleABPolarizability(self,uu):

        levels = self.populated_levels 
        N = self.getSize()
        
        dbl_orbitals = filter(lambda orb: abs(orb[2] - 2.)<settings.eps,levels)
        M = len(dbl_orbitals) #number of doubly occupied orbitals

        #loop over each atom and calculate
        ab_polar = numpy.mat(numpy.zeros((N,N),float))

        for ss in range(N):
            for tt in range(ss+1):
                abp = 0.
                for jj in range(M):
                    c_sj,c_tj,c_uj = [self.eigen_vecs[jj][idx] for idx in [ss,tt,uu]]
                    e_j = self.eigen_vals[jj]
                    tmp = 0
                    for kk in range(M,N):
                        c_sk,c_tk,c_uk = [self.eigen_vecs[kk][idx] for idx in [ss,tt,uu]]
                        e_k = self.eigen_vals[kk]
                        tmp = tmp+c_uk*(c_sj*c_tk+c_tj*c_sk)/(e_j-e_k)
                        
                    abp =abp+c_uj*tmp

                ab_polar[ss,tt] = abp
                ab_polar[tt,ss] = abp

        return 2*ab_polar
        
#------------------------------------------------------------------------------------------------
    def _calcBondOrders(self):
        
        """Calculate the bond order matrix"""
        size = self.getSize()
        self.bond_orders = numpy.mat(numpy.zeros((size,size),float))
        self.spin_densities = numpy.zeros(size)
        if len(self.populated_levels)>0:
            for ii in range(size):
                spin_density = 0.
                for jj in range(ii+1):
                    
                    bo = 0.
                    for val,vec,num_e in self.populated_levels:
                        bo += num_e*vec[ii]*vec[jj]
    
                        #spin density for a given atom is calculated as occupation weighted average of all SOMO coefficients squared. 
                        if num_e < 2 and ii == jj:
                            spin_density += num_e*vec[ii]*vec[ii]
                            
                    if abs(bo)<settings.eps:
                        bo =  0.0
                    self.bond_orders[ii,jj] = bo
                    self.bond_orders[jj,ii] = bo
                    
                    self.spin_densities[ii] = spin_density
    
        

#------------------------------------------------------------------------------------------------    
    def _calcNetCharges(self):
        size = self.getSize()
        self.net_charges = numpy.zeros(size)
        if self.bond_orders.any():
            for ii in range(size):
                self.net_charges[ii] = 1. - self.bond_orders[ii,ii]
        
        
#------------------------------------------------------------------------------------------------
    def _doCalcs(self):
        
        self._solveHuckel()
        
        self._populateLevels()
        
        self._calcBondOrders()
        
        self._calcNetCharges()
        
        #self._calcAAPolarizability()
        
        #self._calcABPolarizability()
        
        self.update_required = True

  
#------------------------------------------------------------------------------------------------        
    def _findUniqueLevels(self):
        eigen = self.getEigenValsVecs()
        uniques = []
        while eigen:
            eigenval,eigenvec = eigen.pop(0)

            if len(uniques) == 0:
                uniques.append( { "eigenval":eigenval, "eigenvecs":[eigenvec] })
            elif abs(eigenval-uniques[-1]["eigenval"]) < settings.eps:
                uniques[-1]["eigenvecs"].append(eigenvec)
            else:
                uniques.append( { "eigenval":eigenval, "eigenvecs":[eigenvec] })
                
        self.unique_levels = uniques
        
#------------------------------------------------------------------------------------------------        
    def _populateLevels(self):
        
        uniques = self.unique_levels

        num_e = self.num_e

        lvl_idx = 0
        populated_levels = []
        
        if len(uniques)>0:
            while lvl_idx<len(uniques):
                
                lvl = uniques[lvl_idx]
                eigenval = lvl["eigenval"]
                degeneracy = len(lvl["eigenvecs"])

                frac = min(2,num_e/float(degeneracy))
                num_e -= frac*degeneracy

                for eigenvec in lvl["eigenvecs"]:
                    populated_level = [eigenval,eigenvec,frac]
                    populated_levels.append(populated_level)
                    
                lvl_idx+=1                                                

        self.populated_levels = populated_levels
        

    def bondsExist(self):
        if 0 in self.data.shape:
            return False
        
        for ii in range(self.getSize()):
            for jj in range(ii):
                if abs(self.data[ii,jj])>settings.eps:
                    return True
        return False
    
#------------------------------------------------------------------------------------------------    
    def _solveHuckel(self):
        """Calculate the eigenvalues and vectors of the secular matrix"""
        if self.getSize()>0:#self.bondsExist():
            # get eigenvalues/vectors from lower triangular portion of secular matrix
            
            vals,vecs = numpy.linalg.eigh(self.data,'L')
    
            #switch eigenvecs to be row vectors
            vecs = vecs.T
    
            self.eigen_vals,self.eigen_vecs = self._sortEigens(vals,vecs)
        
            self._findUniqueLevels()
            
#------------------------------------------------------------------------------------------------        
    def _sortEigens(self,vals,vecs):
        """helper function that convertes the vectors to lists so that they can be zipped together with the values and then sorted by eigenvalue"""
        eigen = []
        
        for idx,val in enumerate(vals):
            vec = list(numpy.array(vecs[idx])[0])
            eigen.append([val,vec])
            
        eigen.sort()
        
        new_vals = []
        new_vecs = []
        for val,vec in eigen:
            new_vals.append(val)
            new_vecs.append(numpy.array(vec))

        return new_vals,new_vecs
    
#import molecule
#mol = molecule.Molecule(name="napthalene")
#huc = HuckelSolver(mol.huckel_matrix)
#huc._calcAAPolarizability()
#huc._calcABPolarizability()

#print huc.aa_polar[0,:]
#print ""
#print huc.ab_polar[0]


        
