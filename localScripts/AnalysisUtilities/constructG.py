import json
from numpy import *
import os
import sys


def split_in_2D(array, nrows, ncols):
    """Split a matrix into sub-matrices. Broadcasts on all dimensions except the 2 last ones"""
    n = array.shape[-2]//nrows
    r, h = array.shape[-2:]
    return (array.reshape(array.shape[:-2] + (h//nrows, nrows, -1, ncols))
                 .swapaxes(-2, -3)
                 .reshape(array.shape[:-2] + (n,n, nrows, ncols)))


#Construct the link matrix from what is available in the kind of simulation
def construct_Link(params,input_dir):
    if "LINKA" in params: #We create the total Link array
        Link = []
        with open(os.path.join(input_dir,params["LINKA"])) as file:
            LinkA = json.load(file)
        with open(os.path.join(input_dir,params["LINKN"])) as file:
            LinkN = json.load(file)

        for i in range(len(LinkN)):
            Link.append(LinkN[i].copy())
            Link[i].extend(LinkA[i])

        for i in range(len(LinkN)):
            Link.append(LinkA[i])
            Link[i+len(LinkN)].extend(LinkN[i])

    else: #The Link array is already created
        with open(os.path.join(input_dir,params["LINK"])) as file:
            Link = json.load(file)
    return Link


#This class is a template that allows one to construct the superLattice Green's function (getFullG), the periodized Green's function (PeriodizedG) and the superfluid stiffness term at (kx,ky) (get_stiffness). 
#Those quantites also depend on matsubara_frequency. In order to get the right matsubara_frequency, you have to call the function Integrand.setZ(index) where index is the number of the frequency.
#As is the Integrand class is not usable as it doesn't load the self energy. 
#In order to use the Integrand class, one should derive from it and respect a few things : 
#1. Rewrite the __init__ function. It should FIRST call the Integrand constructor with the right parameters and then define at least one variable : self.sz 
## self.zs is a 1D numpy array of the imaginary part of the matsubara frequencies (for example self.zs = array([(2*n+1)*pi/beta for n in range(0,nmax)])
#2. Write a computeSelf(self,z) function where z is the matsubara_frerquency index. It should assign to the self.selfEnergy variable the 8*8 copper self-energy (4 sites, 2 spins in Nambu Space) for the matsubara_frequency 1j*self.zs[z]


class Integrand:
    
    #The child class has to implement the conputeSelf(self,z) function (z being the matsubara frequency number)
    #It also has to create the self.zs array of matsubara frequencies
    def __init__(self,params):
        self.mu = params["mu"];self.tpd = params["tpd"]
        self.ep = params["ep"];self.tpp = params["tpp"]
        self.tppp = params["tppp"];self.beta = params["beta"]
        self.n_saves=0
        self.n_computes=0
        
    def setZ(self,z):
        self.all_data = dict()
        self.computeSelf(z)
        self.iomegan = self.zs[z]*1j + self.mu

    def set_susceptibility_q(self,q):
        self.susc_q = q
    
    def get_total_factor(self,kx,ky):
        expkx = cos(kx) + 1j*sin(kx)
        expky = cos(ky) + 1j*sin(ky)
        v = [1 + 0*kx,expkx,expkx*expky,expky]
        fact = repeat(v,3,axis=0)
        fact = moveaxis(fact,0,-1)
        fact1 = moveaxis(tile(fact,12).reshape(fact.shape[:-1]+(12,12)),-1,-2)
        fact2 = moveaxis(fact1,-1,-2).conj()
        return fact1*fact2
    
    def get_total_factor_big(self,kx,ky):
        fact1 = self.get_total_factor(kx,ky)
        fact1 = tile(fact1,(2,2))
        return fact1
   
    def get_single_epsilon(self,kx,ky):
        #kx and ky must have the same shape. Returns the hamiltonian at kx[i],ky[i] for each index i in the kx array
        expkx = cos(kx) + 1j*sin(kx)
        expmkx = cos(-kx) + 1j*sin(-kx)
        expky = cos(ky) + 1j*sin(ky)
        expmky = cos(-ky) + 1j*sin(-ky)
        G_k =[[0*kx,self.tpd*(1-expmkx),self.tpd*(1-expmky)],
              [self.tpd*(1-expkx),self.ep - 2*self.tpp+ 2*self.tppp*cos(kx),self.tpp*(1-expmky)*(1-expkx)],
              [self.tpd*(1-expky),self.tpp*(1-expky)*(1-expmkx),self.ep - 2*self.tpp+ 2*self.tppp*cos(ky)]]
        #We transpose in order to get the band indices as the last indices(simplifies the inverse after that). This works whatever the shape of kx and ky
        G_k = moveaxis(moveaxis(array(G_k),0,-1),0,-1)
        return G_k
    
    def get_epsilon(self,kx,ky):
        e_k = self.get_single_epsilon(kx,ky)
        total_factor = self.get_total_factor(kx,ky)
        e_k = tile(e_k, (4,4))*total_factor
        return e_k

    def get_epsilon_tilde(self,kx,ky):
        result =self.get_epsilon(kx,ky)
        result+=self.get_epsilon(kx+pi,ky)
        result+=self.get_epsilon(kx,ky+pi)
        result+=self.get_epsilon(kx+pi,ky+pi)
        result/=4
        return result

    def get_invG0_up(self,kx,ky):
        G_k_small = eye(3)*self.iomegan - self.get_single_epsilon(kx,ky)
        total_fact = self.get_total_factor(kx,ky)
        G_k = tile(G_k_small, (4,4))*total_fact
        return G_k

    def get_invG0_tilde_up(self,kx,ky):
        result = self.get_invG0_up(kx,ky)
        result+=self.get_invG0_up(kx+pi,ky)
        result+=self.get_invG0_up(kx,ky+pi)
        result+=self.get_invG0_up(kx+pi,ky+pi)
        result/=4
        return result

    def get_invG0_down(self,kx,ky):
        G_k_small = eye(3)*self.iomegan - self.get_single_epsilon(-kx,-ky)
        G_k = tile(G_k_small, (4,4))
        total_fact = self.get_total_factor(kx,ky)
        return -G_k.conj()*total_fact

    def get_invG0_tilde_down(self,kx,ky):
        result = self.get_invG0_down(kx,ky)
        result+=self.get_invG0_down(kx+pi,ky)
        result+=self.get_invG0_down(kx,ky+pi)
        result+=self.get_invG0_down(kx+pi,ky+pi)
        result/=4
        return result

    def getFullG(self,kx,ky):
        return linalg.inv(self.getFullInverseG(kx,ky))
    
    def getFullInverseG(self,kx,ky):
        invG0_up =  self.get_invG0_tilde_up(kx,ky)
        invG0 = zeros(invG0_up.shape[:-2] + (24,24),dtype=complex)
        invG0[...,:12,:12] = self.get_invG0_tilde_up(kx,ky)
        invG0[...,12:,12:] = -conj(self.get_invG0_tilde_up(-kx,-ky))
        invG0[...,::3,::3] -= self.selfEnergy
        return invG0
         
    def gUC(self,A,i,j):
        return A[i:i+3,j:j+3]
    
    def PeriodizedG(self,kx,ky):
        kx = array(kx)
        ky = array(ky)
        expkx = cos(kx) + 1j*sin(kx)
        expky = cos(ky) + 1j*sin(ky)
        matrix = conj(self.get_total_factor_big(kx,ky))
        
        #The formula tels us we have to put the kx and ky vectors in the reduced brillouin zone 
        kx[logical_or(kx < -pi/2,kx > pi/2)] += pi
        ky[logical_or(ky < -pi/2,ky > pi/2)] += pi
        
        fullG = self.getFullG(kx,ky)
        matrix = matrix*fullG
        matrix = split_in_2D(matrix,3,3)
        PerioG = zeros((6,6),dtype=complex)
        PerioG = zeros(fullG.shape[:-2] + (6,6),dtype=complex)
        PerioG[...,:3,:3] = sum(matrix[...,:4,:4,:,:],axis=(-4,-3))
        PerioG[...,3:,:3] = sum(matrix[...,4:,:4,:,:],axis=(-4,-3))
        PerioG[...,:3,3:] = sum(matrix[...,:4,4:,:,:],axis=(-4,-3))
        PerioG[...,3:,3:] = sum(matrix[...,4:,4:,:,:],axis=(-4,-3))
        return PerioG/4
   
              
    def construct_first_d_up(self,kx,ky):
        expkx = cos(kx) + 1j*sin(kx)
        expmkx = cos(-kx) + 1j*sin(-kx)
        expky = cos(ky) + 1j*sin(ky)
        expmky = cos(-ky) + 1j*sin(-ky)
        
        v_k = [[0*kx,1j*self.tpd*expmkx,0*kx],
               [-1j*self.tpd*expkx,-2*self.tppp*sin(kx),-1j*self.tpp*(1-expmky)*expkx],
               [0*kx,1j*self.tpp*(1-expky)*expmkx,0*kx]]
        v_k = moveaxis(moveaxis(array(v_k),0,-1),0,-1)
        return v_k

    def construct_first_d_total(self,kx,ky):
        kx = array(kx)
        ky = array(ky)
        first_d_up = self.construct_first_d_up(kx,ky)
        first_d_down = self.construct_first_d_up(-kx,-ky)
        first_d = zeros(kx.shape + (6,6),dtype=complex)
        first_d[...,:3,:3] = first_d_up
        first_d[...,3:,3:] = conj(first_d_down)
        return first_d
    
    def construct_vertex_correction(self,kx,ky):
        h=1e-4
        G_kx_ph = linalg.inv(self.PeriodizedG(kx + h,ky))
        G_kx_mh = linalg.inv(self.PeriodizedG(kx - h,ky))
        diffinvG = (G_kx_ph - G_kx_mh)/(2*h);
        diffinvG[...,:3,3:] = diffinvG[...,3:,:3] = 0
        return diffinvG
    
    def get_stiffness(self,kx,ky):
        green = self.PeriodizedG(kx,ky)
        first_d = self.construct_first_d_total(kx,ky)
        diffinvG = self.construct_vertex_correction(kx,ky)
        tau3 = diag([1]*3 + [-1]*3)

        result = -tau3 @ first_d @ green @ tau3 @ diffinvG @ green + first_d@green@diffinvG@green
        return result.trace(axis1=-2,axis2=-1)/self.beta

    
    def get_susceptibility(self,kx,ky):
        green_k = self.PeriodizedG(kx,ky)
        green_k_q = self.PeriodizedG(kx + self.susc_q[0],ky + self.susc_q[1])
        tau3 = diag([1]*3 + [-1]*3)
        result = tau3@green_k_q@green_k@tau3
        return -result[0,0]/self.beta
    
    def getStiffnessFunction(self,real):
        def inner(kx,ky):
            this_data = self.get_stiffness(kx,ky)
            if(real):
                to_return = this_data.real
            else:
                to_return = this_data.imag
            return to_return
        return inner
    
    
    def getSusceptibilityFunction(self,real):
        def inner(kx,ky):
            return self.get_susceptibility(kx,ky)
        return inner
        
    def getLen(self):
        return len(self.zs)

    
#This class is a derivative of Integrand used to load self-energies of the type used by Partick Sémon's Segment impurity Solver
#It uses external Link files to describe the names of the self Energy components and selfi.json files to save the self energy accordingly

class Patrick_Integrand(Integrand):
    
    def __init__(self,params,input_dir,selffile):
        Integrand.__init__(self,params)
        
        
        self.Link = construct_Link(params,input_dir)
        
        #List all the components of the Link file
        all_components = set([item for sublist in self.Link for item in sublist])
        with open(selffile) as file:
            json_self = json.load(file)
        
        #Store each component of the selfEnergy file in the self.indexs dict.
        self.indexs = {}
        components_length = -1
        for index in all_components:
            if index != "empty":
                if index not in json_self:
                    raise Exception("All link indices should be in in the self-energy file : " + str(index))
                else:
                    self.indexs[index] = array(json_self[index]["real"]) + 1j*array(json_self[index]["imag"])
                    if params["beta"] != json_self[index]["beta"]:
                        print("Error, beta mismatch between self and params file")
                    components_length = len(self.indexs[index])
          
        #The empty component is just a bunch of zeros
        self.indexs["empty"] = zeros(components_length,dtype=complex)
    
        #Initialization of the Matsubara frequencies
        self.zs = pi*(2*arange(0,components_length)+1)/params["beta"]  
        
        
    def computeSelf(self,z):
        self.selfEnergy = zeros((8,8),dtype=complex)
        for i in range(8):
            for j in range(8):
                self.selfEnergy[i,j] = self.indexs[self.Link[i][j]][z]
                if i>=4 and j>=4:
                    self.selfEnergy[i,j] = -self.selfEnergy[i,j].conj()

      
#This was used to load a selfarray which is already a numpy array and results were obtained at 0 temperature (thus the computation of beta inside the __init__ function.
class Sid_Integrand(Integrand):
    
    def __init__(self,mu,tpd,ep,tpp,tppp,selfarray):
        params = {"mu":mu,"tpd":tpd,"ep":ep,"tpp":tpp,"tppp":tppp}
        wmax = 2
        nmax = selfarray.shape[0]
        beta=(2*nmax+1)*pi/wmax
        params["beta"] = beta
        Integrand.__init__(self,params)
        self.zs=array([(2*n+1)*pi/beta for n in range(0,nmax)])
        self.selfarray = selfarray
        
    def computeSelf(self,z):
        self.selfEnergy = self.selfarray[z]
  
#This derivative uses a 0 selfEnergy (too complicated to solve the problem at U=0 but can be used as a test for future fonction developments)

class Zero_Integrand(Integrand):
    
     def __init__(self,mu,tpd,ep,tpp,tppp,beta,Nmax):
        params = {"mu":mu,"tpd":tpd,"ep":ep,"tpp":tpp,"tppp":tppp,"beta":beta}
        Integrand.__init__(self,mu,tpd,ep,tpp,tppp,beta)
        self.selfEnergy = zeros((8,8),dtype=complex)
        self.zs = (2*arange(Nmax) + 1)*pi/beta
    
     def computeSelf(self,z):
        pass
 