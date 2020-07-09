import json
from numpy import *

LinkAfile = "scripts/LinkA.json"
LinkNfile = "scripts/LinkN.json"

def split_in_2D(array, nrows, ncols):
    """Split a matrix into sub-matrices."""
    n = array.shape[0]//nrows
    r, h = array.shape
    return (array.reshape(h//nrows, nrows, -1, ncols)
                 .swapaxes(1, 2)
                 .reshape(n,n, nrows, ncols))


#This class is a template that allows one to construct the superLattice Green's function (getFullG), the periodized Green's function (PeriodizedG) and beta*(the superfluid stiffness term) at (kx,ky) (get_stiffness). 
#Those quantites also depend on matsubara_frequency. In order to get the right matsubara_frequency, you have to call the function Integrand.setZ(index) where index is the number of the frequency.
#As is the Integrand class is not usable as it doesn't load the self energy. 
#In order to use the Integrand class, one should derive from it and respect a few things : 
#1. Rewrite the __init__ function. It should FIRST call the Integrand constructor with the right parameters and then define at least one variable :  
## - self.zs which is a 1D numpy array of the matsubara frequencies (for example self.zs = array([(2*n+1)*pi/beta for n in range(0,nmax)])
#2. Write a computeSelf(self,z) function where z is the matsubara_frerquency index. It should assign to the self.selfEnergy variable the 6*6 self-energy (3 bands, 2 spins in Nambu Space) for the matsubara_frequency self.zs[z]


class Integrand:
    
    #The child class has to implement the conputeSelf(self,z) function (z being the matsubara frequency number)
    #It also has to create the self.zs array of matsubara frequencies
    def __init__(self,mu,tpd,ep,tpp,tppp,beta):
        self.mu = mu;self.tpd = tpd;self.ep = ep;self.tpp = tpp;self.tppp = tppp;
        self.n_saves=0
        self.n_computes=0
        
    def setZ(self,z):
        self.all_data = dict()
        self.computeSelf(z)
        self.iomegan = self.zs[z]*1j + self.mu

        
    def get_total_factor(self,kx,ky):
        expkx = cos(kx) + 1j*sin(kx)
        expky = cos(ky) + 1j*sin(ky)
        v = [1,expkx,expkx*expky,expky]
        fact = repeat(v,3)
        fact1 = zeros((12,12),dtype=complex)
        fact2 = zeros((12,12),dtype=complex)
        fact1[:,0] = fact
        fact2[0,:] = fact.conj()
        total_fact = dot(fact1,fact2)
        return total_fact
    
    def get_total_factor_big(self,kx,ky):
        fact1 = self.get_total_factor(kx,ky)
        c = zeros((24,24),dtype=complex)
        c[:12,:12] =   c[:12,12:] = c[12:,:12] = c[12:,12:] = fact1
        return c

    def get_single_epsilon(self,kx,ky):
        expkx = cos(kx) + 1j*sin(kx)
        expmkx = cos(-kx) + 1j*sin(-kx)
        expky = cos(ky) + 1j*sin(ky)
        expmky = cos(-ky) + 1j*sin(-ky)
        G_k = zeros((3,3),dtype=complex); G_k[0,1] = self.tpd*(1-expmkx) ; G_k[0,2] = self.tpd*(1-expmky)
        
        G_k[1,0] =  self.tpd*(1-expkx) ; G_k[1,1] = self.ep - 2*self.tpp+ 2*self.tppp*cos(kx);
        G_k[1,2] = self.tpp*(1-expmky)*(1-expkx)
        
        G_k[2,0] = self.tpd*(1-expky); G_k[2,1] = self.tpp*(1-expky)*(1-expmkx) ;
        G_k[2,2] =  self.ep - 2*self.tpp+ 2*self.tppp*cos(ky)
        return G_k

    def get_epsilon(self,kx,ky):
        e_k = self.get_single_epsilon(kx,ky)
        total_fact = self.get_total_factor(kx,ky)
        e_k = tile(e_k, (4,4))*total_fact
        return e_k

    def get_epsilon_tilde(self,kx,ky):
        result = zeros((12,12),dtype=complex)
        result+=self.get_epsilon(kx,ky)
        result+=self.get_epsilon(kx+pi,ky)
        result+=self.get_epsilon(kx,ky+pi)
        result+=self.get_epsilon(kx+pi,ky+pi)
        result/=4
        return result

    def get_invG0_up(self,kx,ky):
        G_k_small = eye(3)*self.iomegan - self.get_single_epsilon(kx,ky)
        G_k = tile(G_k_small, (4,4))
        total_fact = self.get_total_factor(kx,ky)
        return G_k*total_fact

    def get_invG0_tilde_up(self,kx,ky):
        result = zeros((12,12),dtype=complex)
        result+=self.get_invG0_up(kx,ky)
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
        result = zeros((12,12),dtype=complex)
        result+=self.get_invG0_down(kx,ky)
        result+=self.get_invG0_down(kx+pi,ky)
        result+=self.get_invG0_down(kx,ky+pi)
        result+=self.get_invG0_down(kx+pi,ky+pi)
        result/=4
        return result

    def getFullG(self,kx,ky):
        invG0 = zeros((24,24),dtype=complex)
        invG0[:12,:12] = self.get_invG0_tilde_up(kx,ky)
        invG0[12:,12:] = -conj(self.get_invG0_tilde_up(-kx,-ky))
        invG0[::3,::3] -= self.selfEnergy
        return linalg.inv(invG0)
    
    def gUC(self,A,i,j):
        return A[i:i+3,j:j+3]
    
    def PeriodizedG(self,kx,ky):
        expkx = cos(kx) + 1j*sin(kx)
        expky = cos(ky) + 1j*sin(ky)
        PerioG = zeros((6,6),dtype=complex)
        fullG = self.getFullG(kx,ky)
        if kx < pi/2 or kx > pi/2:
            kx+=pi
        if ky < pi/2 or ky > pi/2:
            kx+=pi
        v = [1,expkx,expkx*expky,expky]
        matrix = conj(self.get_total_factor_big(kx,ky))*fullG
        matrix = split_in_2D(matrix,3,3)
        PerioG[:3,:3] = sum(matrix[:4,:4],axis=(0,1))
        PerioG[3:,:3] = sum(matrix[4:,:4],axis=(0,1))
        PerioG[:3,3:] = sum(matrix[:4,4:],axis=(0,1))
        PerioG[3:,3:] = sum(matrix[4:,4:],axis=(0,1))
        return PerioG/4
            
    def construct_first_d_up(self,kx,ky):
        expkx = cos(kx) + 1j*sin(kx)
        expmkx = cos(-kx) + 1j*sin(-kx)
        expky = cos(ky) + 1j*sin(ky)
        expmky = cos(-ky) + 1j*sin(-ky)
        v_k = zeros((3,3),dtype=complex); v_k[0,1] = 1j*self.tpd*expmkx ;v_k[0,2] = 0
        v_k[1,0] =  -1j*self.tpd*expkx ; v_k[1,1] = -2*self.tppp*sin(kx);v_k[1,2] = -1j*self.tpp*(1-expmky)*expkx
        v_k[2,0] = 0; v_k[2,1] =  1j*self.tpp*(1-expky)*expmkx ;v_k[2,2] =  0
        
        """
        v_k[0,0] +=0; v_k[0,1] += 1/2*1j*self.tpd*(1-expmkx) ;v_k[0,2] += 0
        v_k[1,0] += -1/2*1j*self.tpd*(1-expkx) ; v_k[1,1] +=0;v_k[1,2] += -1/2*1j*self.tpp*(1-expmky)*(1-expkx);
        v_k[2,0] += 0; v_k[2,1] +=  1/2*1j*self.tpp*(1-expky)*(1-expmkx);v_k[2,2] +=  0
        """
        
        return v_k

    def construct_first_d_total(self,kx,ky):
        first_d_up = self.construct_first_d_up(kx,ky)
        first_d_down = self.construct_first_d_up(-kx,-ky)
        first_d = zeros((6,6),dtype=complex)
        first_d[:3,:3] = first_d_up
        first_d[3:,3:] = conj(first_d_down)
        return first_d
    
    def construct_vertex_correction(self,kx,ky):
        h=1e-4
        G_kx_ph = linalg.inv(self.PeriodizedG(kx + h,ky))
        G_kx_mh = linalg.inv(self.PeriodizedG(kx - h,ky))
        diffinvG = (G_kx_ph - G_kx_mh)/(2*h);
        diffinvG[:3,3:] = diffinvG[3:,:3] = 0
        return diffinvG
    
    def get_stiffness(self,kx,ky):
        green = self.PeriodizedG(kx,ky)
        first_d = self.construct_first_d_total(kx,ky)
        tau3 = diag([1]*3 + [-1]*3)
        diffinvG = self.construct_vertex_correction(kx,ky)
        
        result = first_d @ tau3 @ green @ diffinvG @tau3 @ green - first_d@green@diffinvG@green

        return -result.trace()
     
    def getStifnessFunction(self,real):
        def inner(kx,ky):
            try:
                if(real):
                    to_return = self.all_data[kx,ky].real
                else:
                    to_return = self.all_data[kx,ky].imag
                self.n_saves+=1;
            except:
                self.all_data[kx,ky] = self.get_stiffness(kx,ky)
                if(real):
                    to_return = self.all_data[kx,ky].real
                else:
                    to_return = self.all_data[kx,ky].imag
                self.n_computes+=1;
            return to_return
        return inner
    
    def getLen(self):
        return len(self.zs)


class Patrick_Integrand(Integrand):
    
    def __init__(self,mu,tpd,ep,tpp,tppp,beta,selffile):
        Integrand.__init__(self,mu,tpd,ep,tpp,tppp,beta)
        file = open(LinkAfile)
        self.LinkA = json.load(file)
        file.close()
        file = open(LinkNfile)
        self.LinkN = json.load(file)
        file.close()
        selfEnergytxt = loadtxt(selffile)
        pphi = selfEnergytxt[:,7] + 1j*selfEnergytxt[:,8]
        mphi = selfEnergytxt[:,9] + 1j*selfEnergytxt[:,10]
        s00 = selfEnergytxt[:,1] + 1j*selfEnergytxt[:,2]
        s01 = selfEnergytxt[:,3] + 1j*selfEnergytxt[:,4]
        s11 = selfEnergytxt[:,5] + 1j*selfEnergytxt[:,6]
        self.zs = selfEnergytxt[:,0]
        self.indexs = dict()
        self.indexs["empty"] = len(pphi)*[1j*0.]
        self.indexs["pphi"] = (pphi - mphi)/2
        self.indexs["mphi"] = (mphi - pphi)/2
        self.indexs["00"] = s00
        self.indexs["01"] = s01
        self.indexs["11"] = s11
    
    def computeSelf(self,z):
        self.selfEnergy = zeros((8,8),dtype=complex)
        for i in range(4):
            for j in range(4):
                self.selfEnergy[i,j] = self.indexs[self.LinkN[i][j]][z]
                self.selfEnergy[i+4,j] = self.indexs[self.LinkA[i][j]][z]
                self.selfEnergy[i,j+4] = self.indexs[self.LinkA[i][j]][z]
                self.selfEnergy[i+4,j+4] = -self.indexs[self.LinkN[i][j]][z].conj()


class Sid_Integrand(Integrand):
    
    def __init__(self,mu,tpd,ep,tpp,tppp,selfarray):
        wmax = 2
        nmax = selfarray.shape[0]
        beta=(2*nmax+1)*pi/wmax
        print("beta = " + str(beta)) 
        Integrand.__init__(self,mu,tpd,ep,tpp,tppp,beta)
        self.zs=array([(2*n+1)*pi/beta for n in range(0,nmax)])
        self.selfarray = selfarray
        
    def computeSelf(self,z):
        self.selfEnergy = self.selfarray[z]
        
