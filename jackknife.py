import numpy as np
import xi_2D




class jack_knife:
    def __init__(self, Box, DIM, L_Box, **kwargs):
        self.Box = Box
        self.DIM = DIM
        self.L_Box = L_Box
        if 'ndim' in kwargs:
            self.ndim = int(kwargs.get('ndim'))
        else:
            self.ndim = 2
    

    @staticmethod
    def jack_usage():
        return 'jack_knife(box, DIM, L_Box)'
    
    @staticmethod
    def make_subfields_usage():
        return 'provide subfields =  #subfields/DIM'

    def make_subfields(self, **kwargs):
        try:
            #make a new instance
            subfields_per_dim = float(kwargs.get('subfields_per_dim'))
            if self.DIM % subfields_per_dim != 0:
                print('We are here')
                raise SyntaxError('subfields/DIM must not have a remainder')
            self.subfields_per_dim = int(subfields_per_dim)
            print('we have made' , subfields_per_dim, self.subfields_per_dim )
        except:
            raise SyntaxError(self.make_subfields_usage())

        #now we know how many times we can cut up each dimension
        #if the array is 2D then there will be subfields_per_dim**2 subfields
        self.pixels_per_subset = int(self.DIM/self.subfields_per_dim)
        Box_subsets = np.zeros((int(self.subfields_per_dim), int(self.subfields_per_dim), int(self.pixels_per_subset), int(self.pixels_per_subset)))
        for i in range(self.subfields_per_dim):
            for j in range(self.subfields_per_dim):
                Box_subsets[i][j] = self.Box[i * self.DIM/self.subfields_per_dim: (i + 1) * self.DIM/self.subfields_per_dim, j * self.DIM/self.subfields_per_dim : (j + 1) * self.DIM/self.subfields_per_dim]

        #lets just make this subset a box an official attribute
        self.Box_subsets = Box_subsets
        return Box_subsets

            
    def compute_xi(self,  **kwargs):
        
        if 'n_bins' in kwargs:
            n_bins = kwargs.get('n_bins')
        else:
            n_bins = 100
    
        if 'bins' in kwargs:
            acfbins = kwargs.get('bins')
        else:
            acfbins = np.linspace(0.5, self.L_Box, n_bins)

        if 'Box' in kwargs:
            compute_this_box = kwargs.get('Box')
        else:
            compute_this_box = self.Box

        if 'return_extras' in kwargs:
            return_extras = kwargs.get('return_extras')
        else:
            return_extras = False

        xi, xi_bins , binn_counter, countbox = xi_2D.ps_k(1.5, 20, temperature = compute_this_box, mode = 'custom', bins = acfbins)
        
        if return_extras == True:
            return xi, xi_bins , binn_counter, countbox
        else:
            return xi
    
    

    def compute_jackknife_stat(self, **kwargs):
        #make sure whether the user made the subfields first
        try:
            subfields = self.Box_subsets
        except:
            return 'Run make_subfields method first'
        
        #sort through user changeable definitions (bins)
        if 'nbins' in kwargs:
            nbins = kwargs.get('nbins')
        else:
            nbins = 50

        if 'bins' in kwargs:
            acfbins = kwargs.get('bins')
            #how does acfbins know about the subfields DIM and L_box????
            acfbins = xi_2D.create_k_boundaries(0.008, 1.3 , mode = 'custom', bins = acfbins)[1]
        else:
            print('bins will go from 0.5 to ' + str(self.L_Box/self.subfields_per_dim) +' in ' + str(nbins) + ' steps' )
            acfbins = np.linspace(0.5, self.L_Box/self.subfields_per_dim, nbins)
            acfbins = xi_2D.create_k_boundaries(0.008, 1.3 , mode = 'custom', bins = acfbins)[1]

   
        #make sure the last two dimensions are the actual max length of the subfield boxes
        xi_subfields = np.zeros((subfields.shape[0], subfields.shape[1] , len(acfbins)))
        xi_bar = np.zeros((len(acfbins)))


        #compute the  correlation fuction for each of the subfields
        for i in range(subfields.shape[0]):
            for j in range(subfields.shape[1]):
                xi, xi_bins , binn_counter, countbox = xi_2D.ps_k(1.5, 20, temperature = subfields[i][j], mode = 'custom', bins = acfbins)
                xi_subfields[i][j] = xi
                xi_bar += xi

        #find the average at each bin
        xi_bar = np.true_divide(xi_bar, subfields.shape[0]*subfields.shape[1])
        


        #do the correlation matrix
        C_ave = np.zeros(( len(acfbins) , len(acfbins)))
        for sx in range(subfields.shape[0]):
            for sy in range(subfields.shape[1]):
                #print('Doing subfield' , sx, sy)
                #print('the ACF here is', xi_subfields[sx][sy])
                C = np.zeros(( len(acfbins) , len(acfbins)))
                for i in range(len(acfbins)):
                    for j in range(len(acfbins)):
                        #for each subfield, compute the covariance matrix
                        C[i][j] = (xi_subfields[sx][sy][i] - xi_bar[i])*(xi_subfields[sx][sy][j]- xi_bar[j])
                C_ave += C
        C_ave = np.true_divide((subfields.shape[0]*subfields.shape[1] - 1) * C_ave , (subfields.shape[0]*subfields.shape[1] - 1 ) )
        
        #normalization for 'correlation matrix'
        r = np.zeros_like(C_ave)
        for i in range(C_ave.shape[0]):
            for j in range(C_ave.shape[1]):
                r[i][j] = np.true_divide(C_ave[i][j],  np.sqrt(C_ave[i][i]*C_ave[j][j]) )

        return r


