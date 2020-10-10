import time
import numpy as np
import sys
from decimal import *
from numpy.fft import fft, fftfreq, ifft, fft2, fftshift, fftn, ifft, ifftn, ifftshift
from subprocess import *
import os


###############################
#      Useful Functions     #
###############################

def create_r_boundaries(DIM, L_box ,*args, **kwargs):
    r_x = L_box*fftshift(fftfreq(DIM))
    rmax = np.sqrt(2)*np.max(np.abs(r_x))
    
    mode = kwargs.get('mode')
    if mode == 'power':
        r_end_of_first_bin = args[0]
        r_growth_factor = args[1]
    if mode == 'linear':
        r_end_of_first_bin = args[0]

    numbins = 0
    r_floor = 0
    r_range = []
    
    if mode == 'power':
        #make the bins for k with growth factor
        r_end_of_first_bin = float(r_end_of_first_bin)
        r_ceil = r_end_of_first_bin
        while(r_ceil< rmax):
            numbins += 1
            r_range.append(r_ceil)
            r_floor = r_ceil
            r_ceil = r_ceil*r_growth_factor

    if mode == 'linear':
        
        #make the bins for k with linear addition
        r_end_of_first_bin = float(r_end_of_first_bin)
        r_ceil = r_end_of_first_bin
        while(r_ceil< rmax):
            numbins += 1
            r_range.append(k_ceil)
            #update floor and ceiling
            r_floor = r_ceil
            r_ceil = r_ceil+r_end_of_first_bin

    if mode == 'custom':
        custom_r_bins = kwargs.get('bins')
        #print('custom_k_bins are ', custom_k_bins)
        #print('there are this many bins:', len(custom_k_bins))
        #print('kmax and bins are :', kmax, custom_k_bins)
        for i in range(len(custom_r_bins)):
            if custom_r_bins[i] < rmax:
                r_range.append(custom_r_bins[i])
                numbins += 1
            #print('k_range:', k_range)
            #print('now k range is', k_range)
            #add while here
            #while(np.max(k_range) < kmax):
            #execute
        if len(r_range) < len(custom_r_bins) and np.max(r_range) < rmax:
            #print('max of k range is', np.max(k_range))
            #print('len of k_range and custom_k_bins is ', len(k_range), len(custom_k_bins))
            #print('appending ', custom_k_bins[len(k_range)])
            r_range.append(custom_r_bins[len(r_range)])
            numbins += 1

    return numbins, r_range




def z2f(z):
    F21 = 1.42040575177
    return F21/(z + 1)


def load_binary_data(filename, dtype=np.float32):
    f = open(filename, "rb")
    data = f.read()
    f.close()
    _data = np.fromstring(data, dtype)
    if sys.byteorder == 'big':
        _data = _data.byteswap()
    return _data

def oops_we_are_at_the_edge(y,z,HII_DIM):
    if y == HII_DIM and z == HII_DIM:
        return int(HII_DIM - 1), int(HII_DIM - 1)
    else:
        if y == HII_DIM:
            return int(HII_DIM - 1) , z
        if z == HII_DIM:
            return y, int(HII_DIM - 1)


#Check to see if we are running from console

if __name__ == "__main__":
    #let us check to make sure the arguments are properly inputted
    if len(sys.argv) < 4:
        raise IOError('Did not provide either of Box, L_Box, DIM, bins')
    else:
    #the right number of arguments were provided. Sort them through
        Box = sys.argv[1]
        L_box = sys.argv[2]
        DIM = sys.argv[3]
        bins = sys.argv[4]
        #make sure they were provided in the right order
        if Box.ndim != 2:
            raise ValueError('The box is in the wrong format!')
        if type(bins) == list:
            bins = np.array(bins)
        if bins.dim != 1 or len(bins) < 1:
            raise ValueError('fourth arg (bins) are incorrectly provided')

                    
def import_usage():
    return 'DIM = __, L_box = __, Box = __ , mode = __, if mode is linear, provide r_end_of_first_bin , if mode is power, additionally provide r_growth_factor, if mode = custom, provide bins '
                      
                      
                      

                      
        
def  xi_avg(**kwargs):
    #lets time this, because why not?!
    start = time.time()
    #if we're running this from the console then lets skip this part since it has already been done
    if __name__ != "__main__":
    #get the keys from the dictionary
        try:
            keys = kwargs.keys()
            #take in all user specified data and store in a dictionary for now
            user_input = {}
            for i in range(len(keys)):
                user_input[keys[i]] = kwargs.get(keys[i])
            DIM = user_input['DIM']
            L_box = user_input['L_box']
            Box = user_input['Box']
            mode = user_input['mode']
            if mode == 'custom':
                bins = user_input['bins']
                custom_r_range = bins = user_input['bins']
                custom_r_max = np.max(custom_r_range)
            if mode == 'linear':
                r_end_of_first_bin = user_input['r_end_of_first_bin']
            if mode == 'power':
                r_end_of_first_bin = user_input['r_end_of_first_bin']
                r_growth_factor = user_input['r_growth_factor']
                
        except:
            import_usage()
            raise SyntaxError('provide kwargs, see usage')

        
    
    #r axes for the binning phase
    r_x = L_box*fftshift(fftfreq(DIM))
    r_y = r_x
    rmax = np.sqrt(2)*np.max(np.abs(r_x))
    print('Rmax for this box is ', rmax)



    def populate_bins(i,j, r):
        if r >= np.max(r_range):
            #print('we are throwing out', k)
            return
        if mode == 'power':
            #loop through, creating all k's up until kmax
            r_floor = 0
            r_ceil = r_end_of_first_bin
            binn = 0
            while(r_ceil < rmax):
                if ((r >= r_floor) and (r < r_ceil)):
                    acf[binn] += xi[i][j]
                    bin_counter[binn] += 1
                    break
                r_floor = r_ceil
                r_ceil = r_ceil*r_growth_factor
                binn += 1
                
        if mode == 'linear':
            #here the 'growth' factor input is useless. Fix this in version 3.0
            r_floor = 0
            r_ceil = r_end_of_first_bin
            binn = 0
            while(r_ceil < kmax):
                if ((r >= k_floor) and (r < r_ceil)):
                    acf[binn] += xi[i][j]
                    bin_counter[binn] += 1
                    break
                r_floor = r_ceil
                r_ceil = r_ceil+r_end_of_first_bin
                binn += 1

                
        if mode == 'custom':
            r_floor = 0
            r_ceil = custom_r_range[0]
            if (r > custom_r_max):
                #print('oops, we have to cut this one out!')
                #print('k and custom k max are' , k, custom_k_max)
                return
            binn = 0
            while(r_ceil <= custom_r_max):
                if ((r >= r_floor) and (r < r_ceil)):
                    acf[binn] += xi[i][j]
                    bin_counter[binn] += 1
                    break
                r_floor = r_ceil
                r_ceil = custom_r_range[binn+1]
                binn += 1



    
    print('loading the box')
    if type(Box) == str:

        arr = load_binary_data(Box)
        arr.shape = (DIM,DIM,DIM)
        Box_data = arr.reshape((DIM,DIM, DIM), order = 'F')
        Box_data = np.subtract(Box_data, np.full((DIM, DIM, DIM), np.mean(Box_data)))
        raw_transform = fftshift(fftn(fftshift(Box_data[0])))
        ps = raw_transform*np.conjugate(raw_transform)
        xi = np.real(fftshift(ifftn(fftshift(ps))))
    if type(Box) == np.ndarray:
        if Box.shape[1] != DIM:
            print('looks like we got a list array of shape ', Box.shape )
            Box_data = np.zeros((DIM,DIM))
            print('Num halos in this box : ', Box.shape[0])
                
            for i in range(Box.shape[0]):
                x , y = int(np.round(Box[i][0],0)) ,int(np.round(Box[i][1],0))
                if x == DIM or y == DIM:
                    x, y  = oops_we_are_at_the_edge(x,y, DIM)
                Box_data[x][y] += 1
        else:
            Box_data = Box
            #should put a check here to make sure it is the right shape (DIM, DIM)

    if type(Box) == list:
        Box_data = np.zeros((DIM,DIM))
        for i in range(len(Box)):
            x , y = int(np.round(Box[i][0],0)) ,int(np.round(Box[i][1],0))
            if x == DIM or y == DIM:
                x, y  = oops_we_are_at_the_edge(x,y, DIM)
            Box_data[x][y] += 1
                

    #here we divide by the mean

    Nd = np.sum(Box_data)
    Np = float(Nd*(Nd-1))/float(2)
    np_density = float(np.sum(Box_data))/float(DIM**2)

    n_density = float(np.sum(Box_data))/float(DIM**2)
    n_density = np_density

    print('n_density is ' , n_density)
    test_data = Box_data
    test_data = np.true_divide(test_data , n_density)
    test_data = np.subtract(test_data , 1)

    #temperature_data = np.divide( temperature_data * Rmax**2, np.sum(temperature_data))
    #temperature_data = np.subtract(temperature_data, np.full((DIM, DIM), np.mean(temperature_data)))
        
    Box_data = np.true_divide( Box_data, float(np.sum(Box_data))/float(DIM**2) )
    Box_data = np.subtract(Box_data, np.full((DIM, DIM),1) )


    #checks
    print('Checking that the field has been properly normalized - should be zero sum')
    print(np.sum(Box_data))
    #print(np.max(temperature_data), np.min(temperature_data))

        
    raw_transform = fftshift(fftn(fftshift(Box_data)))
    ps = raw_transform*np.conjugate(raw_transform)

    #remove the shot noise
#print('Check that we have subtracted correctly : ps[3][3] n_density',ps[3][3], n_density )
    ps = np.subtract( np.true_divide(ps, 1), float(1.)/float(n_density))
#print('Check that we have subtracted correctly : ps[3][3] n_density',ps[3][3], 1./n_density )

    xi = np.real(fftshift(ifftn(fftshift(ps))))
 
 
 
    #determine number of bins
    if mode == 'custom':
        num_bins, r_range = create_r_boundaries(DIM, L_box, mode = mode , bins = custom_r_range )
    else :
        if mode == 'linear':
            num_bins, r_range = create_r_boundaries(DIM, L_box, r_end_of_first_bin, mode = mode )
        if mode == 'power':
            num_bins, r_range = create_r_boundaries(DIM, L_box, r_end_of_first_bin, r_growth_factor, mode = mode )

    
    #initialize arrays
    acf = np.zeros(num_bins)
    bin_counter = np.zeros(num_bins)

    
    #cycle through all magnitude of k values determine which bin this k falls into
    rsum=0
    rsum2=0
    rstuff = []
    for i in range(DIM):
        rsum = r_x[i]**2
        for j in range(DIM):
            rsum2 = np.sqrt(rsum+ r_y[j]**2)
            rstuff.append(rsum2)
            populate_bins(i,j, rsum2)

    print('max of r range is :, ' , np.max(rstuff))


    
    #average out over shells of constant r
    for i in range(num_bins):
        if bin_counter[i] == 0:
            acf[i] = 0
        else:
            acf[i] = float(acf[i])/(float(bin_counter[i]*DIM**2))

    #not a great coding practice but this will make things easier for me to plot
    global r_bins
    r_bins = r_range
    
    print('Runtime (s): ', (time.time() - start))
    return acf , r_bins , bin_counter , Box_data
    

