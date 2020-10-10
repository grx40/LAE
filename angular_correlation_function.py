import numpy as np
from numpy.fft import fft, fftfreq, ifft, fft2, fftshift, fftn , ifftn, ifftshift
import sys
from astropy.cosmology import Planck15 as p15
from astropy.cosmology import z_at_value
from astropy.cosmology import WMAP9 as cosmo

###############################################################
##
##
##   this code takes in a halo position  box, of a path to a box, and sorts
##	the corresponding points according to their radial distances
##	 from one another. 
##
############################################################



#Useful functions for this code

def usage():
    print('Usage : halo_pos = halo_pos, DIM = DIM , Box_length = Box_length , mode = linear/power, N = N')
    print('If using power, please specify power_law = power_law')
    print('Halo_pos must be in the format of (Halo_number, [LAE count, x , y ]')


#Sort arguments
	
def main(input , DIM, Box_length, mode , N,  **kwargs):
    try:
        if input.shape[1] == 3:
            halo_pos = input
    except:
        SyntaxError('Halo Box needs to be in shape Nhalos X 3')
        
    if mode != ('linear' or 'power' ):
        raise SyntaxError('Usage: mode is either linear or power')
    else:
        mode = str(mode)
        Rmax = np.sqrt(2)*int(Box_length)
        N = int(N)
        DIM  = int(DIM)
        inter_pixel_distance = float(Box_length)/float(DIM)
        if mode == 'linear':
            bins = np.linspace(inter_pixel_distance, Rmax, N)[1:N]
        if mode == 'power':
            try:
                power_law = sys.argv[6]
            except:
                SyntaxError('You must provide power law spacing as the last cmnd line argument')


if __name__ == "__main__":
	if len(sys.argv) < 5:
		raise SyntaxError('Usage: python angualar_correlation_function.py halo_pos DIM Box_length mode N')
	else:
		print('Looks good, moving to main')
		main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])




def ACF(**kwargs):

    if __name__ != "__main__":
        #parse arguments. This is really only needed if we've imported the module
        if 'DIM' not in kwargs:
            raise SyntaxError('Usage: provide DIM = ___')
        else:
            DIM = kwargs.get('DIM')
 
        if 'Box_length' not in kwargs:
            #raise SyntaxError('Usage: provide Box_length = ___')
            print('well then')
        else:
            Box_length = kwargs.get('Box_length')
            Rmax = np.sqrt(2)*int(Box_length)
            #The distance from adjacent pixels is Box_length/DIM and represents the minimum distance
            inter_pixel_distance = float(Box_length)/float(DIM)

        if 'mode' not in kwargs:
            raise SyntaxError('Usage: provide mode = linear or power')
        else:
            mode = kwargs.get('mode')
            if mode != 'linear' and mode != 'power' and mode != 'custom':
                raise SyntaxError('Usage: mode is either linear or power or custom')
            else:
                if mode == 'linear':
                    if 'N' not in kwargs:
                        raise SyntaxError('Usage: provide number of bins N  = ___')
                    else:
                        N = int(kwargs.get('N'))
                        bins = np.linspace(inter_pixel_distance, Rmax, N)[1:N]
                if mode == 'power':
                    try:
                        kwargs.get('power_law')
                    except:
                        SyntaxError('Usage: if setting mode = power you must provide power_law = something')
                    power_law = kwargs.get('power_law')
                if mode == 'custom':
                    try:
                        kwargs.get('bins')
                    except:
                            SyntaxError('Usage: if setting mode = custom you must provide bins = something')
                    bins = kwargs.get('bins')
	
	if 'halo_pos' not in kwargs:
            usage()
            raise SyntaxError('See usage')
	else:
            halo_pos = kwargs.get('halo_pos')
            if halo_pos.shape[1] != 2:
                raise SyntaxError('The Halo_position array must be')



                 
	
	################################# INITIALIZATION OF ARRAYS ###################################
    high_resolution_interpixel_distance = 0.375

    if mode == 'power':
        bins = []
        if inter_pixel_distance > 1:
            r = inter_pixel_distance
        else:
            r = inter_pixel_distance
            while( r  <= 1):
                r += r
        while(r < Rmax):
            r = (r**power_law)
            bins.append(r)
        bins.append(Rmax)
        bins = np.array(bins)
    if mode == 'custom':
        proposed_r_bins = kwargs.get('bins')
        if np.min(proposed_r_bins) <= high_resolution_interpixel_distance:
            raise ValueError('Bins are smaller than the big box resolution')
        bins = []
#print('custom bins are', proposed_r_bins)
        for i in range(len(proposed_r_bins)):
            if proposed_r_bins[i] < Rmax:
                bins.append(proposed_r_bins[i])
#print('appending', proposed_r_bins[i] )
#print('your bins are now ' , bins)

#print('np.max(bins) and Rmax are' , np.max(bins), Rmax)
        if np.max(bins) < Rmax:
            #print('the proposed bins are the following for which we will append Rmax', bins)
            bins.append(Rmax)
        bins = np.array(bins)
            #print('Done massaging bins', bins)




    acf = np.zeros((len(bins)))
    DD = np.zeros((len(bins)))
    RR = np.zeros((len(bins)))
    DR = np.zeros((len(bins)))
    bin_counter = np.zeros_like(acf)
    

#repeats_random = 0
#repeats_data = len(halo_pos[:,0][halo_pos[:,0] > 1 ])
#print('There are ' , repeats_data , ' repeats points in the Halo Field')
    num_DD = 0
    num_RR = 0

    #lets populate the list for the random field
    if 'random_pos' in kwargs:
        RR_pos = kwargs.get('random_pos')
    else:
        #make a random poisson distribution
        RR_pos = []
        RR_Box = np.random.poisson(float((2*halo_pos.shape[0]))/float((DIM**2)),(DIM,DIM))
        for i in range(DIM):
            for j in range(DIM):
                    if RR_Box[i][j] !=0 :
                        RR_pos.append((i,j))
#RR_pos = []
#       for i in range(2*halo_pos.shape[0]):
#           xprime , yprime = (int(np.random.uniform(0,DIM)), int(np.random.uniform(0,DIM)))
#           RR_pos.append((xprime, yprime))
        RR_pos = np.array(RR_pos)
    print('all done placing', RR_pos.shape[0],  ' into a random field and ', halo_pos.shape[0], ' in the data field')

    Nd = halo_pos.shape[0]
    Nr = RR_pos.shape[0]
    print('LENGTH OF DATA ARRAY IS ' , halo_pos.shape[0] ,' and so the sum of DD should be' , float(Nd*(Nd-1))/float(2))
    print('LENGTH OF RANDOM ARRAY IS ' , Nr ,' and so the sum of DD should be' , float(Nr*(Nr-1))/float(2))
    print('LENGTH OF DATA-RANDOM ARRAY should be' , float(Nd*Nr))

    #make sure everything went to plan
    #if halo_pos.shape[0] != RR_pos.shape[0]:
    print('halo_pos shape and RR_count', halo_pos.shape[0], RR_pos.shape)
#raise ValueError('something didnt go according to plan')


################################# END INITIALIZATION OF ARRAYS ###################################

    def populate_bins(r, array):
        if mode == 'power':
            if r > bins[-1]:
                print('whoops, r', r, ' is larger than the highest bin')
                return
            #loop through, creating all r's up until Rmax
            r_floor = inter_pixel_distance
            r_ceil = bins[0]
            binn = 0
            while(r_ceil <= Rmax):
                if ((r >= r_floor) and (r < r_ceil)):
                    array[binn] += 1
                    bin_counter[binn] += r
                    return
                r_floor = r_ceil
                r_ceil = bins[binn +1]
                binn += 1

        if mode == 'linear':
            #here the 'growth' factor input is useless. Fix this in version 3.0
            r_floor = inter_pixel_distance
            delta_r = bins[0] - inter_pixel_distance
            r_ceil = bins[0]
            binn = 0
            while(r_ceil <= Rmax):
                if ((r >= r_floor) and (r < r_ceil)):
                    array[binn] += 1
                    bin_counter[binn] += r
                    return
                r_floor = r_ceil
                r_ceil = r_ceil+delta_r
                binn += 1
            if r_ceil > Rmax:
                print('r rceil Rmax' , r, r_ceil, Rmax)

                    
        if mode == 'custom':
            r_floor = 0
            r_ceil = bins[0]
            if (r > Rmax):
                print('oops, we have to cut this one out!')
                return
            binn = 0
            while(r_ceil <= np.max(bins)):
                if ((r >= r_floor) and (r < r_ceil)):
                    #if r_ceil == Rmax:
                    array[binn] += 1
                    bin_counter[binn] += r
                    return
                r_floor = r_ceil
                r_ceil = bins[int(binn)+int(1)]
                binn += 1
        print('r has tapped out! r r_floor r_ceil' , r, r_floor, r_ceil)
        raise ValueError('Something broke in populate bins!')


    #Data Data pairs
    data_pairs = []
    for i in range(halo_pos.shape[0]):
        xi, yi = ((inter_pixel_distance * halo_pos[i][0]), (inter_pixel_distance * halo_pos[i][1]) )
        for j in range(	i, halo_pos.shape[0]):
            xf, yf = ((inter_pixel_distance * halo_pos[j][0]) , (inter_pixel_distance * halo_pos[j][1]) )
            r = np.sqrt( (xf - xi)**2 + (yf - yi)**2 )
            populate_bins(r, DD)
            data_pairs.append(r)


    print('DD is ' , DD, ' but we have to subtract ' , halo_pos.shape[0] , '  zeros which shouldnt be counted. Those are the number of unique LAEs from the smallest bin')
    DD[0] = DD[0] -  halo_pos.shape[0]
    print('DD is now' , DD)
#print('DD adjusted is ' , DD)
#print('we have placed' , np.sum(DD) , ' that many points')
    avg_r_bins = np.divide(bin_counter, DD)



                

    #Data Random pairs
    print('halo_pos.shape[0] and RR_pos.shape[0]' , halo_pos.shape[0] , RR_pos.shape[0])
    for i in range(halo_pos.shape[0]):
        xi, yi =  (inter_pixel_distance * halo_pos[i][0], inter_pixel_distance * halo_pos[i][1])
        for j in range(RR_pos.shape[0]):
            xf, yf =  (inter_pixel_distance*RR_pos[j][0], inter_pixel_distance*RR_pos[j][1])
            r = np.sqrt( (xf - xi)**2 + (yf - yi)**2 )
            populate_bins(r, DR)


#print('DR is ' , DR, ' we have placed a total number of ', np.sum(DR))


    #Random Random pairs
    for i in range(RR_pos.shape[0]):
        xi, yi =  (inter_pixel_distance*RR_pos[i][0], inter_pixel_distance*RR_pos[i][1])
        for j in range(i, RR_pos.shape[0]):
            xf, yf  =  (inter_pixel_distance*RR_pos[j][0] , inter_pixel_distance*RR_pos[j][1])
            r = np.sqrt( (xf - xi)**2 + (yf - yi)**2 )
            populate_bins(r, RR)

    print('RR is ' , RR , ' we have placed a total number of ', np.sum(RR))
    print('we need to adjust RR0 by this many points: ',  RR_pos.shape[0], Nr  )
    RR[0] = RR[0] - RR_pos.shape[0]
    print('RR is now' , RR)
#print('RR is ' , RR , ' we have placed a total number of ', np.sum(RR))


    Sum_RR = np.sum(RR)
    Sum_DD = np.sum(DD)

    print('last check before finishing np.sumdr, Nr*Nd SummRR SummDD Nr Nd', np.sum(DR), Nr*Nd, Sum_RR, Sum_DD, Nr, Nd )




    #this is the old way of computing the acf
    #acf = np.divide(np.divide(DD, np.sum(DD)) - 2*np.divide(DR, np.sum(DR)) + np.divide(RR, np.sum(RR)) , np.divide(RR, np.sum(RR)))
    for i in range(len(acf)):
        if bins[i] < 5:
            print('---------------- doing bin ----------------', bins[i])
            print('DD[i] floatDD[i]/floatnp.sumDD :' , DD[i], np.sum(DD), float(float(DD[i])/float(np.sum(DD))))
            print('DR[i] -2*floatDR[i]/floatnp.sumDR :' , DR[i], np.sum(DR), float(float(-2*DR[i])/float(np.sum(DR))))
            print('RR[i] floatRR[i]/floatnp.sumRR :' , RR[i], np.sum(RR), float(float(RR[i])/float(np.sum(RR))))
            print('Nr*Nd and np.sumdr are ', Nr*Nd, np.sum(DR))
        numerator =  (float(DD[i])/float(Sum_DD)) + (float(RR[i])/float(Sum_RR)) - (float(2*DR[i]))/(float(np.sum(DR)))
        denominator = float(RR[i])/float(Sum_RR)
        if bins[i] < 5:
            print('numerator is ' , numerator)
            print('denominator is ' , denominator)
        if denominator == 0:
            #if bins[i] < 5:
            #   print('UHOH!!!!! DENOMINATOR IS 0, set RRi equal to 1/Sum(Nr) in bin ', bins[i])
            denominator = float(float(1)/float(Sum_RR))
#print('denominator is now', denominator)
            #if numerator < 0:
            #acf[i] = 0
            #continue

        acf[i] = (float(numerator)/float(denominator))
#print('acf is ', acf[i])

    bin_widths = np.zeros_like(acf)
    bin_widths[0] = (bins[0] - 0)
    for i in range(len(bin_widths)-1):
        bin_widths[i+1] = (bins[i+1] - bins[i])

    integral_acf = 0
    for i in range(len(acf)):
        integral_acf += bin_widths[i]*acf[i]
#print('bin widths are ' , bin_widths)


    return np.divide(acf, 1), bins, RR, DD, avg_r_bins, np.array(data_pairs)


def create_R_bins(mode, **kwargs):
    high_resolution_interpixel_distance = 0.375
    Rmax = np.sqrt(2)*Box_length
    if mode == 'power':
        bins = []
        if inter_pixel_distance > 1:
            r = inter_pixel_distance
        else:
            r = inter_pixel_distance
            while( r  <= 1):
                r += r
        while(r < Rmax):
            r = (r**power_law)
            bins.append(r)
            bins.append(Rmax)
            bins = np.array(bins)
    if mode == 'custom':
        proposed_r_bins = kwargs.get('bins')
        if np.min(proposed_r_bins) <= high_resolution_interpixel_distance:
            raise ValueError('Bins are smaller than the big box resolution')
        bins = []
        print('custom bins are', proposed_r_bins)
        for i in range(len(proposed_r_bins)):
            if proposed_r_bins[i] < Rmax:
                bins.append(proposed_r_bins[i])
                print('appending', proposed_r_bins[i] )
        print('your bins are now ' , bins)
            
        print('np.max(bins) and Rmax are' , np.max(bins), Rmax)
        if np.max(bins) < Rmax:
            print('the proposed bins are the following for which we will append Rmax', bins)
            bins.append(Rmax)
        bins = np.array(bins)
        print('Done massaging bins', bins)
    return bins
