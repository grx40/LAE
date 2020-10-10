import numpy as np
import xi_2D
from numpy.fft import fft, fftfreq, ifft, fft2, fftshift, fftn, ifftn, ifftshift


def top_hat(*args,**kwargs):
    ar0 = args[0]
    ar1 = args[1]
    if type(ar0) == np.ndarray:
        box = ar0
        k_cutoff = ar1
    else:
        box = ar1
        k_cutoff = ar0

    try:
        L_box = kwargs.get('L_box')
    except:
        print('provide L_box in kwargs')

    if 'Fourier_field' in kwargs:
        Fourier_field = np.get('Fourier_field')
    else:
        Fourier_field = False

    DIM = box.shape[0]
    ndim = box.ndim

    box_tilde = fftshift(fftn(fftshift(box)))
    box_top_filtered = np.zeros_like(box_tilde, dtype = complex)
    box_bot_filtered = np.zeros_like(box_tilde, dtype = complex)

    k_x = fftshift(fftfreq(DIM, d = float(L_box)/float(2*np.pi*DIM)))

    if ndim == 3:
        k_y = k_x
        k_z = k_x
        for a in range(DIM):
            for b in range(DIM):
                for c in range(DIM):
                    if np.sqrt(k_x[a]**2 +k_y[b]**2 + k_z[c]**2) > k_cutoff:
                        box_top_filtered[a][b][c] = box_tilde[a][b][c]
                    else:
                        box_bot_filtered[a][b][c] = box_tilde[a][b][c]

    if ndim == 2:
        k_y = k_x
        for a in range(DIM):
            for b in range(DIM):
                if np.sqrt(k_x[a]**2 +k_y[b]**2) > k_cutoff:
                    box_top_filtered[a][b] = box_tilde[a][b]
                else:
                    box_bot_filtered[a][b] = box_tilde[a][b]

    else:
        raise ValueError('This code only handles 2 or 3 dimensional boxes')

    if Fourier_field:
        return box_top_filtered, box_bot_filtered
    else:
        return ifftshift(ifftn(ifftshift(box_top_filtered))), ifftshift(ifftn(ifftshift(box_bot_filtered)))

