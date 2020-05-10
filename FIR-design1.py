#coding:utf-8

#  approximate minimum-phase FIR filter design from specified frequency characteristic by use Hilbert transform
# 

# Check version
#  Python 3.6.4 on win32 (Windows 10)
#  numpy 1.16.3
#  scipy 1.4.1
#  matplotlib 2.1.1

import argparse
import csv
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal




def get_portion(FIR_real, ratio=0.99):
	# return FIR_real portion of the ratio
    #
    sum_all= np.sum( np.abs(FIR_real))
    sum=0.
    i=0
    while sum < (sum_all * ratio):
        sum += np.abs(FIR_real[i])
        i +=1
    
    print ('ratio ', ratio, ' index ',i)
    return FIR_real[0:i].copy()


def save_txt(FIR_real, text_path='FIR_output.txt'):
    # save to text file
    with open( text_path, mode='w') as f:
        for x in FIR_real:
            print("{:.9e}".format(x), file=f)
    
    print ('save to text file ', text_path)


def plot_waveform(FIR_complex, sr, plt_show=False):
    # plot waveform
    fig = plt.figure()
    plt.subplot(211)
    tlist= np.arange( len(FIR_complex) ) * (1 /sr)
    plt.xlabel('mSec')
    plt.ylabel('level')
    plt.title('REAL part of impulse waveform')
    plt.plot( tlist, FIR_complex.real )
    
    plt.subplot(212)
    tlist= np.arange( len(FIR_complex) ) * (1 /sr)
    plt.xlabel('mSec')
    plt.ylabel('level')
    plt.title('IMAGINARY part of impulse waveform')
    plt.plot( tlist, FIR_complex.imag)
    
    fig.tight_layout()
    plt.grid()
    plt.axis('tight')
    
    if plt_show:  # At last, do plt.show()
        plt.show()


def f_show(freq, gain, title=None, phase_show=True, plt_show=False):
    # plot frequency characteristic
    
    # remove DC characteristic
    if freq[0]== 0.:
        freq=freq[1:]
        gain=gain[1:]
    
    fig = plt.figure()
    if title is not None:
        plt.title('frequency characteristic of ' + title)
    else:
        plt.title('frequency characteristic')
    ax1 = fig.add_subplot(111)
    plt.semilogx(freq, 20 * np.log10(abs(gain)), 'b') 
    plt.ylabel('Amplitude [dB]', color='b')
    plt.xlabel('Frequency [Hz]')
    
    if phase_show:
        ax2 = ax1.twinx()
        angles = np.unwrap(np.angle(gain))
        angles = angles / ((2.0 * np.pi) / 360.0)
        plt.semilogx(freq, angles, 'g')  # plt.plot(flist, angles, 'g')
        plt.ylabel('Angle(deg)', color='g')
        
    plt.grid()
    plt.axis('tight')
    
    if plt_show:  # At last, do plt.show()
        plt.show()


def is_number(s):
    # check if number float
    try:
        float(s)
    except ValueError:
        return False
    else:
        return True


def read_freq_gain_table( file_path, sr ):
    # read freq vs gain list
    f=[]
    H=[]
    with open( file_path ) as fcsv:
        reader = csv.reader(fcsv)  # header skip
        for line in reader:
            if is_number( line[0]):  # skip non number line
                if float( line[0])  > sr/2.0:
                    print ('warning: input data was ignored due to it is over than Nyquist frequency', float( line[0]) )
                else:
                    f.append( float( line[0]) )
                    H.append( float( line[1]) )
    return np.array(f), np.array(H)


if __name__ == '__main__':
    #
    parser = argparse.ArgumentParser(description=' approximate minimum-phase FIR filter design from specified frequency characteristic')
    parser.add_argument('--freq_gain_table_csv', '-c', default='freq_gain_table.csv', help='specify frequency gain table (csv format)')
    parser.add_argument('--sampling_rate', '-s', type=int, default=44100, help='sampling rate frequency (default 44100Hz)')
    args = parser.parse_args()
    
    # sampling rate frequency should be even number.
    sr= args.sampling_rate
    print('sampling rate frequency is ', sr)
    
    # freq_gain_table.csv should be
    #   frequency data points order must be increasing.
    #   DC (0 Hz) and Nyquist frequency(half of sampling rate frequency) gain data must be included.
    #   avoid 0(complete zero data) to gain data. Or use very small value instead of 0
    file_path = args.freq_gain_table_csv
    freq_raw,gain_raw= read_freq_gain_table( file_path, sr )
    print('load freq gain table from ', file_path)
    
    # show frequency characteristic of freq_gain_table_csv
    f_show(freq_raw, gain_raw, title='freq_gain_table_csv',phase_show=False)
    
    # gain is linear interpolated to number of sampling rate points
    sr_half= int(sr/2)   # sampling rate frequency sr should be even number.
    sr_new= int(sr_half * 2) # new even sampling rate
    freq_interp = np.linspace(0,sr_half, sr_half+1) 
    gain_interp= np.interp(freq_interp, freq_raw, gain_raw, left=gain_raw[0], right=gain_raw[-1]) # line interpolate
    
    # show frequency characteristic by linear_interpolate
    #f_show(freq_interp, gain_interp, title='linear interpolated', phase_show=False)
    
    # extend frequency characteristic symmetrically via Nyquist frequency
    gain_interp_reverse= gain_interp[::-1]
    gain_extend= np.concatenate( [gain_interp, gain_interp_reverse[:-2] ], 0)
    #freq_extend = np.linspace(0,sr_new-1, sr_new)
    
    # hilbert transfer  by FFT
    # H(z)=ln(|H(z)|) + j arg(H(z))
    # log and ifft
    lnH_list=np.log(gain_extend) # ln(|H(z)|) 
    IFFT0=np.fft.ifft(lnH_list)
    
    # h is like hilbert filter
    h= np.zeros(sr_new)
    h[0]=1.0
    h[int(sr_new/2)]=1.0
    h[1:int(sr_new/2)]=2.0
    
    # get minimum phase response
    FFT0=np.fft.fft( IFFT0 * h )  # imag(FFT0) is  approximate minimum phase
    
    gain_extend_complex = gain_extend * np.exp( 1j * FFT0.imag)
    
    # show frequency characteristic with approximate minimum phase
    f_show(freq_interp, gain_extend_complex[:int(sr_half+1)], title='approximate minimum phase')
    
    # get impulse (= FIR filter coefficient, only real part) from frequency response with approximate minimum phase via ifft
    Impulse_complex=np.fft.ifft(gain_extend_complex) # complex impulse response
    FIR_real=Impulse_complex.real.copy() # get only real part
    
    # show impulse waveform
    plot_waveform( Impulse_complex, sr)
    
    # get frequency characteristic of FIR filter
    wlist, fres = signal.freqz(FIR_real, worN=int(sr_half+1))
    # show frequency characteristic of FIR filter
    f_show(wlist * sr_half, fres, title='FIR filter' + ' length ' + str(len( FIR_real)) )
    # save FIR filter coefficient as text file
    save_txt(FIR_real, text_path='FIR_output_' + str(len( FIR_real)) + '.txt')
    
    
    # get some portion (ratio) of FIR filter
    FIR_real_portion= get_portion(FIR_real, ratio=0.99)
    
    # get frequency characteristic of some portion (ratio) of FIR filter
    wlist, fres = signal.freqz(FIR_real_portion, worN=int(sr_half+1))
    # show frequency characteristic of some portion (ratio) of FIR filter
    f_show(wlist * sr_half, fres, title='FIR filter' + ' length ' + str(len( FIR_real_portion)) )
    # saves ome portion (ratio) of FIR filter coefficient as text file
    save_txt(FIR_real_portion, text_path='FIR_output_' + str(len( FIR_real_portion)) + '.txt')
    
    
    # At last, do plt.show() to show all figures
    plt.show()
    
