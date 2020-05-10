# minimum-phase FIR filter design  
  
An approximate minimum-phase FIR filter design from specified frequency characteristic by use Hilbert transform.  
  
[github repository](https://github.com/shun60s/Python-minimum-phase-FIR-design/)  
  
  
## usage  
  
python FIR-design1.py  
draw frequency phase characteristic of approximate minimum-phase FIR filter from the frequency characteristic, freq_gain_table.csv.  
save FIR filter coefficient to a text file.  

## output sample  

frequency characteristic of freq_gain_table_csv  
![figure1](docs/1_freq_gain_table_csv__frequency_characteristic.png)  


frequency phase characteristic with approximate minimum phase by use Hilbert transform  
![figure2](docs/2_approximate_minimum_phase_frequency_characteristic.png)  


impulse waveform  
![figure3](docs/3_impulse_waveform.png)  
  
  
frequency phase characteristic of FIR length 44100 by use scipy.signal.freqz  
![figure4](docs/4_FIR_filter_length_44100_frequency_characteristic.png)  
  
  
frequency phase characteristic of FIR length 8253 by use scipy.signal.freqz  
![figure5](docs/5_FIR_filter_length_8253_frequency_characteristic.png)  
  
  
These figures are in the docs folder.  


## requirements package  
  
python3  
Please see  Check version in FIR-design1.py.  

## License  
MIT  
