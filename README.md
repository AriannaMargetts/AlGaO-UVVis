# AlGaO-UVVis
Calculating bandgap and aluminium concentration of (-201) β-(AlxGa1-x)2O3 grown on c-plane sapphire using UV-Vis data

UV-Vis Analysis code
12/7/23

Code was written in Spyder IDE 5.1.5, included in Anaconda. 
Python 3.8.12

Code is designed to calculate bandgap and aluminium concentration of AlGaO samples through UV-Vis data, as well as plot Tauc plots. Code can be used in theory to plot Tauc plots and find bandgaps of any sample, but ignore Al% outputs. 

Data input:

- Code will take .csv file of wavelength (nm) and transmittance (T%) data taken from UV-Vis measurements.
- The .csv file name will need to be manually inputted into the code between the imports and variables.
- The .csv file must have the same filepath as the program.
- The .csv file content must have the following format: wavelength (nm) in the first column, transmittance (T%) in the second column, with NO COLUMN TITLES OR UNITS, just data.
- You must know the thickness of the thin-film. Take measurement using reflectometry if needed.

Follow console prompts to input thin-film thickness and sample name. Program will output two plots - the first is a Tauc plot of all data points, the second is a Tauc plot of the linear section of Tauc plot used to calculate bandgap and Al%. 

Console will print the bandgap, and two Al% values, based on two sources: one from consultation with Agnitron, the other from a paper by Miller et al: Miller, Ross & Alema, Fikadu & Osinsky, Andrei. (2018). Epitaxial β-Ga2O3 and β-(AlxGa1-x)2O3/β-Ga2O3 Heterostructures Growth for Power Electronics. IEEE Transactions on Semiconductor Manufacturing. PP. 1-1. 10.1109/TSM.2018.2873488. 
