# moogplusplus
Hi. This is very much a readme-in-progress. Good Luck.

Named for its expansion on the functionality of MOOG, as well as the primary language in which it was written, `moogplusplus` is designed to be a highly accurate line-by-line stellar spectrum fitting code fit for use with as many wavelength ranges and stellar spectral types as possible.

Currently, work is being done to add multithreading support, and future work will focus on adding more line lists, especially for the near-infrared wavelength range, as support is extended to late-K and M type stars.

# Installation

## 1. Manage Dependencies

### Python
- astroquery
- astropy
- scipy
- numpy

If you you are missing any of these, simply run
`pip install astroquery astropy scipy numpy`

### C++

- **GNU Scientific Library (`libgsl`)**:
  - On my Ubuntu-based machine, I was simply able to run `sudo apt install libgsl-dev` and have everything done for me. On the other hand, if your favorite package manager does not have a GSL package ready, or you do not have write access to `/usr/include`, you will need to compile this package from source. Doing so in most cases should install GSL in the default (correct) location


### External data and code
**MOOG++ WILL NOT WORK without each of these having been installed properly.** None of these are strictly required for compiling the C++ or even running the Python scripts, but they are absolutely necessary if you want to actually do science

- **MOOG**:
  - MOOG is now included in moogplusplus. Enter the "moog" directory and simply run `make` from the command line. This is all the setup that will be needed for MOOG itself. The particular version of moog used here was adapted slightly from that provided with [PyMOOGi](https://github.com/madamow/pymoogi).

- **linemake**:
  - The ([Linemake](https://github.com/vmplacco/linemake)) source code itself is now included with moogplusplus. Compile linemake by entering the following `f95` command into the command line while inside the linemake directory: `f95 linemake.f -o linemake -O2 -ffixed-line-length-none`

- **MIST stellar isochrone models**:
  - These files were downloaded and extracted from [MIST](https://waps.cfa.harvard.edu/MIST/interp_isos.html). A direct download link can be found [here](https://drive.google.com/file/d/1WLZmoJvy1Zaznw7W1Q54UEueMSIFL2zE/view?usp=sharing) (If this link is broken, please contact me or open an issue). Extract the contents of the "isochrones" folder inside to the "isochrones" folder inside the MOOG++ directory

- **PHOENIX model atmospheres and limb-darkening profiles**:
  - These files were downloaded and extracted from the [Goettingen Spectral Library](https://phoenix.astro.physik.uni-goettingen.de/). A direct download link can be found [here](https://drive.google.com/file/d/1XWYMXebAfvxdYTqsoHW4pzrw9-OXVhOq/view?usp=sharing) (If this link is broken, please contact me or open an issue). Extract the contents of the "PHOENIX-new" folder inside to the "PHOENIX" folder inside the MOOG++ directory

- **Linemake's "MOOGlists" database of atomic and molecular lines**:
  - A direct download link can be found [here](https://drive.google.com/file/d/1agxDQRdBYV_SX4kKDTBqFkb3DIVOjnoU/view?usp=sharing) (If this link is broken, please contact me or open an issue). Extract the contents inside to the "mooglists" folder inside the MOOG++ directory


Once you've done everything listed above, 

## 2. Compile the C++ code
I have only used g++ for this. I make no promises as to how using other compilers will work
### Default case: GSL installed normally
- Assuming nothing special was done, GSL should have installed itself in /usr/local/lib, where all other C/C++ libraries are located. In this case, run the following:
  - `g++  RunAbundanceOnGoodLines.cpp -o RunAbundanceOnGoodLines -lgsl -lgslcblas -O2`

### Special case: GSL not installed in default location
- In this case you'll have to tell your compiler where exactly GSL is located, as indicated with the `-I` and `-L` compiler flags:
  - `g++  RunAbundanceOnGoodLines.cpp -o RunAbundanceOnGoodLines -I[/path/to/your/libaries]/include -L[/path/to/your/libaries]/lib -lgsl -lgslcblas -O2`




# Running MOOG++

## 1. Prepare the Stellar Spectrum
The observed stellar spectrum file **absolutely must:**
- Be normalized (continuum flux = 1 at all wavelengths)
- Be corrected for stellar RV and barycentric velocity 
- Be stored in plain-text format with one data point per line, in the format of `wavelength flux errorbar`. (it's fine if you don't have error bars, just make that column all zeros)
- Have wavelengths in air wavelengths (not vacuum wavelengths) and in units of Angstroms


As of now, if you name the spectrum file with the string "flattened" in it, the python script will automatically pick it up and add it to the parameters file.

**It is highly recommended** that you create a unique folder for each stellar spectrum, as MOOG++ will dump output files to the same directory the spectrum is in.

## 2. Choose Line Lists
In its default configuration, MOOG++ is built to fit lines of Fe, Ca, Ti, Mg, Si, C, O, Na, Al, and K, in that order, with a range of wavelengths chosen to coincide with KPF. These line lists are stored in the appropriately named `linelists` folder. This folder contains an additional file, `defaults.txt` which tells MOOG++ what files should be used to fit which elemental abundances.

If this all sounds good to you (you're interested in fitting the stated elements using KPF spectra), great! You can stop reading the rest of this section and go to the next step.

...

Otherwise, you'll have to pick what lines you want to fit. I have some helper scripts for this which will be available on written request (they're ugly and not ready to be public just yet, but they work). You're also welcome to simply pick out your favorite lines, though the wavelength data must match the wavelength data provided in linemake's repository of lines, otherwise you risk MOOG++ just not finding your lines. I recommend combing through linemake's `mooglists` directory and choosing the lines from there. This way, you can make sure that a) linemake actually has data for the lines you want, and b) the wavelengths you're telling MOOG++ the lines are located at are actually correct.


## 3. Run `./runStar.sh`
First, open up a terminal in the MOOG++ directory

The bash script requires a few command line arguments to kick everything off, but will carry the code from there to its completion. In order, they are as follows:

1. The name of your star, in a format resolvable by SIMBAD. This is used by the python script to get photometry and thus get stellar parameters
2. The folder which contains the spectrum file. This will be where all output files are dumped, and they are indeed dumped, without care for what else might have been in said folder, hence the recommendation to create a unique directory for each spectrum you wish to get abundances for.
3. An initial guess at the metallicity of your star (\[M/H\]). The solar value of 0.0 is always a safe place to start, but you might wish to modify it, as a better guess can reduce the number of iterations required to fit the stellar parameters, and can significantly affect execution time

Here is an example of what those command line arguments look like:

`./runStar.sh "HD 10700" /home/jared/Documents/Spectra/KPF/CAP4/10700/ -0.41`


# Understanding the Output
MOOG++ creates a directory for each element, named after its atomic number. In these folders, a text file is created for each line which was successfully fit. 

Generally, however, you will be most interested in the "params.txt" file. This file is used as both an input and an output file, and contains all the information you're likely to want out of MOOG++ at a glance:

- Stellar parameters:
  - $T_{\textrm{eff}}$, $\log{g}$, $[M/H]$, $[\alpha/M]$, $M/M_\odot$, $R/R_\odot$, and $L/L_\odot$. Each is listed with their values and 1-sigma uncertainties.

- Limb darkening coefficients:
  - Reported for a linear limb darkening model, with the name of the photometric band, wavelength in angstroms, and the coefficient itself.

- Elemental Abundances:
  - Labeled according to the atomic number of their element, and output in the $[X/H]$ format (abundance relative to solar, where $X$ is any given element e.g. Fe, C, etc.). The rightmost column is again 1-sigma error bars.

### Note:
If one of these error bars seems unusually large, you can try to go into the corresponding element's folder, open "abundanceSummary.txt" and remove lines from it that are significantly and obviously poor fits. This text file may be difficult to parse as of now, but contains best-fit data for each spectral line fit in the format 'wavelength abundance v_broad chi^2_nu' with one data point per line. Removing offensively bad lines from this file and then running 

`python computeParamFile.py "HD 10700" /home/jared/Documents/Spectra/KPF/CAP4/10700/ params.txt -0. 0. 4`

should improve things. 

**IMPORTANT: Before you do the above** As for runStar.sh, replace "HD 10700" with the name of your star, and the directory path with the path to your star's directory. Leave params.txt as is, but replace the two zeros in the above line with the \[M/H\] and \[alpha/M\] values respectively, from params.txt. 
###### I'm aware this is an extremely clunky way to perform manual sigma clipping of the abundances. I am certainly open to suggestions for improvement





