# concreteTSD
These scripts and functions are the companion code to my Ph. D. Dissertation on TSD data analysis and interpretation for jointed pavements.

# Link to Dissertation
My Dissertation can be retrieved from: https://vtechworks.lib.vt.edu/handle/10919/111606 

# Bsaic use instructions
Matlab (or GNU-Octave [see below for compatibility concerns]) is needed to execute this code! None of these software tools is provided here. Matlab(R) can be acquired from The Mathworks Inc., GNU-Octave can be downloaded from https://www.gnu.org/software/octave/index
You may only need to execute the provided "front-end" scripts to analyze the TSD data. Both front-ends provided here are ready to work with the TSD data chunks in the TSD data folder. If you're using another set of TSD measurements, modify the front-end scripts accordingly before executing.
Observe as well that the processed outcome is exported as Excel spreadsheets. Bear that in mind when exectuing the front-ends that you are sure you may not be overriding a previous results file!.
## A comment on excecution time: 
The BPD+RWWL1 example case study using the National Mall data may take 45-60 minutes to complete [analyzing roughly 80000 data points from 57 road semgents]. The back-calculation code with one MnROad Loop example may take around 1/2 hour to complete. These tests were done on a personal computer with a 3rd-Gen i5 CPU running Matlab on Windows (single-threadedly), these front-ends running on a newer PC and the same 'environment' may take 2/3 that time.

# Dependencies
1) This code needs the Wavelab package to perform Basis Pursuit denoising and RWL1. You can retrieve Wavelab from: https://statweb.stanford.edu/~wavelab/Wavelab_850/download.html
2) The BPD+RWL1 front-end script is capable of doing batch-analysis of several segments within a single TSD data file. If you'd like the resulting plots to be exported to PDF, the code can do, but before exectuing you need to retrieve a copy of the EXPORTFIG package. EXPORTFIG is available through:  https://www.mathworks.com/matlabcentral/fileexchange/23629-export_fig/

# License Terms and Disclaimer 
This code is distributed under the terms of the Creative Commons BY-SA 4.0 International. The terms of the license can be found in: https://creativecommons.org/licenses/by-sa/4.0/
Also, this code is distributed AS-IS. If you agree to utilize this code, you become responsible of whatever decisions you make out of the code's outcome. Yet anyway, if you encounter a big bug in the code or find something that may warrant a revision, please feel welcome to reach out so that the issue can be tackled (I can be found over linkedIn).
Please respect the license terms, and may this code be useful to you. Enjoy!!.

# GNU-Octave compatibility
This code might work on GNU-Octave; as of May 2022, I haven't tested it myself. However, since it both reads and writes to MS Excel files, it may be very slow to initialize and export the results.
