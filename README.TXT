
This package contains model files described in:

	Barela et al. An Epilepsy Mutation in the Sodium Channel SCN1A 
	That Decreases 	Channel Excitability.  J. Neurosci. 26(10): 
	p. 2714-2723 

        Spampanato J, Aradi I, Soltesz I, Goldin AL (2004) Increased
        neuronal firing in computer simulations of sodium channel
        mutations that cause generalized epilepsy with febrile
        seizures plus. J Neurophysiol 91:2040-50

The Spamanato et al. 2004 paper used slightly different parameters in
the same model. However, we believe that the Barela et al. 2006 model
parameters as included in this collection of files is more accurate.

README.TXT 	    - this file
ichanWT2005.mod	    - Nav1.1 wild-type channel model	
ichanR859C1.mod     - R859C mutation (w/ B1)
APthreshold.hoc	    - generates data to find action potential thresholds
                      

To generate action potential data using APthreshold.hoc:

1. Install latest pre-compiled package of NEURON (if not already installed)
   and allow it to associate .HOC files.

2. Copy these files to a working directory (eg. C:\R859C)

3. Run the NEURON mknrndll utility to compile the model files.  
   mknrndll is likely found under    Start / Programs / NEURON 
   Use the mknrndll directory selector to navigate to the folder
   where you copied these files and while in that directory choose 
   "Make nrnmech.dll".  A shell window will open and you should 
   get a message saying nrnmech.dll was created.  
   
4. Double click APthreshold.hoc to execute the program using the 
   default parameters (wild-type data).  The data will be saved
   in a new file called APthreshold.txt.  (Any existing file with 
   this name will be overwritten).
   
5. To generate similar data for the R859C channels change the value 
   of the useCellType variable set near the top of the APthreshold.hoc 
   file to "R859C".  

   

   
   
   


  

