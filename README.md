# Edge diffraction Matlab toolbox

This is a Matlab toolbox for computing sound reflections and diffractions for external and internal scattering problems, in the time- or frequency-domain, for problems with Neumann boundary conditions. Some limited possibilities for Dirichlet problems exist as well. 
As of version 0.1, only external, convex Neumann scattering problems can be handled, in the frequency-domain.

## Getting Started

Copy all these m-files to a folder ("EDtoolbox" could be a good name for that folder). Add the path to this folder in Matlab, e.g., by giving the Matlab command

```
>>addpath <path to the folder where you installed the m files>
```

Get two files from Mathworks (Matlab File Exchange): lgwt.m and DataHash.m. Store them somewhere, where Matlab finds them.

Get the package EDexamples, which has a few examples to get you started.

Inside the folder EDexamples, you find some example script files which execute various simple examples. In Matlab, you can, e.g., execute the file called EDexample_LspKessel_minimal.m, which should run the function EDmain_convexESIE and present a resulting frequency response in a plot window, and the geometrical model in another plot window.

## Prerequisites

You'll need Matlab, and those two files mentioned above, from Matlab File Exchange.


## Running an example

Inside the folder EDexamples, you find some example script files which execute various simple examples.

## Documentation

The file EDtoolbox_manual.pdf gives some details about the EDtoolbox. Also, the m-file EDversionhistory displays, in the Matlab command window, some documentation on changes in the toolbox.

## Version numbering

A new version number is given whenever some change is made that affects the numerical results, for all or some cases. If there are just minor changes, e.g., to what is printed out on the screen, or some internal code improvement, then there will just be a date/time change. 

## Authors

* **Peter Svensson** - *Initial work* - (https://github.com/upsvensson)


## License

This project is licensed under the BSD License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

Many people have contributed on a smaller and larger scale during more than twenty years of work. The last year, the collaboration with Sara Martin, Jan Slechta and Jason Summers is acknowledged in particular.
Parts of this work have been financed through Sara Martin’s project ”Hybrid methods” funded by the Research Council of Norway, and through Jan Slechta’s stipend from the ERCIM Alain Bensoussan Fellowship Programme.


