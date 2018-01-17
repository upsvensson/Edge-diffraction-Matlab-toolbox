# Edge diffraction Matlab toolbox

This is a Matlab toolbox for computing sound reflections and diffractions for external and internal scattering problems, in the time- or frequency-domain, for problems with Neumann boundary conditions. Some limited possibilities for Dirichlet problems exist as well. 
As of version 0.1, only external, convex Neumann scattering problems can be handled, in the frequency-domain.

## Getting Started

Copy all these m-files to a folder ("EDtoolbox" could be a good name for that folder). Add the path to this folder in Matlab, e.g., by giving the Matlab command

```
>>addpath <path to the folder where you installed the m files>
```
Get the package EDexamples, which has a few examples to get you started.

Inside the folder EDexamples, you find some example script files which execute various simple examples. In Matlab, you can, e.g., execute the file called EDexample_LspKessel_minimal.m, which should run the function EDmain_convexESIE and present a resulting frequency response in a plot window, and the geometrical model in another plot window.

### Prerequisites

You'll need Matlab.


## Running an example

Inside the folder EDexamples, you find some example script files which execute various simple examples.

## Authors

* **Peter Svensson** - *Initial work* - (https://github.com/upsvensson)


## License

This project is licensed under the GNU License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments



