# Edge diffraction Matlab toolbox

This is a Matlab toolbox for computing sound reflections and diffractions for external and internal scattering problems, in the time- or frequency-domain, for problems with Neumann boundary conditions. Some limited possibilities for Dirichlet problems exist as well. 
As of version 0.1, only external, convex scattering problems can be handled.

## Getting Started

Install all these m-files in a folder. Add the path to this folder in Matlab, e.g., by giving the Matlab command

```
>>addpath <path to the folder where you installed the m files>
```

Now you could execute the main program/function EDmain_convexESIE from the command window - but you will get an error message since that function expects some input data. See Examples below.

Inside the folder EDexamples, you find some example script files which execute various simple examples. In Matlab, you can, e.g., execute the file called EDexample_LspKessel_minimal.m, which should run the function EDmain_convexESIE and present a resulting frequency response on the screen.

### Prerequisites

You'll need Matlab.

```
Give examples
```

### Installing

A step by step series of examples that tell you have to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone who's code was used
* Inspiration
* etc

