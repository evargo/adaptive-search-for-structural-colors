# adaptive-search-for-structural-colors #

## Package Contents ##
READ_ME: Mathematica documentation for the package<br>
StructuralColorFramework: general functions used by every adaptive mesh search<br>
StructuralColor1D: designing for a 1D stack<br>
StructuralColor2D: designing for a 2D lattice<br>
StructuralColor3D: designing for a 3D semi-ordered glass<br>
StructuralColorTemplate: template to define a custom search over a different geometry<br>
Design_From_An_Image: example of choosing geometries needed to create a color image<br>
ONYX_Modified_Fortran_Code: Fortran files used for the 2D calculation

## Documentation ##

### General Start-Up ###
To use the package, open a new Mathematica notebook in the same location as the package, and evaluate these lines:<br>
`Get["StructuralColorFramework", Path -> NotebookDirectory[]]`<br>
This loads the general Structural Color package.<br>

### Colors in Mathematica ###
Colors can be visualized in Mathematica with the XYZColor function, which can be used to visually choose colors.<br>
`XYZColor[0.27, 0.26, 0.74]`<br>
![Sample Mathematica Color Object](https://github.com/evargo/adaptive-search-for-structural-colors/blob/master/color_object.png)<br>

The Color object can be converted back into a set of XYZ coordinates with List@@

### Design of 1D Photonic Stacks ###
Make sure you first load the Framework package, then run the following code:<br>
`Get["StructuralColor1D", Path -> NotebookDirectory[]]`<br>

This is the simplest of the packages, but they're all organized similarly.<br>

If you simply wish to know the spectrum or color of a certain photonic stack, use the Spectrum1D function:
`Spectrum1D[l1, l2, {n1, n2, layerNum, theta}]`, where the arguments are: width 1 (nm), width 2 (nm), refractive index 1, refractive index 2, number of layers, viewing angle (degrees). `XYZColor[SpectrumToColor[Spectrum1D[300, 200, {1.49, 1, 5, 0}]]]` visualizes the color.

This code was originally published via the Wolfram Demonstrations Project, http://demonstrations.wolfram.com/MultilayerPhotonicBandgap

To search for parameters given a desired color, use the `PhotonicSystem1D` function to hold fixed parameters. The default is searching over layer thicknesses, but this can be changed in the package. A 5-layer system with of PMMA and air, viewed from above, has the following entry: `PhotonicSystem1D[1.49, 1, 5, 0]`.

Now you're ready to iterate! Stick your `PhotonicSystem` into the function `Iteration1D`, along with the color you're searching for:  `Iteration1D[{0.2, 0.3, 0.4}, PhotonicSystem1D[1.49, 1, 5, 0]]`

Your output will look like this:<br>
![1D Output](https://github.com/evargo/adaptive-search-for-structural-colors/blob/master/1D_output.png)<br>

If you want to keep iterating, evaluate the same code again. Once you're ready to switch colors, evaluate the function `IterationReset1D`.

### Design of 2D Photonic Lattices: Requires Fortran ###

`Get["StructuralColor2D`", Path -> NotebookDirectory[]]`

This calculation is prohibitively slow in Mathematica, so the package integrates with a Fortran compiler. 

The Fortran code was written by Ward and Pendry, and the original documentation is found here: http://www.cmth.ph.ic.ac.uk/photonics/Newphotonics/pdf/CompPCom128_590.pdf
The version of ONYX you downloaded has been modified, but the underlying calculation is the same.

The README Mathematica Notebook contains detailed instructions for setting up the Fortran code.

The `PhotonicsSystem2D`function again holds fixed parameters. The default is searching over hole radii and lattice spacings, but this can be changed in the package. A system made of air holes in silicon has the following entry: `PhotonicSystem2D[1, 3.4]`.

Now you're ready to iterate! Stick your `PhotonicSystem` into the function `Iteration2D`, along with the color you're searching for:  `Iteration2D[{0.2, 0.3, 0.4}, PhotonicSystem2D[1,3.4]]`Your output will look like this:<br>
![2D Output](https://github.com/evargo/adaptive-search-for-structural-colors/blob/master/2D_output.png)<br>

### Design of 3D Photonic Glass ###

`Get["StructuralColor3D", Path -> NotebookDirectory[]]`

The 3D code runs exactly like the 1D code, so no sweat!

The code for this calculation was originally written in Python by the Manoharan Group. Their code can be found here: https://github.com/manoharan-lab/structural-color
The translated Mathematica code is found as part of the StructuralColor3D package.

The relevant `PhotonicSystem` here is `PhotonicSystemGlass`. The default is searching over sphere radii and packing fraction, but this can be changed in the package. A system made of PMMA spheres in water has the following entry: `PhotonicSystemGlass[1.49, 1]`.

Again, plug your `PhotonicSystem` into an iterating function: `Iteration3D[{0.2, 0.4, 0.3}, PhotonicSystemGlass[1.49, 1.]]`. To continue iterating, evaluate the function again. To switch to a new color, use `IterationReset3D`.

Your output will look like this:<br>
![3D Output](https://github.com/evargo/adaptive-search-for-structural-colors/blob/master/3D_output.png)<br>

### Design of Your Own Geometry ###

If you'd like to design for a completely new geometry, all you need is code that will output a reflectance spectrum given a set of parameters. It's easiest if this code is written in Mathematica, and it's often worth your time to translate into Mathematica from another language if the code is short enough. If you don't want to translate, look at the 2D sub-package for an example of interfacing with an external compiler.

The `StructuralColorTemplate` sub-package is meant to guide you through the process of making your own sub-package. (It may be useful to open `StructuralColor1D` as well as a reference.)

### Using an Image as a Template ###

The notebook `Design_From_An_Image` will walk you through the process. An example is shown below:
![Design from an Image](https://github.com/evargo/adaptive-search-for-structural-colors/blob/master/Design_from_an_image.png)<br>
