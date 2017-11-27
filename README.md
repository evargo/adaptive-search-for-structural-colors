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
To use the package, open a new Mathematica notebook in the same location as the package, and evaluate these lines:
`Get["StructuralColorFramework", Path -> NotebookDirectory[]]`
This loads the general Structural Color package.

### Colors in Mathematica ###
Colors can be visualized in Mathematica with the XYZColor function, which can be used to visually choose colors.
`XYZColor[0.27, 0.26, 0.74]`

The Color object can be converted back into a set of XYZ coordinates with List@@
