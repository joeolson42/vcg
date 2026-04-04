# A vertical coordinate generator for WRF and MPAS.

## A slick tool to generate sigma/height coordinates for use in numerical weather prediction models. This tool allows for:

1. the specification of lowest model level height.
2. number of layers in a stable boundary layer and a linear
3. number of layers in the between the boundary layer top and the lower stratosphere
4. number of layers at very high levels (𝜎 < 0.1) for the gravity wave sponge layer and for satellite radiance data assimilation (not yet mature)
5. has the capability to allow for skewing the density of levels to lower or upper levels, if desired.
6. can output both sigma and height coordinates (for WRF and MPAS)

Depending which compiler you use, compile as:
* ifort -o vcg vcg.f90
* ifx -o vcg vcg.f90

etc, then simply execute the program, i.e., “./vcg”

But before you run the program, edit the namelist *vcg.nml*. For example:

* sig1       = 0.9977  !first sigma level above the surface (top of first model layer)
* nsig_cos   = 63      !number of layers in the analytical (cosine) sigma levels
* nsig_pbl   = 10      !number of layers within the stable pbl (from z=0 to z_spbl meters)
* z_spbl     = 300.    !height of stable boundary ayer (meters)
* nsig_p1    = 12      !minimum number of layers above 𝜎 = 0.1
* alfa1      = 1.0      !skewness factor (> 1 increases clustering near the surface; <1
* ptop       = 100.0    !model pressure top (Pascals), needed for the 𝜎-to-z transformation

With this flexibility, users have the potential to do both good and harm. Please refer to some of the literature for scientific guidance and use this tool wisely.

