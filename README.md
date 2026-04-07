# A vertical coordinate generator for WRF and MPAS.

## A slick tool to generate sigma/height coordinates for use in numerical weather prediction models. This tool allows for:

1. the specification of lowest model level height.
2. number of layers in a stable boundary layer (specified with a linear increase in dz with height)
3. approximate number of layers in the between the boundary layer top and the lower stratosphere
4. number of layers at very high levels (𝜎 < 0.1) for the gravity wave sponge layer and for satellite radiance data assimilation (not yet mature--it works but the type of blending needs to be improved. For now, I suggest making *nsig_cos* large if you want a lot of layers in the stratosphere.)
5. has the capability to allow for skewing the density of levels to lower or upper levels, if desired.
6. can output both sigma and height coordinates (for WRF and MPAS, respectively)

Depending on which compiler you use, compile as:
* ifort -o vcg vcg.f90
* ifx -o vcg vcg.f90

etc, then simply execute the program, i.e., “./vcg”

But before you run the program, edit the namelist *vcg.nml*. For example:

* sig1       = 0.9977  ${\color{red}!first\space sigma\space level\space above\space the\space surface\space (top\space of\space first\space model\space layer)}$
* nsig_cos   = 66      ${\color{red}!number\space of\space layers\space in\space the\space analytical\space (cosine)\space sigma\space levels}$
* nsig_pbl   = 7        ${\color{red}!number\space of\space layers\space within\space the\space stable\space pbl\space (from\space z=0\space to\space zpbl,\space  meters)}$
* z_spbl     = 210.    ${\color{red}!height\space of\space stable\space boundary\space layer\space (meters)}$
* nsig_p1    = 9        ${\color{red}!minimum\space number\space of\space layers\space above\space 𝜎 = 0.1}$
* alfa1      = 1.0      ${\color{red}!skewness\space factor\space (> 1\space increases\space clustering\space near\space the\space surface; <1\space increases\space clustering\space aloft}$
* ptop       = 200.0    ${\color{red}!model\space pressure\space top\space (Pascals),\space needed\space for\space the\space 𝜎-to-z\space transformation}$
* ddelz_max  = 300.     ${\color{red}!maximum\space allowable\space increase\space in\space delta-z\space (meters)}$

I have not (and can not) fully tested all possible configurations, so I can not guarantee you won't find some irregularities with some configurations. I recommened plotting the ∆𝜎 and ∆(∆𝜎) [or ∆z and ∆(∆z)] to check for irregulaties before you use the profiles in a numerical weather prediction model.

With this flexibility, users have the potential to do both good and harm. Please refer to some of the literature for scientific guidance and use this tool wisely.

