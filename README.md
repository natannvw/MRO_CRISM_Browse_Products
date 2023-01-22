# MRO_CRISM_Browse_Products

`browse_product()` can load Browse Products[^1] of pre-calculated Summary Parameters[^1] from SU or SR cubes, or calculate them on its own. Or, it can be used also to calculate a mineral lab reference color that expected to appear in the same browse product as did Viviano et al. ([2014]([^1])) in the following figure:


Usage:
1. Get scene files:
For example [FRT00009326](https://ode.rsl.wustl.edu/mars/productsearch) scene.

2. Run:
````matlab
browse_product()
````
Just follow the GUI instructions. I recommend using the default settings for the first time by hitting ‘Enter’ (and clicking on ‘OK’ in stretch parameters choose) through the process.

[^1]: Viviano et al. (2014)  Revised CRISM spectral parameters and summary products based on the currently detected mineral diversity on Mars (https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2014JE004627)
