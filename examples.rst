Examples
============

In Python, run:

::


  import Mshpy


Then you can run Mshpy.main as explained in the documentation.

::


  Mshpy.main(path)

* path: File directory. The directory should includes

1. SW/IMF conditions file
2. sc trace file
For 1 and 2, see 'Datasets' for the detail.

3. MP/BS offset file (Optional)
File containing following information
text format:
mpoff bsoff
mpoff: Manual magnetopause offset along x axis.
    Magnetopause model for the MHD-based MSH model is Shue et al. 1998,
    except for the Romashets et al. [2019] model which uses Jellinek et al. 2012 model.
bsoff: Manual bow shock offset along x axis
    Bow shock model for the MHD-based MSH model is Jelinek et al. 2012.
