Datasets
=============

MHD-based model data files (*.npy) can be found in ``/modeling/cn`` and
``/modeling/cs``.

SW/IMF conditions file and spacecraft (sc) trace file are needed.

1. SW/IMF conditions file

* text format: (without header)
* Year DOY HR MN Bx By Bz Vx Vy Vz n Pd Ma Mm, as in OMNIweb
* time is for the matching with sc trace.
* One line data can be used for constant SW cond.
* SW file should be named omni.lst

2. sc trace file

* text format: (without header)
* yy/mm/dd x y z, as in SSCWeb (x,y,z in GSE)
* Notice that you need to choose yy/mm/dd output option instead of yyyy ddd.
* sc trace file should be named orbit.txt

3. MP/BS offset file (Optional) named 'Msh_offset.txt'

*   text format:
*   mpoff bsoff

1.mpoff: Manual magnetopause offset along x axis.
      Magnetopause model for the MHD-based MSH model is Shue et al. 1998,
      except for the Romashets et al. [2019] model which uses Jellinek et al. 2012 model.
   
2.bsoff: Manual bow shock offset along x axis
      Bow shock model for the MHD-based MSH model is Jelinek et al. 2012.
