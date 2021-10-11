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
