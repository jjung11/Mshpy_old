Usage
-----------------------------------------

Main module is Msh_Nstep_master.py

Args:
    1. path: File directory. For input files, see ``Datasets``
    2. mpoff: Optional; Manual magnetopause offset along x axis.
        Magnetopause model for the MHD-based MSH model is Shue et al. 1998.
    3. bsoff: Optional; Manual bow shock offset along x axis
        Bow shock model for the MHD-based MSH model is Jelinek et al. 2012

Returns:
    Following files are saved in the path.
        1. Msh_MHD_out.txt: MHD-based model result.
        2. Msh_Soucek_out.txt: Soucek & Escoubet [2012] model plasma velocity results.
        3. Msh_Romashtes_out.txt: Romashets et al. [2019] model magnetic field results.
