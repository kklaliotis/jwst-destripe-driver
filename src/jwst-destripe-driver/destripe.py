import copy
import glob
import json
import os
import sys

import numpy as np
from astropy.io import fits
from pyimcom import imdestripe
from pyimcom.config import Settings, JWST


def destripe(cfg_file, noiseid=None, verbose=False, max_files=None):
    """
    Destripes one layer from the indicated set of files.

    The data in the files are *overwritten* and the ``processinfo`` leaf gets a new
    field: ``file["processinfo"]["destripe"]`` gives the number of
    noise layers that have been destriped, and ``file["processinfo"]["destripe_complete"]``
    is set to True if everything has been destriped.

    Parameters
    ----------
    cfg_file : str
        The configuration file.
    noiseid : int or None, optional
        If specified, destripes that particular noise layer.
        Otherwise, does the science layer.
    verbose : bool, optional
        Whether to talk a lot to the output.

    Returns
    -------
    int or None
        Number of noise layers. None if no files found.

    """

    # first get the file prefix and information
    with open(cfg_file, "r") as file:
        cfg = json.load(file)

    if JWST: Settings.jwst()
    filter = Settings.RomanFilters[cfg["FILTER"]]

    # add tails as needed
    file_prefix = cfg["DSOBSFILE"]
    if verbose:
        print("File prefix =", file_prefix)


    # (file, identifier)
    fdir, fileprefix = os.path.split(file_prefix)
    n = len(fileprefix)
    numus = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "_"}

    use_files = []
    for f in os.listdir(fdir):
        if len(use_files) == max_files:
            break
        # What does this do? It checks if the file starts with the prefix, ends with _crf.fits, and has only numbers and underscores in the middle part. If so, it adds it to the list of files to use.
        if f[:n] == fileprefix and f[-9:] == "_crf.fits" and all(c in numus for c in f[n:-9]):
            use_files.append((os.path.join(fdir, f), f[n:-9])) # adds the full path and the identifier (the part between the prefix and _crf.fits)

    if verbose:
        print("Files selected:", use_files)

    # if there aren't any files, return None
    if len(use_files) == 0:
        print("No files found; aborting.")
        return None

    # cleanup output directory (except for overlap matrices)
    clearfiles = glob.glob(os.path.join(cfg["DSOUT"][0] + "/masks", "*_mask.fits"))
    clearfiles += glob.glob(os.path.join(cfg["DSOUT"][0] + "/masks", "*_mask.fits.lock"))
    clearfiles += glob.glob(os.path.join(cfg["DSOUT"][0], "*.fits"))
    if noiseid is None:
        # files to clear only the first time; these are reused for each noise layer
        clearfiles += glob.glob(os.path.join(cfg["DSOUT"][0], "ovmat.npy"))
        clearfiles += glob.glob(os.path.join(cfg["DSOUT"][0], "SCA_list.txt"))
        clearfiles += glob.glob(os.path.join(cfg["DSOUT"][0], "*.out"))
    if verbose:
        print("Clearing files:")
        for f in clearfiles:
            print("    ", f)
        print("")

    for p in clearfiles:
        os.remove(p)

    # main destriping
    dsout = imdestripe.main(cfg_file, overlaponly=False)
    if verbose:
        print("Output -->", dsout)
        print("")
        sys.stdout.flush()

