# JWST-destripe-driver

This repo is a simple driver (with a few utils) to run image destriping via pyimcom.imdestripe on JWST images.

# Get code

`git clone git@github.com:Roman-HLIS-Cosmology-PIT/jwst-destripe-driver.git`

# Install code

`pip install .`

If you are developing locally, editable install is recommended:

`pip install -e .`

# Run

`python run.py <config>`


# Notes

- Set `INSTRUMENT=NIRCAM` in your environment before running so pyimcom uses JWST settings.
- `DSOBSFILE` in config should be the full path prefix up to and including `jw` (for example: `/path/to/data/jw`).
- Input files are discovered with pattern `jw*_nrcb*_crf.fits` and parsed like
	`jw04793001001_02101_00001_nrcb1_crf.fits`.

This repo is developed and maintained by Katherine Laliotis and Devisree Tallapaneni.
