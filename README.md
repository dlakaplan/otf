# otf
On-the-fly pulsar mapping


## Requirements:
* `astropy`, `pyyaml`
* [`presto`](https://github.com/scottransom/presto) is desired (and needed for most things) but it won't cause `import` to fail if it's not available.

## Example:

* Extract the pointing info from PSRFITS files: 
  * `extract_pointing vegas_59328_55419_J2035+36_0001_0001.fits` generates `vegas_59328_55419_J2035+36_0001_0001.pointing.ecsv`
  * `extract_pointing vegas_59328_56190_J2035+36_0002_0001.fits` generates `vegas_59328_56190_J2035+36_0002_0001.pointing.ecsv`
* Then get the position from the `pfd` files:
  * `fit_otf --pfd vegas_59328_55419_J2035+36_0001_0001_ACCEL_Cand_13.pfd vegas_59328_56190_J2035+36_0002_0001_ACCEL_Cand_2.pfd --pointing vegas_59328_55419_J2035+36_0001_0001.pointing.ecsv vegas_59328_56190_J2035+36_0002_0001.pointing.ecsv --plot=test.png`
  * Should return `20h35m23.7488s +36d56m12.1434s` and make the plot below (see `example/` directory for files):
   [Output plot](examples/test.png)
