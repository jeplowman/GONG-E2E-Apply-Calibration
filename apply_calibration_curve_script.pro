; Run this example script by typing @apply_calibration_curve_script.pro in the IDL command line
.compile apply_calibration_curve3.pro

; Change these to point to your location for the cal savefile, the magnetogram to be 
; calibrated (zqafile) and correspinging GONG intensity file (iqa), respectively:
cal_savefile = 'calibration_data_070119.sav'
zqafile = './uncal_fits/bbzqa100608t2004.fits.gz'
iqafile = './uncal_fits/bbiqa100608t2004.fits.gz'

; These will put the calibrated fits file and the example plot both in a subdirectory called cal_output
; of the current working directory:
out_dir='./cal_output/'
plot_dir='./cal_output/'

prange = 200 ; For output plot range of += 200 Gauss

; Apply the calibration. On return, cal will contain the calibrated magnetogram, uncal the uncalibrated one:
cal = apply_calibration_curve3(cal_savefile, zqafile, iqafile, out_dir=out_dir, plot_dir=plot_dir, prange=200, zqa_uncal=zqa_uncal)