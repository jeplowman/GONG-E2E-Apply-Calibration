;+
; function: apply_calibration_curve3
;
; Apply a calibration curve to a GONG magnetogram file, optionally with Monte Carlo. Inputs:
;
;	cal_savfile: Save file containing the data used by the calibration.
;	zqafile: GONG magnetogram fits file (these files often contain 'zqa' in their name)
;	iqafile: GONG pseudocontinuum intensity file (often contains iqa in its name). Used
;			to determine AR and non-AR boundaries -- regions (size: dsmall) that are dark
;			compared to the local (size:dlarge) median intensity are considered active regions.
;	dlarge: Large median smoothing diameter for determining AR boundaries. Default: 30
;	dsmall: Small median smoothing diameter for determining AR boundaries. Default: 3
;	out_dir: Directory to place calibrated magnetograms.
;	plot_dir: Directory to place comparison figures.
;	prange: Magnetogram intensity range for comparison figures.
;	ondisk_frac: Fraction of solar radius considered to be 'on disk'. Default: 0.99
;	zqa_uncal: Optionally, output uncalibrated magnetogram as well.
;	do_mc: Do Monte Carlo (see below).
;	mc_max_sunspot_angle: Maximum angle to include in MC for sunspot calibration (default: 0.0)
;
; Monte Carlo: For a given measured value to be calibrated, the scatter in the `ground truth' 
; for all points in the calibration data set at that value is referred to and a random scatter
; is chosen from them.
;
; Note that I personally do not recommend using the Monte Carlo for the calibration. On large scales, 
; it does the same thing as an averaged curve but on a per-pixel level it adds scatter equal to that 
; of all values in the calibration data set near the given measured value. Because it includes all 
; such points, not just the ones that are in a similar context to the point in question, and it does 
; so in a way that is incoherent from pixel to pixel, this exaggerates the pixel to pixel variation 
; in the resulting calibrated data set in a way that is not consistent with the input magnetogram, or
; (for that matter) the actual variation in the calibration data set: application of this Monte Carlo
; method will result in a larger pixel-to-pixel variation than that of the ground truth in the
; calibration data set. I have included it as an option on the insistence of other members of this
; project.
;
; The calibration save file results from the instrument simulation and its comparison with the 
; MURaM ground truth. It has the following components:
;
;	cal_measimages: The `measurement' images produced by the GONG simulator, one for each angle.
;	cal_gtimages: The `ground truth' images (one for each angle) corresponding to cal_measimages.
;	qsflags: Pixels flagged as QS/non-sunspot in the process of making calibration curves
;			(note: only reliable for tile_levels=1, i.e., the GONG resolution tile level).
;			One image for each angle.
;	arflags: Like qsflags, but for sunspots. See paper.
;	tile_levels: Images giving the tiling levels in the simulation -- 1=GONG resolution, 
;			2=Twice GONG resolution, etc. One image for each angle.
;	ondisks: Images giving which pixels are ondisk, one image for each angle.
;	gongitots: GONG pseudo-continuum intensity images produced by the GONG simulator, one for each angle.
;	angles: The angles at which the radiative transfer and GONG simulator were run. E.g., 
;			0, 25, 45, 60, 70, and 75 degrees. In radians.
;	angles_deg: Same as angles but in degrees.
;	anglestrs: Same as angles_deg, but as strings.
;	qs_obsfits: The measurement (obs) axis of the calibration curve fit for QS/non-sunspot, one set of
;			values for each angle.
;	qs_modfits: The `ground truth' (mod) axis of the calibration curve fit for QS/non-sunspot, one
;			set of values for each angle. Calibration can be applied for a single angle (index i) by
;			interpolating: e.g., mag_cal = interpol(qs_modfits[i,*],qs_obsfits[i,*],mag).
;	ar_obsfits, ar_modfits: Same as qs_obsfits and qs_modfits but for sunspot/active regions.
;
;	Note that qs_obsfits, qs_modfits and ar_obsfits, ar_modfits are fit at low resolution but
;	reinterpolated to high resolution so that they can all be stored in the same array. The curves
;	are explicitly constructed with the expectation that they will be applied by linear interpolation
;	(see paper) and this reinterpolation has negligible effect as long as the new resolution is much
;	higher than the original one and linear interpolation is used throughout.
;
;	07-01-19 -- JEP
;
;-

;+
; This function applies median smoothing to an image. Unlike the built-in IDL median function
; it uses a circular, not square, window, which avoids some of the artifacts that a
; retangular window can produce (because it is anisotropic). Edge treatment is equivalent
; to IDL's edge_truncate for convolution/smoothing. Inputs:
;
;	datain: The image to be median smoothed.
;	widthin: The width of the median smoothing window. This is forced to be an odd integer, which ensures
;			that the window is centered. 
;
; Returns a median smoothed image.
;
;-
function median2, datain, widthin
		
	nx = n_elements(datain[*,0])
	ny = n_elements(datain[0,*])
	
	pad = floor(widthin*0.5)
	width = pad*2+1
	
	data = dblarr(nx+2*pad,ny+2*pad)
	data[pad:-pad-1,pad:-pad-1]=datain

	for i=0,pad-1 do begin
		data[i,*] = data[pad,*]
		data[nx+2*pad-i-1,*] = data[-pad-1,*]
	endfor
	for i=0,pad-1 do begin &$
		data[*,i] = data[*,pad]
		data[*,ny+2*pad-i-1] = data[*,-pad-1]
	endfor
	
	xa = (lindgen(width))#transpose(1+lonarr(width))-pad
	ya = (1+lonarr(width))#transpose(lindgen(width))-pad
	rada = sqrt(xa*xa+ya*ya)
	xa_in = xa[where(floor(rada) le pad)]+pad
	ya_in = ya[where(floor(rada) le pad)]+pad
	dataout = dblarr(nx,ny)
	for i=0,nx-1 do begin
		for j=0,ny-1 do begin
			dataout[i,j] = median(data[xa_in+i,ya_in+j])
		endfor
	endfor

	return,dataout
	
end


;+
; This subroutine is the inner section of the loop for applying the calibration curves with Monte Carlo:
;-
function do_mc_cal_inner, img, imgflags, cal_meas, cal_gt, calflags

	imgvals = img[imgflags]
	n_imgvals = n_elements(imgvals)
	imgvals += 0.1*randomn(seed,n_imgvals)
	cal_measvals = [-cal_meas[calflags],cal_meas[calflags]]
	cal_gtvals = [-cal_gt[calflags],cal_gt[calflags]]
	n_calvals = n_elements(cal_measvals)
	msort = sort(cal_measvals)
	indices = floor(interpol(lindgen(n_calvals),cal_measvals[msort],imgvals)+randomu(seed,n_imgvals))
	meas_cal = cal_gtvals[msort[indices]]

	print,n_elements(cal_meas),n_elements(cal_gt),n_elements(imgvals),min(meas_cal),max(meas_cal),min(cal_gtvals),max(cal_gtvals)

	nbins = 40
	binbs = floor((n_calvals-1)*findgen(nbins+1)/nbins)
	binlocs = 0.5*(binbs[0:nbins-1]+binbs[1:nbins])
	deltas = fltarr(nbins)
	for i=0,nbins-1 do deltas[i] = 0.1+ 1.0*(stdev(cal_gtvals[msort[binbs[i]:binbs[i+1]]])^2+stdev(cal_measvals[msort[binbs[i]:binbs[i+1]]])^2)/sqrt(n_calvals)
	
	return, meas_cal+(interpol(deltas,binlocs,imgvals) < 1.0)*randomn(seed,n_imgvals)

end

;+
; Outer loop of the MC calibration subroutine of apply_calibration_curve3. Not meant to be called
; independently of apply_calibration_curve3, see comments on that subroutine for additional
; description/caveats.
;-
function do_mc_cal, image, qsflags, arflags, angles, cal_measimages, cal_gtimages, cal_qsflags, cal_arflags, cal_angles, qs_obsfits, qs_modfits, ar_obsfits, ar_modfits, angles_mc

	n_angles = n_elements(cal_angles)
	nx = n_elements(image[*,0])
	ny = n_elements(image[0,*])

	img_cal = fltarr(nx,ny)
	ang_dithimg = randomu(seed,nx,ny)
	angle_indices0 = interpol(indgen(n_angles),cal_angles,angles_mc) < n_angles
	angle_indices = floor(ang_dithimg + angle_indices0) < ceil(angle_indices0)

	qs_px = where(qsflags)
	ar_px = where(arflags)

	zqa_cals = fltarr(nx,ny,n_angles)
	for i=0,n_angles-1 do begin
		zqa_cal = fltarr(nx,ny)	
		zqa_cal[qs_px] = interpol(qs_modfits[i,*],qs_obsfits[i,*],image[where(qsflags)])		
		zqa_cal[ar_px] = interpol(ar_modfits[i,*],ar_obsfits[i,*],image[where(arflags)])
		zqa_cals[*,*,i] = zqa_cal
	endfor

	for i=0,nx-1 do begin
		for j=0,ny-1 do begin
			img_cal[i,j] = interpol(zqa_cals[i,j,*],cal_angles,angles[i,j])
		endfor
	endfor
	
	for i=0,n_angles-1 do begin
		calimg_mc = fltarr(nx,ny)
		resid_img = fltarr(nx,ny)
		cal_measimg = cal_measimages[*,*,i]
		cal_gtimg = cal_gtimages[*,*,i]
		cal_wqs = where(cal_qsflags[*,*,i])
		cal_war = where(cal_arflags[*,*,i])
		resid_img[cal_wqs] = cal_gtimg[cal_wqs]-interpol(qs_modfits[i,*],qs_obsfits[i,*],cal_measimg[cal_wqs])
		resid_img[cal_war] = cal_gtimg[cal_war]-interpol(ar_modfits[i,*],ar_obsfits[i,*],cal_measimg[cal_war])
		calimg_mc[where(qsflags)] = do_mc_cal_inner(image, where(qsflags), cal_measimg, resid_img, cal_wqs)
		calimg_mc[where(arflags)] = do_mc_cal_inner(image, where(arflags), cal_measimg, resid_img, cal_war)
		img_cal[where(angle_indices eq i)]+=calimg_mc[where(angle_indices eq i)]
		print,n_elements(where(qsflags)),n_elements(where(arflags)),n_elements(where(angle_indices eq i)),min(calimg_mc),max(calimg_mc)
	endfor

	return, img_cal

end

function apply_calibration_curve3, cal_savefile, zqafile, iqafile, dlarge=dlarge, dsmall=dsmall, $
		out_dir=out_dir, plot_dir=plot_dir, prange=prange, ondisk_frac=ondisk_frac, zqa_uncal=zqa_uncal, $
		do_mc=do_mc, mc_max_sunspot_angle=mc_max_sunspot_angle
	
	restore,cal_savefile

	if(n_elements(ondisk_frac) eq 0) then ondisk_frac=0.99

	if(n_elements(dlarge) eq 0) then dlarge = 30
	if(n_elements(dsmall) eq 0) then dsmall = 3

	; Don't apply QS fit for fluxes larger than what are in the calibration curve array:
	if(n_elements(bmax_qs) eq 0) then bmax_qs = max(min(abs(qs_obsfits),dimension=1))

	; Read the GONG fits files:
	iqa = readfits(iqafile,iqahdr)
	zqa = readfits(zqafile,zqahdr)
	zqa_uncal = zqa
	zqa0 = zqa

	; Replace zeros with a very small placeholder value for plotting purposes
	; (they will be restored to zero at the end).
	zerovals = where(zqa0 eq 0)
	zqa[zerovals] = 1.0e-10
	
	; Large diameter median smoothing to estimate local pseudocontinuum intensity:
	iqa_ms_large = median2(iqa,dlarge)

	; Small diameter median smoothing:
	iqa_ms_small = median2(iqa,dsmall)

	rsun_asec = sxpar(iqahdr,'RADIUS')*3600.0D0*180.0D0/!PI
	
	nx = n_elements(iqa[*,0])
	ny = n_elements(iqa[0,*])

	xa = lindgen(nx)#(1+lonarr(ny))
	ya = (1+lonarr(nx))#lindgen(ny)

	; To find the image center, we begin by assuming that the offdisk section
	; of the image is contiguous, less than half the maximum
	; of the median smoothed pseudo-continuum intensity, and includes the lower
	; left corner of the image (index [2,2]). The edges of the image are also
	; explicitly included.
	regions = label_region(iqa lt 0.5*max(iqa_ms_large),/all_neighbors,/ulong)
	offdisk = regions eq regions[2,2] or (xa eq 0) or (xa eq nx-1) or (ya eq 0) or (ya eq ny-1)
	ondisk = offdisk eq 0

	; Compute the image center based on the above criterion:
	xc = mean(xa[where(ondisk)])
	yc = mean(ya[where(ondisk)])
	dx = 2.5*1.021;sxpar(iqahdr,'CDELT1')
	dy = 2.5*1.021;sxpar(iqahdr,'CDELT2')

	; To refine the offdisk region, we use the solar radius value in the header:
	ra = sqrt((xa-xc)^2+(ya-yc)^2)
	rsun = rsun_asec/sqrt(dx*dy)
	ondisk = ra lt rsun*ondisk_frac
	offdisk = ra gt rsun*ondisk_frac
	
	; AR/Sunspot is where small diameter median smoothed image is less than 98% of large diameter 
	; median smoothed image, or the measured field values is greater than bmax_qs and (in either
	; case) its radius is less than 0.95 rsun.
	arflg = (iqa lt 0.98*iqa_ms_large)*(ra lt 0.95*rsun) or (ra lt 0.95*rsun)*(abs(zqa0) gt bmax_qs)

	; QS/Non-sunspot is everywhere else on disk:
	qsflg = (arflg eq 0)*(ra lt ondisk_frac*rsun)
	qs_px = where(qsflg)
	ar_px = where(arflg)
	n_angles = n_elements(angles)

	zqa_cals = fltarr(nx,ny,n_angles)

	; Compute inclunation angles for each pixel:
	img_angles = (abs(asin(ra/rsun))) < max(angles)
	img_angles_mc = img_angles

	if(n_elements(mc_max_sunspot_angle) eq 1) then begin
		img_angles_mc[ar_px] = img_angles_mc[ar_px] < mc_max_sunspot_angle
	endif
	if(keyword_set(do_mc)) then begin
		zqa_cal = do_mc_cal(zqa, qsflg, arflg, img_angles, cal_measimages, cal_gtimages, cal_qsflags, cal_arflags, angles, qs_obsfits, qs_modfits, ar_obsfits, ar_modfits, img_angles_mc)
	endif else begin
		; Apply all curves at every angle:
		for i=0,n_angles-1 do begin
			zqa_cal = fltarr(nx,ny)	
			zqa_cal[qs_px] = interpol(qs_modfits[i,*],qs_obsfits[i,*],zqa[qs_px])		
			zqa_cal[ar_px] = interpol(ar_modfits[i,*],ar_obsfits[i,*],zqa[ar_px])
			zqa_cals[*,*,i] = zqa_cal
		endfor
		
		; Interpolate the all curves/every angle calibrated images down to the specific angle
		; for each pixel. Thus only the two nearest angles in the cal data are used for any
		; given pixel. A bit inefficient, but makes the code less cumbersome:
		zqa_cal = fltarr(nx,ny)
		testimg = fltarr(nx,ny)
		for i=0,nx-1 do begin
			for j=0,ny-1 do begin
				zqa_cal[i,j] = interpol(zqa_cals[i,j,*],angles,img_angles[i,j])
				testimg[i,j] = interpol(indgen(n_angles),angles,img_angles[i,j])
			endfor
		endfor
	endelse
	
	; Restore offdisk pixels to their original values:
	zqa_cal[where(offdisk)] = zqa[where(offdisk)]
		
	; Write calibrated images to output directory, using original headers:
	writefits, out_dir+file_basename(zqafile), zqa_cal, zqahdr, /compress
	
	; Diagnostic plots and printouts:
	print,'Calibrated net flux = ',total(zqa_cal[where(ondisk)]), ' Uncalibrated net flux = ',total(zqa0[where(ondisk)])
	print,'Calibrated AR net flux = ',total(zqa_cal[ar_px]), ' Uncalibrated AR net flux = ',total(zqa0[ar_px])
	print,'Calibrated QS net flux = ',total(zqa_cal[qs_px]), ' Uncalibrated QS net flux = ',total(zqa0[qs_px])

	if(n_elements(plot_dir) eq 1) then begin
		if(n_elements(prange) eq 0) then begin
			prange = max([mean(abs(zqa0[where(ondisk)]))+4.0*stdev(abs(zqa0[where(ondisk)])), $
						mean(abs(zqa_cal[where(ondisk)]))+4.0*stdev(abs(zqa_cal[where(ondisk)]))])
		endif
		set_plot,'z'
		device,set_resolution=round([3.9*nx,3.9*ny])
		!p.multi=[0,2,2]
		plot_image,zqa0*ondisk,min=-prange,max=prange,title='Original GONG image',charsize=2.5
		plot_image,zqa_cal*ondisk,min=-prange,max=prange,title='Curve-corrected GONG image',charsize=2.5
		cal_ratio = abs(zqa_cal)/abs(zqa)
		plot_image,cal_ratio,title='Ratio of corrected to uncorrected',min=0,max=2,charsize=2.5
		plot_image,qsflg,title='QS flags',charsize=2.5
		write_png,plot_dir+file_basename(zqafile,'.fits.gz')+'.png',tvrd()
	endif
	
	zqa_cal[zerovals] = 0.0

	return,zqa_cal
	
end
