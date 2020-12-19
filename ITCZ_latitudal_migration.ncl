;------------------------------------------
; Seasonal ITCZ migration: Maximum Precipitation region
; latitude x calendar month
; CCSM4
;------------------------------------------

begin

name			=	"CCSM4"
subject		=	(/"piControl","midHolocene"/)			; experiments

rad				=	4.*atan(1.)/180.									; radian = pi/180

; Region
sstS	=	-15.
sstN	=	 15.
sstL	=	-60.
sstR	=	 20.	

; ITCZ location storage to plot both experiments at once
itcz	=	new((/dimsizes(subject),12/),float)

do subj	= 0,1						; loop for subject

; Don't remove the climatology!
; We are finding seasonal mean.

;-------------------------------------------------
	;--------data handling
	; Precipitation
	variable	=	"pr"	
	file_name3	=	systemfunc("ls /storage/CMIP5_raw/"+subject(subj)+"/"+name+"/"+variable+"/*regrid*.nc")	
	file_in3	=	addfile(file_name3,"r")
	print("File: "+file_name3)
	var_pivot	=	lonFlip(file_in3->pr(:,{sstS:sstN},:))	
	prc_raw		=	var_pivot(:,:,{sstL:sstR})	
	delete([/var_pivot/])

	; dtrend
	prc_raw		=	dtrend_msg_n(prc_raw&time,prc_raw,False,False,0)

	; land mask
	in		=   addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
	lsdata  =   in->LSMASK
	lsm 	=   landsea_mask(lsdata,prc_raw&lat,prc_raw&lon)
	prc_raw	=   mask(prc_raw,conform(prc_raw,lsm,(/1,2/)).eq.1,0)
	delete([/in,lsdata,lsm/])	
	
	;--------get climatology (monthly mean)
	; Climatology of calendar month
	prc_clm		=	clmMonTLL(prc_raw)
	
	; zonal mean of climatology
	prc_zon		=	dim_avg_n_Wrap(prc_clm,2)	
 
	; assign time coordinate
	prc_zon!0	=	"month"
	prc_zon&month	=	ispan(1,12,1)

	delete([/prc_raw,prc_clm/])

;-------------------------------------------------
	;--------plot climatological precipitation
	wks1_name	=	"~/"+name+"_"+subject(subj)+"_zonal_mean_precipitation"
	wks1			=	gsn_open_wks("x11",wks1_name)
	print("plot= "+wks1_name)
	
	res     	=   True
	
	; shade
	res@cnFillOn	 			=	True
	res@cnFillOn	 			=	True
	res@cnLineLabelsOn	=	False
	res@cnLinesOn				=	False
	res@cnFillPalette   = "MPL_Oranges"	
;	res@cnLineColor			=	"Gray"

	; title
	res@tiMainString	=	"Zonal mean Precipitation"	
	res@gsnLeftString	=	name+" "+subject(subj)
	res@gsnRightString	=	""	

	; axis
	res@tmXBMode			=	"Explicit"
	res@tmXBValues		=	prc_zon&month
	res@tmXBLabels		=	prc_zon&month
	res@tiXAxisString	=	"Calendar month"
	res@tiYAxisString	=	""

	plot1	=	gsn_csm_contour(wks1,prc_zon(lat|:,month|:),res)

;-------------------------------------------------
	; find maximum index along the latitude for each calendar month
	; dim_maxind function can take only dim=0, so we need to exchange lat and month dimension.

	;--------find maximum precipitation location
	prc_zon_inv		=	prc_zon(lat|:,month|:)										; exchange lat and month dimension
	itcz(subj,:)	=	prc_zon&lat( dim_maxind(prc_zon_inv,0) )
	
	print(itcz(subj,:))

	delete([/prc_zon/])
	delete([/prc_zon_inv/])

end do

delete([/wks1_name,wks1,res,plot1/])
	
	
;-------------------------------------------------
;--------plot ITCZs in two experiments at once
itcz!1			=	"month"
itcz&month	=	ispan(1,12,1)

wks1_name	=	"~/ITCZ_migration"
wks1			=	gsn_open_wks("png",wks1_name)
print("plot= "+wks1_name)
	
; plot resource
res			=	True

; Axis

res@tiYAxisString	=	""
res@trYMinF			=	-7.
res@trYMaxF			=	 13.

res@tiXAxisString	=	"Calnedar Month"
res@tmXBMode		=	"Explicit"
res@tmXBValues		=	itcz&month
res@tmXBLabels		=	(/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
res@trXMinF			=	1
res@trXMaxF			=	12
res@tmXBLabelAngleF		=	-45		; tilt x label

; title
res@tiMainString	=	"Seasonality of ITCZ"
;res@gsnRightString	=	"Latitude of Maximum zonal mean convection"
res@gsnRightString	=	"Latitude of Minimum meridional wind"

; Line
res@xyMarkLineMode	=	"MarkLines"
res@xyDashPatterns	=	(/1,0/)
res@xyMarkers		=	(/4,4/)
res@xyLineColors	=	(/"tomato","dodgerblue"/)	
res@xyLineThicknesses	=	(/3.,3./)

; Grid
res@tmXMajorGrid	=	True
res@tmYMajorGrid	=	True
res@tmXMajorGridThicknessF		=	1.0
res@tmYMajorGridThicknessF		=	1.0
res@tmXMajorGridLineDashPattern	=	2
res@tmYMajorGridLineDashPattern	=	2

; Legend
res@pmLegendDisplayMode    = "Always"            ; turn on legend
res@pmLegendSide           = "Top"               ; Change location of 
res@pmLegendParallelPosF   =  0.2                  ; move units right
res@pmLegendOrthogonalPosF = -0.3                ; move units down
res@pmLegendWidthF         = 0.15                ; Change width and
res@pmLegendHeightF        = 0.10                ; height of legend.
res@lgPerimOn              = False               ; turn off box around
res@lgLabelFontHeightF     = .01                 ; label font height
res@xyExplicitLegendLabels	=	(/subject(0),subject(1)/)

plot1	=	gsn_csm_xy(wks1,itcz&month,itcz,res)


end
