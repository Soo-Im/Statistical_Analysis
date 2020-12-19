;------------------------------------------
; Regression map with significance test (t-test)
; Regression of 'monthly anomalous variables' to 'monthly Nino3 index'
; Ex) I want to find out the change of SST anomalies which varies with Nino3. 
; CCSM4
;------------------------------------------

begin

name		=	(/"CCSM4"/)
subject		=	(/"midHolocene","piControl"/)

; Select Season
Ninosson	=	"Sep"				; Choose September Nino3
Ninomon		=	8				; Index starts from 0 thus September = 8th
Ninoreg		=	"Nino3"
sson		=	(/"Sep","Oct","Nov","Dec"/)	; Choose Sep~Dec Variables
cmon		=	(/8,9,10,11/)

rad		=	4.*atan(1.)/180.		; radian = pi/180

; Region
sstS	=	-40.
sstN	=	 40.
sstL	=	-90.
sstR	=	 20.	

do subj	= 0,1		; loop for two experiments

	descript	=	"SLP and 700-omega"

;---------------------------------------------------
	;------data handling
	; Sea Level Pressure
	variable	=	"psl"
	file_name3	=	systemfunc("ls /storage/"+subject(subj)+"/"+name+"/"+variable+"/*regrid*.nc")	
	file_in3	=	addfile(file_name3,"r")
	print("File: "+file_name3)
	psl_pivot	=	lonFlip(file_in3->psl(:,{sstS:sstN},:))
	psl_raw		=	psl_pivot(:,:,{sstL:sstR})
	delete([/file_name3,file_in3,psl_pivot/])

	; Omega (Vertical Temp) at 700mb
	variable	=	"wap"
	file_name3	=	systemfunc("ls /storage/"+subject(subj)+"/"+name+"/"+variable+"/*regrid*.nc")	
	file_in3	=	addfile(file_name3,"r")
	print("File: "+file_name3)
	wap_pivot	=	lonFlip(file_in3->wap(:,{70000},{sstS:sstN},:))
	wap_raw		=	wap_pivot(:,:,{sstL:sstR})
	delete([/file_name3,file_in3,wap_pivot/])

	; land mask
	in		=   addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
	lsdata  =   in->LSMASK
	lsm 	=   landsea_mask(lsdata,psl_raw&lat,psl_raw&lon)
	psl_raw	=   mask(psl_raw,conform(psl_raw,lsm,(/1,2/)).eq.1,0)
	wap_raw	=   mask(wap_raw,conform(wap_raw,lsm,(/1,2/)).eq.1,0)	
	delete([/in,lsdata,lsm/])

	; detrend
	psl_raw		=	dtrend_msg_n(psl_raw&time,psl_raw,False,False,0)
	wap_raw		=	dtrend_msg_n(wap_raw&time,wap_raw,False,False,0)

	; remove seasonality (monthly anomaly)
	anomalypsl	=	rmMonAnnCycTLL(psl_raw)	
	anomalywap	=	rmMonAnnCycTLL(wap_raw)	
	delete([/psl_raw,wap_raw/])

	; assign time(month) to index and variable field
	dsize	=	dimsizes(anomalywap)
	ntim	=	dsize(0)
	nlat	=	dsize(1)
	nlon	=	dsize(2)
	
	anomalywap!0		=	"month"
	anomalypsl!0		=	"month"
	anomalywap&month	=	ispan(0,ntim-1,1)
	anomalypsl&month	=	ispan(0,ntim-1,1)

;---------------------------------------------------
	;------Nino index
	; Nino3: area averaged SST anomalies over 150W-90W, 5S-5N 
	ninoS	=	-5.
	ninoN   = 	 5.
	ninoL	=	-150.			 	
	ninoR	=	-90.
	Ninoreg	=	"Nino3"
	
	; read SST 
	file_name3	=	systemfunc("ls /storage/PMIP3/"+subject(subj)+"/"+name+"/tos/*regrid*.nc")	
	file_in3	=	addfile(file_name3,"r")
	nino_pivot	=	lonFlip(file_in3->tos(:,{ninoS:ninoN},:))
	nino_raw	=	nino_pivot(:,:,{ninoL:ninoR})	

	; detrend
	nino_raw	=	dtrend_msg_n(nino_raw&time,nino_raw,False,False,0)

	; remove seasonality
	anomalynino	=	rmMonAnnCycTLL(nino_raw)	

	; weighted area average
	ninolat		=	file_in3->lat({ninoS:ninoN})
	ninocost	= 	cos(ninolat*rad) 
	Opt		=	1			; wgt_areaave_Wrap Option (0=missingvalue allowed, 1=not allowed)
	nino		=	wgt_areaave_Wrap(anomalynino(:,{ninoS:ninoN},{ninoL:ninoR}),ninocost,1.0,Opt)	

	; extract a month from Nino index
	nino!0		=	"month"
	nino&month	=	ispan(0,ntim-1,1)
	ninoind		=	ind(mod(nino&month,12).eq.Ninomon)
	sson_nino	=	nino(ninoind)
	delete([/ninoind/])

	delete([/file_name3,file_in3,nino_pivot/])
	delete([/nino_raw,anomalynino,nino/])

;---------------------------------------------------
	;------Panel plot of Monthly Regression maps
	wks1_name	=	"~/"+name+"_"+subject(subj)+"_regmap_"+Ninoreg+"_"+Ninosson+"_"+descript
	wks1		=	gsn_open_wks("x11",wks1_name)
	plot1		=	new(dimsizes(sson),graphic)
	print("plot= "+wks1_name)

	; Title
	res     		=   True
	res@gsnRightString	=	"Reg["+Ninoreg+": "+Ninosson+"(0)]"
	
	
	do season	=	0,dimsizes(sson)-1				; loop for months
			
		;------Regression
		; extract variables of each month
		print("month: "+sson(season)+" index: "+cmon(season))
		
		mthind	=	ind(mod(anomalypsl&month,12).eq.cmon(season))
		psl	=	anomalypsl(mthind,:,:)
		wap	=	anomalywap(mthind,:,:)
		delete([/mthind/])

		; Regression variables to Nino3	
		print(sson(season)+" :same year")
		res@gsnLeftString	=	sson(season)+"(0)"
		rc1	=	regCoef_n(sson_nino,psl,0,0)
		rc2	=	regCoef_n(sson_nino,wap,0,0)

		; assign coordinates to regression coefficients
		copy_VarCoords(psl(0,:,:),rc1)
		copy_VarCoords(psl(0,:,:),rc2)

		; significance test (0.95, two tailed t-test)
		df1	=	onedtond(rc1@nptxy,dimsizes(rc1))-2		
		tval1	=	onedtond(rc1@tval,dimsizes(rc1))		
		b1	=	tval1
		b1	=	0.5
		prob1	=	betainc(df1/(df1+tval1^2),df1/2.0,b1)
		
		df2	=	onedtond(rc2@nptxy,dimsizes(rc2))-2		
		tval2	=	onedtond(rc2@tval,dimsizes(rc2))		
		b2	=	tval2
		b2	=	0.5
		prob2	=	betainc(df2/(df2+tval2^2),df2/2.0,b2)

		; Remove regression coefficients where p>0.05
		rc951	=	where(prob1.le.0.05,rc1,rc1@_FillValue)
		rc952	=	where(prob2.le.0.05,rc2,rc2@_FillValue)
		copy_VarCoords(rc1,rc951)
		copy_VarCoords(rc2,rc952)
	
		;------Plot
		; Panel frame
		res@gsnDraw		=	False
		res@gsnFrame		=	False

		; map
		res@mpMinLatF		=	sstS
		res@mpMaxLatF		=	sstN
		res@gsnAddCyclic	= 	False
		res@mpMinLonF		=	sstL
		res@mpMaxLonF		=	sstR	
		res@mpCenterLonF	=	(sstL+sstR)/2
		res@mpFillOn		=	False


		; copy: res for shading, sres for contouring (two resource values)
		res	=	sres

		; resource for shading (res)
		res@cnFillOn	 	=	True
		res@cnLineLabelsOn	=	False
		res@cnLinesOn		=	False
		res@cnFillPalette   	=   "BlueRed"	
		res@cnLevelSelectionMode=	"ManualLevels"
		res@lbLabelBarOn	=	False
		res@cnMaxLevelValF	=	 80.
		res@cnMinLevelValF	=	-80.
		res@cnLevelSpacingF	=	 10.


		; resource for contour (sres)
		sres@cnLevelSelectionMode	=	"ManualLevels"
		sres@cnMaxLevelValF		=	 0.04 
		sres@cnMinLevelValF		=	-0.04
		sres@cnLevelSpacingF		=	 0.004
		sres@cnLineLabelsOn		=	True
;		sres@cnInfoLabelOn		=	False
		sres@gsnContourZeroLineThicknessF	=	2.
		sres@gsnContourNegLineDashPattern	=	1
		sres@cnLineLabelBackgroundColor		=	"white"
		sres@cnConstFLabelFontHeightF		=	2.
		sres@cnConstFLabelFontThicknessF	=	2.


		plot1(season)	=	gsn_csm_contour_map_overlay(wks1,rc951,rc952,res,sres)
		
		delete([/psl,wap,rc1,rc2,prob1,prob2,rc951,rc952/])	

	
	end do 

	; resource for panel (resP)
	resP				=	True
	resP@gsnMaximize		=	True
	resP@gsnPanelLabelBar		=	True
	resP@gsnPanelMainString		=	name(expm)+", "+subject(subj)+" "+descript
	resP@gsnPanelYWhiteSpacePercent	=	5.
	
	; Make a panel by 2*2
	gsn_panel(wks1,plot1,(/2,2/),resP)

	delete([/anomalywap,anomalypsl/])
	delete([/sson_nino/])
	delete([/wks1,plot1,res/])
	
end do

end
