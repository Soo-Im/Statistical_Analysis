;------------------------------------------
; PMIP3
; Regression map (monthly) to Nino3:JJA
; CCSM4
;------------------------------------------

begin

name		=	(/"CCSM4"/)
subject		=	(/"midHolocene","piControl"/)

; season select
Ninosson	=	"Sep"
Ninomon		=	8
Ninoreg		=	"Nino3"
sson		=	(/"Sep","Oct","Nov","Dec"/)
cmon		=	(/8,9,10,11/)

rad		=	4.*atan(1.)/180.		; radian = pi/180

; Region
sstS	=	-40.
sstN	=	 40.
sstL	=	-90.
sstR	=	 20.	; needs lonPivot(240)--> 120W~240E

do subj	= 0,1						; loop for subject
do expm = 0,dimsizes(name)-1		; loop for experiment

;---------------------------------------------------
;	descript	=	"SLP and wind-speed"
	descript	=	"SLP and 700-omega"

	; read tos data
	variable	=	"psl"
	file_name3	=	systemfunc("ls /storage/CMIP5_raw/"+subject(subj)+"/"+name(expm)+"/"+variable+"/*regrid*.nc")	
	file_in3	=	addfile(file_name3,"r")
	print("File: "+file_name3)
	psl_pivot	=	lonFlip(file_in3->psl(:,{sstS:sstN},:))
	psl_raw		=	psl_pivot(:,:,{sstL:sstR})
	delete([/file_name3,file_in3,psl_pivot/])

	variable	=	"tauu"
	file_name3	=	systemfunc("ls /storage/CMIP5_raw/"+subject(subj)+"/"+name(expm)+"/"+variable+"/*regrid*.nc")	
	file_in3	=	addfile(file_name3,"r")
	print("File: "+file_name3)
	var_pivot	=	lonFlip(file_in3->tauu(:,{sstS:sstN},:))
	wsu_raw		=	var_pivot(:,:,{sstL:sstR})	
	delete([/file_name3,file_in3,var_pivot/])
	
	variable	=	"tauv"
	file_name3	=	systemfunc("ls /storage/CMIP5_raw/"+subject(subj)+"/"+name(expm)+"/"+variable+"/*regrid*.nc")	
	file_in3	=	addfile(file_name3,"r")
	print("File: "+file_name3)
	var_pivot	=	lonFlip(file_in3->tauv(:,{sstS:sstN},:))
	wsv_raw		=	var_pivot(:,:,{sstL:sstR})	
	delete([/file_name3,file_in3,var_pivot/])

	; wind speed
	; "Anomaly of wind speed", not "speed of anomaly"
;	wap_raw		=	wind_speed(wsu_raw,wsv_raw)
	delete([/wsu_raw,wsv_raw/])

	; 500-omega
	variable	=	"wap"
	file_name3	=	systemfunc("ls /storage/CMIP5_raw/"+subject(subj)+"/"+name(expm)+"/"+variable+"/*regrid*.nc")	
	file_in3	=	addfile(file_name3,"r")
	print("File: "+file_name3)
	wap_pivot		=	lonFlip(file_in3->wap(:,{70000},{sstS:sstN},:))
	wap_raw			=	wap_pivot(:,:,{sstL:sstR})
	delete([/file_name3,file_in3,wap_pivot/])

	; land mask
	in		=   addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
	lsdata  =   in->LSMASK
	lsm 	=   landsea_mask(lsdata,psl_raw&lat,psl_raw&lon)
;	psl_raw	=   mask(psl_raw,conform(psl_raw,lsm,(/1,2/)).eq.1,0)
;	wap_raw	=   mask(wap_raw,conform(wap_raw,lsm,(/1,2/)).eq.1,0)	
	delete([/in,lsdata,lsm/])

	; detrend
	psl_raw		=	dtrend_msg_n(psl_raw&time,psl_raw,False,False,0)
	wap_raw		=	dtrend_msg_n(wap_raw&time,wap_raw,False,False,0)

	; remove seasonality
	anomalypsl	=	rmMonAnnCycTLL(psl_raw)	
	anomalywap	=	rmMonAnnCycTLL(wap_raw)	
	delete([/psl_raw,wap_raw/])


;---------------------------------------------------
	; assign time(month) to index and variable field
	dsize	=	dimsizes(anomalywap)
	ntim	=	dsize(0)
	nlat	=	dsize(1)
	nlon	=	dsize(2)
	
	anomalywap!0	=	"month"
	anomalypsl!0	=	"month"
	anomalywap&month	=	ispan(0,ntim-1,1)
	anomalypsl&month	=	ispan(0,ntim-1,1)

;---------------------------------------------------
	; Get Nino index
	; Nino 3 ( area average over 150W-90W, 5S-5N )
	ninoS	=	 -5.
	ninoN   =     5.
	ninoL	=    -150.			 	
	ninoR	=	 -90.
	Ninoreg	=	"Nino3"
	
	; read SST - due to different region with readed data
	file_name3	=	systemfunc("ls /storage/PMIP3/"+subject(subj)+"/"+name(expm)+"/tos/*regrid*.nc")	
	file_in3	=	addfile(file_name3,"r")
	print("File: "+file_name3)
	nino_pivot	=	lonFlip(file_in3->tos(:,{ninoS:ninoN},:))
	nino_raw	=	nino_pivot(:,:,{ninoL:ninoR})	

	; detrend
	nino_raw	=	dtrend_msg_n(nino_raw&time,nino_raw,False,False,0)

	; remove seasonality
	anomalynino	=	rmMonAnnCycTLL(nino_raw)	

	; weighted area average
	ninolat		=	file_in3->lat({ninoS:ninoN})
	ninocost	= 	cos(ninolat*rad) 
	nino		=	wgt_areaave_Wrap(anomalynino(:,{ninoS:ninoN},{ninoL:ninoR}),ninocost,1.0,1)	
	; variable, weight y, weight x, Opt(0=missingvalue allowed, 1=not allowed)

	nino!0		=	"month"
	nino&month	=	ispan(0,ntim-1,1)

	; extract a season from Nino index
;	sson_nino	=	month_to_season(nino,Ninosson)

	; extract a month from Nino index
	ninoind		=	ind(mod(nino&month,12).eq.Ninomon)
	sson_nino	=	nino(ninoind)
	delete([/ninoind/])

	delete([/file_name3,file_in3,nino_pivot/])
	delete([/nino_raw,anomalynino,nino/])

;---------------------------------------------------
	; extract seasonal variable fields


	wks1_name	=	"~/TraCE/figure/PMIP/reg_Nino/new/"+name(expm)+"_"+subject(subj)+"_regmap_"+Ninoreg+"_"+Ninosson+"_clt_and_psl"
	wks1		=	gsn_open_wks("x11",wks1_name)
	print("plot= "+wks1_name)
	plot1		=	new(dimsizes(sson),graphic)

	res     				=   True
	res@gsnRightString		=	"Reg["+Ninoreg+": "+Ninosson+"(0)]"
	

	do season	=	0,dimsizes(sson)-1
			
		print("month: "+sson(season)+" index: "+cmon(season))
		
		mthind	=	ind(mod(anomalypsl&month,12).eq.cmon(season))		; extract a month
		psl		=	anomalypsl(mthind,:,:)
		wap		=	anomalywap(mthind,:,:)
		delete([/mthind/])
;		printVarSummary(psl)
			
;		if (season.le.1) then
			print(sson(season)+" :same year")
			res@gsnLeftString	=	sson(season)+"(0)"
			rc1	=	regCoef_n(sson_nino,psl,0,0)
			rc2	=	regCoef_n(sson_nino,wap,0,0)
;		else
;			print(sson(season)+" :post year")
;			res@gsnLeftString	=	sson(season)+"(1)"
;			rc1	=	regCoef_n(sson_nino(:dimsizes(sson_nino)-2),psl(1:,:,:),0,0)
;			rc2	=	regCoef_n(sson_nino(:dimsizes(sson_nino)-2),wap(1:,:,:),0,0)
;		end if

		copy_VarCoords(psl(0,:,:),rc1)
		copy_VarCoords(psl(0,:,:),rc2)

		; significance , 0.95, two-tailed	
		df1		=	onedtond(rc1@nptxy,dimsizes(rc1))-2		
		tval1	=	onedtond(rc1@tval,dimsizes(rc1))		
		b1		=	tval1
		b1		=	0.5
		prob1	=	betainc(df1/(df1+tval1^2),df1/2.0,b1)
		
		df2		=	onedtond(rc2@nptxy,dimsizes(rc2))-2		
		tval2	=	onedtond(rc2@tval,dimsizes(rc2))		
		b2		=	tval2
		b2		=	0.5
		prob2	=	betainc(df2/(df2+tval2^2),df2/2.0,b2)

		rc951	=	where(prob1.le.0.05,rc1,rc1@_FillValue)
		rc952	=	where(prob2.le.0.05,rc2,rc2@_FillValue)

		copy_VarCoords(rc1,rc951)
		copy_VarCoords(rc2,rc952)
	
		; Panel frame
		res@gsnDraw			=	False
		res@gsnFrame		=	False

		; map
		res@mpMinLatF		=	sstS
		res@mpMaxLatF		=	sstN
		res@gsnAddCyclic	= 	False
		res@mpMinLonF		=	sstL
		res@mpMaxLonF		=	sstR	
		res@mpCenterLonF	=	(sstL+sstR)/2
		res@mpFillOn		=	False

		;------	
		; resource for rc1
		; shade
		res@cnFillOn	 	=	True
		res@cnLineLabelsOn	=	False
		res@cnLinesOn		=	False
		res@cnFillPalette   =   "BlueRed"	;MPL_Reds or BlueRed

		res@cnLevelSelectionMode	=	"ManualLevels"
		res@lbLabelBarOn	=	False
		res@cnMaxLevelValF	=	 80.
		res@cnMinLevelValF	=	-80.
		res@cnLevelSpacingF	=	 10.


		;------	
		; resource for rc2
		sres	=	True
		sres@gsnDraw		=	False
		sres@gsnFrame		=	False

		; map
		sres@gsnAddCyclic	= 	False
		sres@mpFillOn		=	False	

		sres@cnLevelSelectionMode	=	"ManualLevels"
		; wind-speed
;		sres@cnMaxLevelValF		=	 0.02 
;		sres@cnMinLevelValF		=	-0.02
;		sres@cnLevelSpacingF	=	 0.002
		; 500-omega
		sres@cnMaxLevelValF		=	 0.04 
		sres@cnMinLevelValF		=	-0.04
		sres@cnLevelSpacingF	=	 0.004

		; contour line
		sres@cnLineLabelsOn		=	True
;		sres@cnInfoLabelOn		=	False
		sres@gsnContourZeroLineThicknessF	=	2.
		sres@gsnContourNegLineDashPattern	=	1

		sres@cnLineLabelBackgroundColor		=	"white"
		sres@cnConstFLabelFontHeightF		=	2.
		sres@cnConstFLabelFontThicknessF	=	2.


		plot1(season)	=	gsn_csm_contour_map_overlay(wks1,rc951,rc952,res,sres)
		
		delete([/psl,wap,rc1,rc2,prob1,prob2,rc951,rc952/])	

	
	end do 	; loop for season

	resP					=	True
	resP@gsnMaximize		=	True
	resP@gsnPanelLabelBar	=	True
	resP@gsnPanelMainString	=	name(expm)+", "+subject(subj)+" "+descript
	resP@gsnPanelYWhiteSpacePercent	=	5.
	
	gsn_panel(wks1,plot1,(/2,2/),resP)

	delete([/anomalywap,anomalypsl/])
	delete([/sson_nino/])
	delete([/wks1,plot1,res/])
	
end do	; loop for experiment
end do		; loop for subject

end


;		; Second, remove r2y*Nino from SST
;		dimsize	=	dimsizes(sst)
;		ntim	=	dimsize(0)
;		nlat	=	dimsize(1)
;		nlon	=	dimsize(2)
;		
;		new_sst	=	new((/ntim,nlat,nlon/),typeof(sst))
;		
;		do t = 0,ntim-1
;			new_sst(t,:,:)	=	sst(t,:,:)-r2y(:,:)*sson_nino(t)
;		end do
;
;		; Third, regression of Atl to NEW SST
;		if (sson(season).eq."DJF") then		; Atl(JJA,0) - SST(DJF,+1)
;			rc		=	regCoef_n(sson_nta(0:ntim-2),new_sst(1:ntim-1,:,:),0,0)
;;			rcu		=	regCoef_n(sson_nta(0:ntim-2),new_u(1:ntim-1,:,:),0,0)
;;			rc2		=	regCoef_n(sson_nta(0:ntim-2),new_v(1:ntim-1,:,:),0,0)
;		else								; Alt(JJA,0) - SST(JJA,SON,0)
;			rc		=	regCoef_n(sson_nta,new_sst,0,0)
;;			rcu		=	regCoef_n(sson_nta,new_u,0,0)
;;			rc2		=	regCoef_n(sson_nta,new_v,0,0)
;		end if
	;		rc	=	regCoef_n(sson_nta,new_sst,0,0)


	
;
;/;	
;	; NTA ( 90W-20E, 0-15N )
;	ntaS	=	 0.
;	ntaN	=	 15.
;	ntaL	=	-90.
;	ntaR	=	 20.
;	Atlreg	=	"NTA"
;;/
;
;	; Atl3 ( 20W-0E, 3S-3N )
;	ntaS	=	-3
;	ntaN	=	 3
;	ntaL	=	-20
;	ntaR	=	 0	
;	Atlreg	=	"Atl3"
;
;	; weigthed area average	
;	atllat	=	file_in1->lat({ntaS:ntaN})
;	atlcost	= 	cos(atllat*rad) 
;		 
;	nta		=	wgt_areaave_Wrap(anomalywap(:,{ntaS:ntaN},{ntaL:ntaR}),atlcost,1.0,1)	
;	; variable, weight y, weight x, Opt(0=missingvalue allowed, 1=not allowed)

;	file_name2	=	"~/index/"+name(expm)+"_"+subject(subj)+"_"+Ninoreg+".nc"	
;	file_name3	=	"~/index/"+name(expm)+"_"+subject(subj)+"_"+Atlreg+".nc"	
;	file_in2	=	addfile(file_name2,"c")
;	file_in3	=	addfile(file_name3,"c")
;	file_in2->var	=	nino
;	file_in3->var	=	nta
;/
