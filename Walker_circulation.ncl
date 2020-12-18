;------------------------------------------
; Seasonal mean Walker circulation (latitudal mean of w, u)
; Use CCSM4 model results
;------------------------------------------

begin

name		  =	(/"CCSM4"/)
subject		=	(/"midHolocene","piControl"/)

; season select
sson		=	(/"DJF","JJA"/)

rad		=	4.*atan(1.)/180.		; radian = pi/180

; Tropical Region
sstS	=	-5.
sstN	=	 5.
sstL	=	-180.
sstR	=	 180.	

do subj	= 0,1						; loop for subject
  ;------------Data Handling------------------
  ;---Read Omega, Zonal wind
	variable	=	"ua"
	file_name3	=	systemfunc("ls /storage/CMIP5_raw/"+subject(subj)+"/"+name+"/"+variable+"/*regrid*.nc")	
	file_in3	=	addfile(file_name3,"r")
	u_raw		=	lonFlip(file_in3->ua(:,{:30000},{sstS:sstN},:))
	
	variable	=	"wap"
	file_name3	=	systemfunc("ls /storage/CMIP5_raw/"+subject(subj)+"/"+name+"/"+variable+"/*regrid*.nc")	
	file_in3	=	addfile(file_name3,"r")
	v_raw		=	lonFlip(file_in3->wap(:,{:30000},{sstS:sstN},:))

	;---average over latitude ( time x lev x lon )
	varlat		=	  u_raw&lat
	varcost		= 	cos(varlat*rad) 
	opt		=	1					; average of non-missing value
	u_lon		=	dim_avg_wgt_n_Wrap(u_raw,varcost,opt,2)	
	v_lon		=	dim_avg_wgt_n_Wrap(v_raw,varcost,opt,2)	
	;printVarSummary(u_lon)
	delete([/varlat,varcost/])
	delete([/u_raw,v_raw/])

	;---detrend over time
	u_lon		=	dtrend_msg_n(u_lon&time,u_lon,False,False,0)
	v_lon		=	dtrend_msg_n(v_lon&time,v_lon,False,False,0)

	;---scale omega by u (omega << u)
	uave		=	abs(avg(u_lon(0,:,:)))
	vave		=	abs(avg(v_lon(0,:,:)))
	scale		=	uave/vave
	v_lon		=	v_lon*scale
;	printVarSummary(v_lon)

	;---assign time to index and variable field
	dsize	=	dimsizes(u_lon)
	ntim	=	dsize(0)
	nlev	=	dsize(1)
	nlon	=	dsize(2)

	u_lon!0		=	"month"
	v_lon!0		=	"month"
	u_lon&month	=	ispan(0,ntim-1,1)
	v_lon&month	=	ispan(0,ntim-1,1)

;----------------Plot-----------------------
	; plot directory
	wks1_name	=	"~/figure/"+name+"_"+subject(subj)+"_Walker"
	wks1		=	gsn_open_wks("x11",wks1_name)
	print("plot= "+wks1_name)
	plot1		=	new(dimsizes(sson)+1,graphic)	; one for annual mean

	res     	=   True

 	 ; Panel frame
	res@gsnDraw	=	False
	res@gsnFrame	=	False

	; shade
	res@cnFillOn		=	True
	res@cnFillPalette  	=   "BlueWhiteOrangeRed"	;MPL_Reds or BlueRed
	res@cnLineColor		=	"Gray"

  	; shade intervals
	res@cnLevelSelectionMode	=	"ManualLevels"
	res@cnMaxLevelValF	=	 3.6
	res@cnMinLevelValF	=	-3.6
	res@cnLevelSpacingF	=	 0.4
	res@lbLabelBarOn	=	False


	; loop for seasonal mean plot
	  do season	=	0, dimsizes(sson)-1
		
		  print("season: "+sson(season))
		  u_mean	=	dim_avg_n_Wrap(month_to_season(u_lon,sson(season)),0)
		  v_mean	=	dim_avg_n_Wrap(month_to_season(v_lon,sson(season)),0)

	  	; title
	  	res@gsnLeftString	=	sson(season)
	  	res@gsnRightString	=	"e-11"

    		  ; plot a season
  		plot1(season+1)		=	gsn_csm_pres_hgt_vector(wks1,v_ann,u_ann,v_ann,res)

	  	delete([/u_mean,v_mean/])	
	
	  end do 	; loop for season


  	; Frame and draw Panels
	resP			=	True
	resP@gsnMaximize	=	True
	resP@gsnPanelLabelBar	=	True
	resP@gsnPanelMainString	=	subject(subj)+" Walker cell"
	resP@gsnPanelYWhiteSpacePercent	=	5.
	
	gsn_panel(wks1,plot1,(/2,2/),resP)

	delete([/u_lon,v_lon/])
	delete([/wks1,plot1,res/])
	
end do		; loop for subject

end
