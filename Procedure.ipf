#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.


function msd_track(w, id) //fps
// Calculation of mean square displacement
	wave w    // trajectory
	variable id    // trajectory id

	variable num=numpnts(w)
	variable num_msd = round(num/2)
	variable fps=29.1

	string msdrname="msd"+num2str(id)	//start to calculate msd
	make/o/n=(num_msd) msdr
	make/o/n=(num) sqr_dis_r
	variable ii, lag
	for (ii=0; ii<num_msd; ii+=1)
		for (lag=0;lag+ii<num-1;lag+=1)
			sqr_dis_r[lag]=(w[lag+ii]-w[lag])^2
		endfor
		msdr[ii]=mean(sqr_dis_r)
	endfor


	setscale/p x, 0, 1/fps, "s", msdr
	setscale/p y, 0, 1, "nm\S2", msdr
	duplicate/o msdr $msdrname
	killwaves sqr_dis_r, msdr

End


function calc_chi(number)
// Calculation of chi
// The name of particle trajectories begin with "trj"

	variable number
	
	make/o/n=(number) v = 0, chi_star = 0, v_fit = 0
	variable num_chi = 11
	make/o/n=(num_chi) temp_w_chi
	variable mean_chi=0
	variable fps=29.1
	variable ic=0, ii, num
	
	string ListOfWaves=WaveList("trj*", ";", "")
	string wname
	for(ic=0; ic<number; ic+=1)

		duplicate/o $StringFromList(ic, ListOfWaves, ";") w
		num = numpnts(w)
 
		for(ii=0; ii<num_chi; ii+=1)
			make/o/n=(num-ii) w_chi = w[p+ii]-w[p]
			temp_w_chi[ii]=2*mean(w_chi)/variance(w_chi)
		endfor
		
		setscale/p x, 0, 1/fps, "", temp_w_chi
		wname = "chi_"+ num2str(ic)
		duplicate/o temp_w_chi $wname
		
		CurveFit/Q/M=2/W=0 line, w/D
		wave W_coef
		v_fit[ic] = W_coef[1]
		v[ic]=(w[num-1]-w[0])/num*fps
		make/o/free/n=5 tmp_c = temp_w_chi[p+8]
		chi_star[ic] = mean(tmp_c)
		
	endfor	

end


Function CalcChi_from_Graph(graphNameStr)
// Calculate chi* from multiple chi plots.

    String graphNameStr         // "" for top graph
   
    String list = TraceNameList(graphNameStr,";",1) // Ignore traces that belong to contour plots
    Variable numTraces = ItemsInList(list)
    make/o/n=(numTraces) chi_star_endpoint = 0, chi_star_endpoint_vel = 0
    make /o/T/n=(numTraces) wname
    Variable i
    variable fps=29.1, num
    string vname
    for(i=0; i<numTraces; i+=1)
        String traceNameStr = StringFromList(i, list)
        Wave w = TraceNameToWaveRef(graphNameStr, traceNameStr )
        make/o/free/n=3 w_t=w[p+8]
        chi_star_endpoint[i] = mean(w_t)
        vname = "trj"+traceNameStr[4,inf]
        duplicate/o $vname trj_tmp
        num = numpnts(trj_tmp)
        chi_star_endpoint_vel[i]=(trj_tmp[num-1]-trj_tmp[0])/(num-0)*fps
        wname[i] = traceNameStr[4,inf]
    endfor
    
    mk_avg_v(chi_star_endpoint, chi_star_endpoint_vel)
End

function mk_avg_v(w_c, w_v)
// Averaging the chi* and corresponding velocity
	wave w_c  // wave of chi*
	wave w_v  // wave of velocity
	
	variable num=numpnts(w_c)
	make/o/n=5 chi_vel_mean=0, nnn=0, chi1_mean=0, chi_var=0, chi_vel_var=0
	
	variable ii, idx
	for (ii=0; ii<num; ii+=1)
		idx = round(10*w_c[ii])	
		if (idx == 0) 
			idx = 1
		endif	
		chi_vel_mean[idx] = (chi_vel_mean[idx]*nnn[idx]+w_v[ii])/(nnn[idx]+1)
		chi_vel_var[idx] = (chi_vel_var[idx]*nnn[idx]+w_v[ii]*w_v[ii])/(nnn[idx]+1)
		chi1_mean[idx] = (chi1_mean[idx]*nnn[idx]+w_c[ii])/(nnn[idx]+1)
		chi_var[idx] = (chi_var[idx]*nnn[idx]+w_c[ii]*w_c[ii])/(nnn[idx]+1)
		nnn[idx] += 1
		
		print w_c[ii], idx
	endfor
	for (idx=0; idx<5; idx+=1)
		if (chi_vel_mean[idx] == 0)
			chi_vel_mean[idx] = NaN
			chi_vel_var[idx] = NaN
		endif
	
	endfor
	duplicate/o chi_vel_var chi_vel_std
	chi_vel_std = sqrt(chi_vel_var-chi_vel_mean^2)
	duplicate/o chi_var chi_std
	chi_std = sqrt(chi_var-chi1_mean^2)
end

