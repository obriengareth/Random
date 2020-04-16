#include <DELPHI_IOc.h>

void least_filterc(float **datin, float **datref, float **datmask, int nt, int nx, int lfilt, int nrot, float eps, float noise, float ampmax, float *filter, float **datout);

static char rcsid[] = "$Id: EPSI_fixed.c,v 1.1 2011/12/09 12:19:46 seistool Exp seistool $";

/*********************** self documentation **********************/
char *sdoc[] = {
" ",
" EPSI for fixed spread input data",
" ",
" EPSI.fixed < stdin > stdout  [optional parameters]",
" ",
" Required parameters:",
" ",
" Optional parameters:",
" ",
"   mode= .................... 0: basic EPSI ; 1 with near offset reconstruction",
"   file_in= ................. input file (assume pipe if not set)",
"   file_prim= ............... estimated primaries (assume pipe if not set)",
"   file_mult= ............... estimated multiples",
"   file_res= ................ residual data",
"   file_x0= ................. estimated primary impulse responses (spikes)",
"   file_wav= ................ estimated wavelet after each iteration ",
"   file_QC = ................ QC-file: QC results of each iteration are saved in this file (assume pipe if not set) ",
"   file_QC_dx0 = ............ save the update of x0 (before sparseness is applied) for the QC-shots to a separate file ",
"   QC = 0.................... = 0: no QC, = 1: QC estimated primaries ",
"                              = 2: QC primaries + residual, ",  
"                              = 3: QC primaries + residual + multiples ",
"                              if mode=1, reconstructed data is added as first entry of QC",
"   QC_shots = nshots/2....... give list of shots for QC (middle shot if not specified) ",
"   niter= ................... number of iterations",
"   v0=1500 .................. velocity in water",
"   key_depth = .............. specify name of header with water depth. This is used to determine",
"                              the lower limit of the time-window which is used to avoid that ",
"                              non-causal events end up in the primary estimation. See it as a mute ",
"                              along the water bottom reflection that will be applied before the picking ",
"   t0=0.25  ................. apex time of first reflection events in the data. ",
"                              If key_depth is not set, t0 is used to determine the lower limit of",
"                              the time-window.  ",
"   tmax_rel = tmax/t0 = 3.0 . tmax  = tmax_rel*t0 ( or tmax = (tmax_rel*key_depth*2)/v0): ",
"                              To avoid that refraction information is filtered out by the time ",
"                              window (muted), tmax limits the muting  ",
"   size_twi = 0.1s .......... size of the initial time-window: upper limit = lower limt + size_twi ",
"   dtw = 0.048s ............. number of seconds that the time-window increases in each iteration     ",
"   lfilt = 41 ............... length (in nr. of samples) of the estimated wavelet",
"   nrot = (lfilt-1)/2 ....... nr. of negative samples in the estimated wavelet",
"   n1=from input ............ number of samples",
"   nxmax=512 ................ maximum number of traces in input file",
"   ntmax=1024 ............... maximum number of samples/trace in input file",
"   ntapoffset = 0 ........... number of the highest offset traces in the data are tapered",
"   ntapt = 41 ............... number of time samples in the data are tapered",
"   ntapx = 0 ................ number of spatial samples in the data are tapered",
"   verbose=0 ................ silent option; >0 display info",
" ",
"   parameters in case near offset reconstruction is included (mode =1):    ",
" ",
"   file_p_recon= ............ data with reconstructed near offset gap",
"   min_offset= .............. minimum offset of the original data (used to create near-offset gap)",
"   tgap= .................... maximum time of the near-offset gap, for t > tgap interpolated data is used",
" ",
"   parameters in case the time variant wavelet extension is included:    ",
"   nwindows = 1 ............. number of time-windows used for the time-variant wavelet extension",
"                              if nwindows = 1: time variant wavelets extension is not included ",
"   size_window = 1s ......... size of time-windows in seconds",
"   Qfactor = 120 ...........  Qfactor ",
"   freq_base = 75 ..........  together with the Qfactor this determines the adaption of the wavelet for deeper parts of the data ",
" ",
" ",
" author  : Rolf Baardman : 16-05-2010 (r.h.baardman@tudelft.nl)",
" product : ME",
" ",
NULL};
/**************** end self doc ***********************************/

void main(int argc, char *argv[])
{
	float **tmp2_x0;			
	intn 	seqnr[MAX_KEYS];
	int	type, dom1, dom2;
	short   trid;
	int     optn, sign, ntmax, nxmax, verbose, first, nkeys, nreceivers, nsamples, mode;
	int     v_shot_start, v_shot_end, ivshot, ipshot, ipsample, QC, iQC, QC_shot;
	int     index_out, index_out_shot, index_freq_start, index_receiver_start;
	int	index_shot_start, index_shot_end, index_sample_start, ntapx, ntapt, ispike, nspikes;
	int     nfreq, error, n1, n2, nshots, maxshots, ret, n1_orig, size, i, j, ishot, ikey;
	int	ifreq, iter, niter, lfilt, nrot, it, ix, ip, itin, itout, shot_start, shot_end;
	int     ireceiver, isample, *mask, nmask, ireceiver_start, ireceiver_end, iwindow;
	int     *dx0_shot, *dx0_receiver, *dx0_sample, ipeak, Npeak, sample_max, nstwi, mstwi, dstwi, Dp, Ngap, Ngapl, Ngapr, var_t0, iwav;
	int     sampspk, shotspk, recspk, itmp, itmppc, istart, istart2, istartpc, npcpcH, nsamp_gap;
	int	*QC_shots, nQC_shots, use_pipe, samp_max, nwindows, nsamp_window, ntapoffset;
	float   *rdata, *datin, *tmpdata, *p, *pc, *Preal, *Pim, *Pcreal, *Pcim, *PPcH, *datout, *corr2D, *xtaper, *ttaper, *dx0, *dx0_amp, **x0, **B, **x0s, **dpminds, **s, *s_tmp, **s_rev, **mask_s, **mask_alfa, *tmp_recgather, **tmp_B, **tmp_x0s;
	float   d1, d2, f1, f2, scale, tstart, t0, dt0, tmax_rel, size_window, *t0mask, t1, t2, t3, offset, v0, tmp, max, alfa, beta, tmp1, tmp2, ampspk, size_twi, dtw, freq_base, Qfactor, *offsettaper, *wortelT;
	float *rdata_pin,*tmpdata_pin,*datout_pinc, **mask_beta, min_offset, tgap, *no_orig, *pin, *pinc, *tmp_dpin, *dpin, *dpinc, *tmpdata_wav, *rdata_wav, *datout_wav;
	complex *cdata, *cdata_pinc, *cdata_wav;
	segyhdr	*hdrs_in, *hdrs_tot;
	char    *file_in, *file_prim, *file_res, *file_p_recon, *file_mult, *file_x0, *file_QC, *file_QC_dx0, *file_wav, *key,  *keys[MAX_KEYS];
	
        String key_depth;	/* key which describes water depth	*/
        String ktype;	/* ... its type			*/
        int kindex;		/* ... its index in hdr.h	*/
        Value kvalkey;	/* ... its value		*/
        float32 depth;
	double tmp100, param;

/* Reading input parameters */

	tstart = cputime();
	initargs(argc, argv);
	requestdoc(0);

	if(!getparint("verbose", &verbose)) verbose=0;
	if(!getparint("mode", &mode)) mode=0;
	if(!getparstring("file_in", &file_in)) {
		if (verbose) sawarn("parameter file_in not given, assume pipe");
		file_in = NULL;
	}
	if(!getparstring("file_prim", &file_prim)){
		if (verbose) sawarn("parameter file_prim not given, assume pipe");
		file_prim = NULL;
		use_pipe=1;
	}
	else {
	 use_pipe = 0;
	 }	
	if(!getparstring("file_res", &file_res)){
		if (verbose) sawarn("parameter file_res not given, residual data not saved");
		file_res = NULL;
	}
	if(!getparstring("file_p_recon", &file_p_recon) ){
		if (verbose && mode==1) sawarn("parameter file_p_recon not given, reconstructed data not saved");
		file_p_recon = NULL;
	}
	if(!getparstring("file_mult", &file_mult)){
		if (verbose) sawarn("parameter file_mult not given, estimated multiples not saved");
		file_mult = NULL;
	}
	if(!getparstring("file_x0", &file_x0)){
		if (verbose) sawarn("parameter file_x0 not given, estimated primary impulse responses not saved");
		file_x0 = NULL;
	}
	if(!getparstring("file_wav", &file_wav)){
		if (verbose) sawarn("parameter file_wav not given, estimated wavelets not saved");
		file_wav = NULL;
	}
	
	if(!getparint("QC", &QC)) QC = 0;
	
	if(!getparstring("file_QC", &file_QC) && QC > 0 ){
		sawarn("parameter file_QC not given, assume pipe");
		file_QC = NULL;
		if (use_pipe == 1 && QC > 0 ) saerr("2 output file to the pipe: need to define file_prim or file_QC");
		
	}
	
	if(!getparstring("file_QC_dx0", &file_QC_dx0)){
		file_QC_dx0 = NULL;
	}		
		
	if(!getparint("n1", &optn)) optn = -1;
	if(!getparfloat("v0", &v0)) v0 = 1500;
	if(!getparfloat("t0", &t0)) t0 = 0.25;
	if(!getparfloat("tmax_rel", &tmax_rel)) tmax_rel = 3.3;
	if(!getparfloat("freq_base", &freq_base)) freq_base = 75.0;
	if(!getparfloat("Qfactor", &Qfactor)) Qfactor = 120;
	if(!getparfloat("size_twi", &size_twi)) size_twi = 0.1;
	if(!getparfloat("dtw", &dtw)) dtw = 0.048;	
	if(!getparint("ntmax", &ntmax)) ntmax = 1024;
	if(!getparint("nxmax", &nxmax)) nxmax = 512;
	if(!getparint("ntapx", &ntapx)) ntapx = 0;
	if(!getparint("ntapt", &ntapt)) ntapt = 21;	
	if(!getparint("lfilt", &lfilt)) lfilt = 41;
	if(!getparint("niter", &niter)) niter = 1;
	if(!getparint("nwindows", &nwindows)) nwindows = 1;
	if(!getparfloat("size_window", &size_window)) size_window = 0.6;
        if(!getparint("maxshots", &maxshots)) maxshots = 201;
/* parameters in case mode=1 */		
	if(!getparfloat("min_offset", &min_offset)) min_offset = 105;
	if(!getparfloat("tgap", &tgap)) tgap = 0.35;	
	
	if (mode==0) samess(" mode =0 : basic EPSI");
	if (mode==1) samess(" mode =1 : offsets up to %.2f m are for the first %.2f s reconstructed",min_offset, tgap);
/* Opening input file for reading */

	error = open_file(file_in, GUESS_TYPE, DELPHI_READ);
	if (error < 0 ) saerr("error in opening file %s", file_in);
	error = get_dims(file_in, &n1, &n2, &type);
	if (error >= 0) {
		if (!getparint("ntmax", &ntmax)) ntmax = n1;
		if (!getparint("nxmax", &nxmax)) nxmax = n2;
		if (verbose>=2 && (ntmax!=n1 || nxmax!=n2))
		    samess("dimensions overruled: %d x %d",ntmax,nxmax);
	}
	else {
		if (verbose>=2) samess("dimensions used: %d x %d",ntmax,nxmax);
		type = SA_TYPE_REAL;
	}
	if(!getparint("nrot", &nrot)) nrot=NINT(lfilt/2 - 0.5);
	//if(!getparint("ntapoffset", &ntapoffset)) ntapoffset=NINT((n2-1)/2);
	if(!getparint("ntapoffset", &ntapoffset)) ntapoffset=0;
	
	size = ntmax * nxmax;
	datin = (float *) malloc(size*sizeof(float));
	hdrs_in = (segyhdr *) malloc(nxmax*sizeof(segyhdr));
	hdrs_tot = (segyhdr *) malloc(n2*nxmax*sizeof(segyhdr));
	if (datin == NULL || hdrs_in==NULL )
		saerr("memory allocation error for input data");

	if (!getparstring("key", &key)) {
		ret = get_sort(file_in);
		if (ret < 0) key = "fldr";
		else key = getkey(ret);
	}
	if (verbose) samess("key used is %s",key);
	set_sukey(key);
	
	
	var_t0=1;
        if (!getparstring("key_depth", &key_depth)) {
	  samess("key_depth not defined, t0 = %.3f and v0 = %.3f are used to determine window", t0,v0);
	  var_t0=0;		
	}
	else {
	 samess("key_depth used is %s",key_depth);
	 ktype = hdtype(key_depth);
         kindex = getindex(key_depth);
	 
	 }	

/* Opening output file for writing */

	error = open_file(file_prim, GUESS_TYPE, DELPHI_CREATE);
	if (error < 0 ) saerr("error on creating output file for primaries");
	if (file_res != NULL) {
		error = open_file(file_res, GUESS_TYPE, DELPHI_CREATE);
		if (error < 0 ) saerr("error on creating output file for residual");
	}
	if (file_p_recon != NULL) {
		error = open_file(file_p_recon, GUESS_TYPE, DELPHI_CREATE);
		if (error < 0 ) saerr("error on creating output file for reconstructed data");
	}
	if (file_mult != NULL) {
		error = open_file(file_mult, GUESS_TYPE, DELPHI_CREATE);
		if (error < 0 ) saerr("error on creating output file for multiples");
	}
	if (file_x0 != NULL) {
		error = open_file(file_x0, GUESS_TYPE, DELPHI_CREATE);
		if (error < 0 ) saerr("error on creating output file for x0");
	}
	if (file_wav != NULL) {
		error = open_file(file_wav, GUESS_TYPE, DELPHI_CREATE);
		if (error < 0 ) saerr("error on creating output file for wavelets");
	}
	if (QC > 0) {
		error = open_file(file_QC, GUESS_TYPE, DELPHI_CREATE);
		if (error < 0 ) saerr("error on creating output file for QC");
	}
	if (file_QC_dx0 != NULL) {
		error = open_file(file_QC_dx0, GUESS_TYPE, DELPHI_CREATE);
		if (error < 0 ) saerr("error on creating output file for QC_dx0");
	}	

/* start reading data from input */

	error = 0;
	first = 1;
	ishot=0;

	for (ikey = 0; ikey < MAX_KEYS; ikey++)
	   	keys[ikey] = (char *) malloc(MAX_KEY_LENGTH);
				
		error = read_data(file_in, datin, size, &n1, &n2, &f1, &f2, 
			&d1, &d2, &type, hdrs_in);

		if (error < 0 ) {
			ret = close_file(file_in);
			if (ret < 0) sawarn("err %d on closing input file",ret);
			ret = close_file(file_prim);
			if (ret < 0) sawarn("err %d on closing output file",ret);
			if (file_res != NULL) {
				ret = close_file(file_res);
				if (ret < 0) sawarn("err %d on closing output file",ret);
			}
			if (file_p_recon != NULL) {
				ret = close_file(file_p_recon);
				if (ret < 0) sawarn("err %d on closing output file",ret);
			}
			if (file_mult != NULL) {
				ret = close_file(file_mult);
				if (ret < 0) sawarn("err %d on closing output file",ret);
			}
			if (file_x0 != NULL) {
				ret = close_file(file_x0);
				if (ret < 0) sawarn("err %d on closing output file",ret);
			}
			if (file_wav != NULL) {
				ret = close_file(file_wav);
				if (ret < 0) sawarn("err %d on closing output file",ret);
			}
			if (QC > 0) {
				ret = close_file(file_QC);
				if (ret < 0) sawarn("err %d on closing output file",ret);
			}						
			if (file_QC_dx0 != NULL) {
				ret = close_file(file_QC_dx0);
				if (ret < 0) sawarn("err %d on closing output file",ret);
			}			
			if (verbose) samess("error reading first gather");
			free(hdrs_in);
			free(datin);
			for (ikey=0; ikey<MAX_KEYS; ikey++) free(keys[ikey]);
			t2 = cputime();
			exit(0);
		}	
		
	 if (verbose) disp_info(file_in,n1,n2,f1,f2,d1,d2,type);
	 
/*  create an offset dependent taper */	 
	 
	offsettaper=(float *) malloc(n2*sizeof(float));
	if(offsettaper == NULL) saerr("memory allocation error offsettaper");
	 
	for (i = 0; i < n2; i++) offsettaper[i] = 1.0;
	if (ntapoffset) {
		for (i = 0; i < ntapoffset; i++) {
			offsettaper[n2-1-i] = cos(PI*(i-ntapoffset)/(ntapoffset*2))*cos(PI*(i-ntapoffset)/(ntapoffset*2)) ;
		}
	}	 
	 
	 
/* Create a taper in both the time- and the x axes */

        xtaper = (float *) malloc(n2*sizeof(float));
	if(xtaper == NULL) saerr("memory allocation error xtaper");

        for (i = 0; i < n2; i++) xtaper[i] = 1.0;
	if (ntapx) {
		for (i = 0; i < ntapx; i++) {
			xtaper[i] = (cos(PI*(i-ntapx)/ntapx)+1)/2.0;
			xtaper[n2-1-i] = xtaper[i];
		}
	}	 
	 
        ttaper = (float *) malloc(n1*sizeof(float));
	if(ttaper == NULL) saerr("memory allocation error ttaper");
	t0mask = (float *) malloc(n2*sizeof(float));
	if(t0mask == NULL) saerr("memory allocation error t0mask");
	

        for (i = 0; i < n1; i++) ttaper[i] = 1.0;
	if (ntapt) {
		for (i = 0; i < ntapx; i++) {
			ttaper[i] = (cos(PI*(i-ntapt)/ntapt)+1)/2.0;
			ttaper[n1-1-i] = ttaper[i];
		}
	} 
	
/* Define a gain factor sqrt(t) which will be used for the picking process */	

	wortelT = (float *) malloc(n1*sizeof(float));
	if(wortelT == NULL) saerr("memory allocation error wortelT");

        for (i = 0; i < n1; i++) wortelT[i] = sqrt(i*d1);
	
/*                               */	

        nsamp_window=NINT(size_window/d1);		 
	 	 		         	 
	while (error >= 0) {
		
		if (!getparint("n1", &optn)) optn = -1;
		if (optn < 0 ) optn = n1;
		n1_orig = n1;
/* Determine from axis information which transform is desired */
		get_axis(&dom1, &dom2);

		nkeys = get_keys(keys,seqnr);
		if (nkeys < 0) {
			nkeys = 1;
			seqnr[0] = 1;
			keys[0] = SA_CSG;
		}

		t1 = cputime();
		//if (first) if (verbose) samess("CPU-time input = %.3f", t1-tstart);
/*save all input data in one file*/
		if (first) {
			p = (float *) malloc(n1*n2*n2*sizeof(float));
			if (p == NULL) saerr("memory allocation error p");
			mask = (int *) malloc(n2*n2*sizeof(int));
			if (mask == NULL) saerr("memory allocation error mask");
		} 
		
		if (var_t0 == 1) {
			gethval((segy *) &hdrs_in[0], kindex, &kvalkey);
			depth = vtof(ktype, kvalkey);			
			t0mask[ishot]=depth*2/v0;
			dt0=nrot*d1;
			samp_max=MIN(NINT(tmax_rel*depth*2/(v0*d1)),n1);
			//samp_max=n1;
		} else {
			t0mask[ishot]=t0;
			dt0=nrot*d1;
			samp_max=MIN(NINT(tmax_rel*t0/d1),n1);
			//samp_max=n1;
		}
		
				   		
		for(i = 0; i < n2; i++) {
			mask[ishot*n2+i]=MIN(MAX(1,NINT(-dt0+sqrt(t0mask[ishot]*t0mask[ishot]+(hdrs_in[i].offset/v0)*(hdrs_in[i].offset/v0))/d1)),samp_max);
		
			tmp=xtaper[i]*offsettaper[abs(i-ishot)];

			hdrs_tot[ishot*n2+i]=hdrs_in[i];
		
			for(j = 0; j < n1; j++) {
				p[ishot*size+i*n1+j]=datin[i*n1+j]*tmp*ttaper[j];

			}  		
		}
											
	        ishot+=1;
		first=0;
			
	        error = read_data(file_in, datin, size, &n1, &n2, &f1, &f2, &d1, &d2, &type, hdrs_in);		

		if (error < 0 ) {
		  ret = close_file(file_in);
		  if (ret < 0) sawarn("err %d on closing input file",ret);
			 if (verbose) samess("end of data reached");
		 }		 
	}	 		
        free(datin);		
        t2 = cputime();

        nshots=ishot;
	nreceivers=n2;
	if (nshots!=n2) saerr("number of shots not equal to number of receivers: %d, %d", nshots, n2);
	
	if (verbose) samess("number of shots read in = %d", nshots);
		
	nQC_shots = countparval("QC_shots");
	if (!nQC_shots > 0) nQC_shots=1;
	QC_shots = ealloc1int(nQC_shots);
	if (!getparint("QC_shots",QC_shots)) QC_shots[0] = NINT(nshots/2);
	
/* In the first step PPcH is determined in the frequency domain and inverse FFT'd to time domain. Same for Pc/pc etc.. */	
        optn = 2*optncr(optn);
	nfreq = optn/2+1;
	if (verbose) samess("nfreq equals %d, optn = %d ", nfreq, optn);	
/* make near offsets zero and save original near offsets in no_orig  */	

	if (mode==1) {
		
	nsamp_gap=NINT(tgap/d1);
	Ngapl=NINT(min_offset/d2 -1);
	Ngap=2*Ngapl + 1;
	Ngapr=Ngap-Ngapl;	
	
	no_orig = (float *) malloc(nshots*Ngap*nsamp_gap*sizeof(float));
	if(no_orig == NULL) saerr("memory allocation error no_orig");	

	for(ishot = 0; ishot < nshots; ishot++) {
	  for(ireceiver = ishot-Ngapl; ireceiver <= ishot+Ngapl; ireceiver++) {
	    ix=ishot*Ngap*nsamp_gap+(Ngapl-ishot+ireceiver)*nsamp_gap;
	    ip=ishot*size+ireceiver*n1;
	    if (ireceiver  < 0 || ireceiver >= nreceivers ) {
	      for(isample = 0; isample < nsamp_gap; isample++) {	     
	        no_orig[ix+isample]=0;
	      }	      
	    } 	      
	    else {
	      for(isample = 0; isample < nsamp_gap; isample++) {
	        no_orig[ix+isample]=p[ip+isample];
	        p[ip+isample]=0;
	      }
	    }	      
	  }		  		
	}
		
/* create estimated near-offsets in time-domain: defined for Nt samples but only untill nsamp_gap will be upd
ated. Rest is needed for FFT to freq. domain. */

	pin = (float *) malloc(nshots*Ngap*nsamp_gap*sizeof(float));
	if(pin == NULL) saerr("memory allocation error pin");
	
	pinc = (float *) malloc(nshots*Ngap*nsamp_gap*sizeof(float));
	if(pinc == NULL) saerr("memory allocation error pinc");	
			
	tmp_dpin = (float *) malloc(Ngap*nsamp_gap*sizeof(float));
	if(tmp_dpin == NULL) saerr("memory allocation error tmp_dpin");
	
	dpin = (float *) malloc(nshots*Ngap*nsamp_gap*sizeof(float));
	if(dpin == NULL) saerr("memory allocation error dpin");
	
	dpinc = (float *) malloc(nshots*Ngap*nsamp_gap*sizeof(float));
	if(dpinc == NULL) saerr("memory allocation error dpinc");
	
	rdata_pin = (float *) malloc(optn*Ngap*sizeof(float));
	if(rdata_pin == NULL) saerr("memory allocation error rdata_pin");
			
	tmpdata_pin = (float *) malloc(2*nfreq*Ngap*sizeof(float));
	if(tmpdata_pin == NULL) saerr("memory allocation error tmpdata_pin");
		
        datout_pinc = (float *) malloc(2*optn*Ngap*sizeof(float));
	if(datout_pinc == NULL) saerr("memory allocation error datout_pinc");		
			
	cdata_pinc = (complex *) malloc(nfreq*Ngap*sizeof(complex));
	if(cdata_pinc == NULL) saerr("memory allocation error cdata_pinc");		
	
	mask_beta = ealloc2float(n2,nshots);   
	if(mask_beta == NULL) saerr("memory allocation error mask_beta");
	 
	 
	 for(i = 0; i < Ngap*nsamp_gap; i++) {	  	     
	        tmp_dpin[i]=0;	 	      		  		
	 }
	 
	 for(i = 0; i < nshots*Ngap*nsamp_gap; i++) {
	        pin[i]=0;
	 	pinc[i]=0;  	     
	        dpin[i]=0;
		dpinc[i]=0;	 	      		  		
	 }
	 	
	} /* end of mode==1*/
	
/* 2D and dipole correction is applied to the data in the freq domain */
/* correction: k*i*exp(-pi/4*i) = exp(pi/4)*sqrt(k) = sqrt(k)*(cos(pi/4)+i*sin(pi/4)) */
/* cos(pi/4)=sin(pi/4)= 0.5*sqrt(2)*/	 
	rdata_wav = (float *) malloc(optn*1*sizeof(float));
	if(rdata_wav == NULL) saerr("memory allocation error rdata_wav");
			
	tmpdata_wav = (float *) malloc(2*nfreq*1*sizeof(float));
	if(tmpdata_wav == NULL) saerr("memory allocation error tmpdata_wav");
		
        datout_wav = (float *) malloc(2*optn*1*sizeof(float));
	if(datout_wav == NULL) saerr("memory allocation error datout_wav");		
			
	cdata_wav = (complex *) malloc(nfreq*1*sizeof(complex));
	if(cdata_wav == NULL) saerr("memory allocation error cdata_wav");
	
			
	 	
        corr2D = (float *) malloc(nfreq*sizeof(float));
	if(corr2D == NULL) saerr("memory allocation error corr2D");
	rdata = (float *) malloc(optn*n2*sizeof(float));
	if(rdata == NULL) saerr("memory allocation error rdata");
	tmpdata = (float *) malloc(2*nfreq*n2*sizeof(float));
	if(tmpdata == NULL) saerr("memory allocation error tmpdata");	
	Preal = (float *) malloc(nfreq*n2*nshots*sizeof(float));
	if(Preal == NULL) saerr("memory allocation error Preal");
	Pim = (float *) malloc(nfreq*n2*nshots*sizeof(float));
	if(Pim == NULL) saerr("memory allocation error Pim");
	Pcreal = (float *) malloc(nfreq*n2*nshots*sizeof(float));
	if(Pcreal == NULL) saerr("memory allocation error Pcreal");
	Pcim = (float *) malloc(nfreq*n2*nshots*sizeof(float));
	if(Pcim == NULL) saerr("memory allocation error Pcim");	
	
        for(i = 0; i < nfreq; i++) {
	  corr2D[i] = 0.5*sqrt(2)*sqrt(i);		
	} 
				
	ishot=0;	
					
	while (ishot<nshots) {	

		sign = -1;
		scale = 1.0;
		for(i = 0; i < n2; i++) {
		  nmask=mask[ishot*n2+i];
		  for(j = 0; j < nmask; j++) rdata[i*optn+j] = 0.0;
		  for(j = nmask; j < n1; j++) rdata[i*optn+j] = p[ishot*size+i*n1+j]*scale;
		  for(j = n1; j < optn; j++) rdata[i*optn+j] = 0.0;
		}

		rcmfft(&rdata[0], (complex *)&tmpdata[0], optn, n2, optn, nfreq, sign);
			
/* save all FFT-data in one file, give number of shots on beforehand */	
								
		for(i = 0; i < nreceivers; i++) {
		  for(j = 0; j < nfreq; j++) {
		            Preal[j*nreceivers*nshots+i*nshots+ishot] = tmpdata[i*nfreq*2+j*2];
			    Pim[j*nreceivers*nshots+i*nshots+ishot] = tmpdata[i*nfreq*2+j*2+1];
			    Pcreal[j*nreceivers*nshots+i*nshots+ishot] = corr2D[j]*(tmpdata[i*nfreq*2+j*2]-tmpdata[i*nfreq*2+j*2+1]);
			    Pcim[j*nreceivers*nshots+i*nshots+ishot] =corr2D[j]*(tmpdata[i*nfreq*2+j*2]+tmpdata[i*nfreq*2+j*2+1]);
			    
		  }
		} 		  			            			
	      ishot+=1;							
	}	
	
          t2 = cputime();
	  //samess("Total CPU-time after getting P and Pc in frequency domain = %.3f", t2-tstart);
	    
	  
/*   P and Pc are determined, now we can get P*PcH   */

	PPcH = (float *) malloc(2*nfreq*n2*nshots*sizeof(float));
	if(PPcH == NULL) saerr("memory allocation error PPcH");
			
/*  correlation  */
	tmp1=0;
	tmp2=0;	
        for(ifreq = 0; ifreq < nfreq; ifreq++) {
	    index_freq_start=ifreq*nreceivers*nshots;
	  for (ishot = 0; ishot < nshots; ishot++) {
	    index_out_shot=2*ishot*nreceivers*nfreq+2*ifreq;
	    index_shot_start=index_freq_start+ishot*nshots;
	    for (ireceiver = 0; ireceiver < nreceivers; ireceiver++) {
	      index_receiver_start=index_freq_start+ireceiver*nshots;
	      index_out=index_out_shot+2*ireceiver*nfreq;
	     for (ix = 0; ix < nshots; ix++) {	                
		tmp1+=Preal[index_receiver_start+ix]*Pcreal[index_shot_start+ix]+Pim[index_receiver_start+ix]*Pcim[index_shot_start+ix];
	        tmp2+=Pim[index_receiver_start+ix]*Pcreal[index_shot_start+ix]-Preal[index_receiver_start+ix]*Pcim[index_shot_start+ix];							              
	      }
	      PPcH[index_out]=tmp1;
	      PPcH[index_out+1]=tmp2;
	      tmp1=0;
	      tmp2=0;	      
	    }
	  }
        }   
	
	t3 = cputime();
	//samess("CPU-time to apply correlation in frequency domain = %.3f", t3-t2);
	
/*  inverse FFT to get ppcH = dx0(1) and pc ; so we get back to the time-domain  */

	cdata = (complex *) malloc(nfreq*n2*sizeof(complex));
	if(cdata == NULL) saerr("memory allocation error cdata");		
	datout = (float *) malloc(2*optn*n2*sizeof(float));
	if(datout == NULL) saerr("memory allocation error datout");
	if (file_mult != NULL || QC == 3 || mode==1) {		
	  pc = (float *) malloc(nshots*n2*n1*sizeof(float));
	  if(pc == NULL) saerr("memory allocation error pc");
	}
	dx0 = (float *) malloc(nshots*n2*n1*sizeof(float));
	if(dx0 == NULL) saerr("memory allocation error dx0");				
        B = calloc2float(n1,nshots*n2);
	if(B == NULL) saerr("memory allocation error B");	
	
        for (i = 0; i < nshots*n2; i++) {
	  for (j = 0; j < n1; j++) {   
	   B[i][j]=p[i*n1+j];	  
	  }
	 }
	 
        ishot=0;
			
	while (ishot<nshots) {		
		
		scale = 1.0/optn;
		
		for(i = 0; i < n2; i++) {
			for(j = 0; j < nfreq; j++) {
			  cdata[i*nfreq+j].r = PPcH[2*ishot*n2*nfreq+2*i*nfreq+2*j]*scale;
			  cdata[i*nfreq+j].i = PPcH[2*ishot*n2*nfreq+2*i*nfreq+2*j+1]*scale;
			}
		}		  
		  
                 sign = 1;

		 crmfft(&cdata[0], &datout[0], optn, n2, nfreq, optn, sign);
			
/* copy data to original trace length and apply mask: following is to take the positive times of output. First get ppcH which is inmediatlly called D */
              
                 if (n1_orig>0 && n1_orig<optn) {
				
		     for (i = 0; i < n2; i++) {
			  for (j = 0; j < n1_orig; j++) {
			    dx0[ishot*n2*n1+i*n1+j]=-datout[i*optn+j];			    
			  }	 				
		       }
		  }
			
/*inverse FFT Pc in the same loop to pc  */ 			
			
		for(i = 0; i < n2; i++) {
			for(j = 0; j < nfreq; j++) {
			  cdata[i*nfreq+j].r = Pcreal[j*n2*nshots+i*nshots+ishot]*scale;
			  cdata[i*nfreq+j].i = Pcim[j*n2*nshots+i*nshots+ishot]*scale;
			}
		  }		  
		  
                 sign = 1;

	         crmfft(&cdata[0], &datout[0], optn, n2, nfreq, optn, sign);
			
/* copy data to original trace length (and apply mask if desired): */

                 if (n1_orig>0 && n1_orig<optn) {
		   if (file_mult != NULL || QC > 2 || mode==1) {		
		     for (i = 0; i < n2; i++) {                      
		       for (j = 0; j < n1_orig; j++) {
			pc[ishot*n2*n1+i*n1+j]=datout[i*optn+j];
			}	 				
		     }
		    }
		    else {
		     for (i = 0; i < n2; i++) {                      
		       for (j = 0; j < n1_orig; j++) {
			p[ishot*n2*n1+i*n1+j]=datout[i*optn+j];
			}	 				
		      }		    
		    } 
		  }	  
		ishot+=1;								
		}
	    	
          t2 = cputime();
	  //samess("CPU-time for inverse FFT to get dx0 and pc n time domain = %.3f", t2-t3);	
	  
/*  in first iteration dx0 = ppcH; copy ppcH to dx0 abd save ppcH as D = ppcH+x0*pcpcH. D can be updated by adding dx0*pcpcH in each iteration.
Same is valid for B = p + x0pc. Note that dx0 is not yet sparse here. After picking the highest values you get the sparse update dx0 which is saved in different tables   The PPcH and Pcreal and Pcim can be freed  */

        if (file_mult == NULL && QC < 3 && mode==0) {
	 pc=p;
	}	 
			
	x0 = calloc2float(n1,nshots*n2);
	if(x0 == NULL) saerr("memory allocation error x0");
	
	dpminds = calloc2float(n1,nshots*n2);
	if(dpminds == NULL) saerr("memory allocation error dpminds");			
		
	mask_s = ealloc2float(n1,n2);
	if(mask_s == NULL) saerr("memory allocation error mask_s");
	
	mask_alfa = calloc2float(n2,nshots);
	if(mask_alfa == NULL) saerr("memory allocation error mask_alfa");
			
	x0s = calloc2float(n1,nshots*n2);
	if(x0s == NULL) saerr("memory allocation error x0s");
	
	//s = (float *) malloc(lfilt*sizeof(float));
	s = calloc2float(lfilt,nwindows);
	if(s == NULL) saerr("memory allocation error s");
	
	//s_rev = (float *) malloc(lfilt*sizeof(float));
	s_rev = calloc2float(lfilt,nwindows);
	if(s_rev == NULL) saerr("memory allocation error s_rev");	
	
	s_tmp = (float *) malloc(lfilt*sizeof(float));
	if(s_tmp == NULL) saerr("memory allocation error s_tmp");

	tmp_recgather = (float *) malloc(nshots*n1*n2*sizeof(float));
	if(tmp_recgather == NULL) saerr("memory allocation error tmp_recgather");
	
	tmp_B = calloc2float(n1,n2);
	if(tmp_B == NULL) saerr("memory allocation error tmp_B");	
	
	tmp_x0s = calloc2float(n1,n2);
	if(tmp_x0s == NULL) saerr("memory allocation error tmp_x0s");		
		
        for (i = 0; i < nshots*n2; i++) {
	  for (j = 0; j < n1; j++) {
	   x0s[i][j]=0;	   	  
	  }
	 }
	
	/* mask_alfa can be used to leave out shots/receivers at the edge of the data during the matching of dpminds to the residual */
	       
        if (mode==1) {
          for (i = 0; i < nshots; i++) {
	    for (j = 0; j < n2; j++) {
	     if (j < i+Ngapl && j > i-Ngapl ) {	  
	      mask_alfa[i][j]=0;	
	      mask_beta[i][j]=1;   	   
	    }
	     else  {	  
	      mask_alfa[i][j]=1;	
	      mask_beta[i][j]=1;   	   
	    }	  
	   }
	  } 
	}
	if (mode==0) {
          for (i = 0; i < nshots; i++) {
	    for (j = 0; j < n2; j++) {	  
	      mask_alfa[i][j]=1;  	   	  
	    }
	  }		
	}
			       
	for (i = 0; i < nshots; i++) {
	  for (j = 0; j < n1; j++) {
	   tmp_recgather[i*n1+j]=0;
	  }
	 } 
 
/*           allocate some more           */	

        Npeak=4;

	dx0_shot= (int *) malloc(Npeak*nshots*n2*sizeof(int));
	if(dx0_shot == NULL) saerr("memory allocation error dx0_shot");
	dx0_receiver= (int *) malloc(Npeak*nshots*n2*sizeof(int));
	if(dx0_receiver == NULL) saerr("memory allocation error dx0_receiver");
	dx0_sample= (int *) malloc(Npeak*nshots*n2*sizeof(int));
	if(dx0_sample == NULL) saerr("memory allocation error dx0_sample");
	dx0_amp= (float *) malloc(Npeak*nshots*n2*sizeof(float));
	if(dx0_amp == NULL) saerr("memory allocation error dx0_amp");	    

/*           start iterations            */

      iter=0; 
      first=1;
      Dp=NINT(0.04/d1);   /*  number of samples mutes around maximum peak  */     
      dstwi = NINT(dtw/d1);
      nstwi = NINT(size_twi/d1);

      while (iter < niter) {	
      
        t2 = cputime();      
            
	 /*   QC after every iteration the x0s, multiples and residual    */
	 
	 for (iQC = 0; iQC < nQC_shots; iQC++) {
	   

	   QC_shot=QC_shots[iQC]-1;
           f2=(float) hdrs_tot[QC_shot*n2].offset;					 

		if (QC > 0) {
				    
		        dom1 = SA_AXIS_TIME;
			type = SA_TYPE_REAL;
			trid = (short) axis2trid(dom1, dom2, type);			


	     for(i = 0; i < n2; i++) {
			hdrs_in[i]=hdrs_tot[QC_shot*n2+i];
			hdrs_in[i].duse = iter;	
		}
		
		if (mode==1) {
			for (i = 0; i < n2; i++) {
			  for (j = 0; j < n1; j++) {
                            datout[i*n1_orig+j]=p[QC_shot*n2*n1+i*n1+j];
			   }	 				
			}					
									
		  ret = write_data(file_QC,datout,n1,n2,f1,f2,d1,d2,type,hdrs_in);
		  if (ret < 0 ) saerr("error on writing output file_QC.");
		}
		
			for (i = 0; i < n2; i++) {
			  for (j = 0; j < n1; j++) {
                            datout[i*n1_orig+j]=x0s[QC_shot*n2+i][j];
			   }					
			}
				
		ret = write_data(file_QC,datout,n1,n2,f1,f2,d1,d2,type,hdrs_in);
		if (ret < 0 ) saerr("error on writing output file_QC.");		
					
		}
				
		if (QC > 1) {
			for (i = 0; i < n2; i++) {
			  for (j = 0; j < n1; j++) {
                            datout[i*n1_orig+j]=B[QC_shot*n2+i][j]-x0s[QC_shot*n2+i][j];
			   }	 				
			}	
		ret = write_data(file_QC,datout,n1,n2,f1,f2,d1,d2,type,hdrs_in);
		if (ret < 0 ) saerr("error on writing output file_QC.");	
		
		}
						
		if (QC > 2) {
		
			for (i = 0; i < n2; i++) {
			  for (j = 0; j < n1; j++) {
                            datout[i*n1_orig+j]=-B[QC_shot*n2+i][j]+p[QC_shot*n2*n1+i*n1+j];
			   }					
			}	
		ret = write_data(file_QC,datout,n1,n2,f1,f2,d1,d2,type,hdrs_in);
		if (ret < 0 ) saerr("error on writing output file_QC.");	
		}
		
/*  write dx0 to file_QC_dx0 for the defined QC-shots */
		
                if (file_QC_dx0 != NULL) {
		
		   if (QC==0) {
		        dom1 = SA_AXIS_TIME;
			type = SA_TYPE_REAL;
			trid = (short) axis2trid(dom1, dom2, type);			
		     }		
		  
    	          for(i = 0; i < n2; i++) {

			hdrs_in[i]=hdrs_tot[QC_shot*n2+i];
			hdrs_in[i].duse = iter+1;	
		    }		
			 	   	     
			for (i = 0; i < n2; i++) {
			  for (j = 0; j < n1; j++) {
                            datout[i*n1_orig+j]=wortelT[j]*dx0[QC_shot*n2*n1+i*n1+j];
			   }	 				
			}
		
		  ret = write_data(file_QC_dx0,datout,n1,n2,f1,f2,d1,d2,type,hdrs_in);
		  if (ret < 0 ) saerr("error on writing output file_QC_dx0.");
	         } 						
	  }      
	  
/*  pick the highest peaks from every trace of x0. Npeak could be set within the loop such that it can be made iteration dependant if desired. */
/*  In that case a npeak_max should be set to allocate the table sizes of dx0_shots, dx0_receiver etc..*/
	 	
        for (i = 0; i < nshots*n2; i++) {
	  for (j = 0; j < n1; j++) {
	   dpminds[i][j]=0;
	   }
	 }	 
	   
	if (iter <= 3) Npeak=1;
	if (iter > 3) Npeak=2;
	if (iter >= 6) Npeak=4;  
	  	          
          ispike=0; 
	   
          for (ishot = 0; ishot < nshots; ishot++) { 
	  
	    mstwi = mask[ishot*n2+ishot]+nstwi;
	    	    
	    ireceiver_start=MAX(ishot-25-iter*7,0); 
	    ireceiver_end=MIN(ishot+25+iter*7,n2);   			
	    for (ireceiver = ireceiver_start; ireceiver < ireceiver_end; ireceiver++) {	
	      index_receiver_start=ishot*n2*n1+ireceiver*n1;
	      index_sample_start=mask[ishot*n2+ireceiver];
	      for (ipeak = 0; ipeak < Npeak; ipeak++) {
	        max=0;
	      	for (isample = index_sample_start; isample < MIN(mstwi + iter*dstwi,n1-lfilt); isample++) {
		  if (fabs(wortelT[isample]*dx0[index_receiver_start+isample]) > fabs(max)) {
		   max=wortelT[isample]*dx0[index_receiver_start+isample];
		   sample_max=isample;
		   }
	        }
		   dx0_shot[ispike]=ishot;
		   dx0_receiver[ispike]=ireceiver;
		   dx0_sample[ispike]=sample_max;
		   dx0_amp[ispike]=max/wortelT[sample_max];
		   ispike++;
		   if (Npeak > 1) {
		     for (i = -Dp; i < Dp+1; i++) {
		      dx0[index_receiver_start+sample_max+i]=0;
		     }
		   }				
	      }
	      	      
	    }
	  }
	  	 
          nspikes=ispike;
	  if (verbose) samess("number of spikes picked from dx0 = %d ", nspikes);  
	        
	  t3 = cputime();
	  //samess("CPU-time to write QC shots and apply the picking to dxo = %.3f", t3-t2);
	  	

/* first calculate dpminds before calculating alfa factor : x0 = alfa*dx0   */
/* in the first iteration step you only need to calculate the convolution of dx0 with pc (first try x0p before x0pc) because wavelet is still zero at this moment*/
      
     for (ireceiver = 0; ireceiver < n2; ireceiver++) {      
       for (ispike = 0; ispike < nspikes; ispike++) {
         if (dx0_receiver[ispike] == ireceiver) {

       
        /* within this loop first calcuate dx0s   */
       
         ampspk=dx0_amp[ispike];
	 sampspk=dx0_sample[ispike];
	 shotspk=dx0_shot[ispike];
	 recspk=dx0_receiver[ispike];	 
	 itmppc=shotspk*n1-sampspk;
	 itmp=shotspk*n1+sampspk-nrot;   
	 
	 iwav=0;
	 if (nwindows > 1) {
	   for (iwindow = 1; iwindow < nwindows; iwindow++) {
	    if (sampspk >= iwindow*nsamp_window && sampspk < (iwindow+1)*nsamp_window) iwav=iwindow;
	    //if (sampspk >= 300 && sampspk < 450) iwav=2;
	    //if (sampspk >= 450 && sampspk < 600) iwav=3;   
	   }
	   if (sampspk >= (nwindows-1)*nsamp_window) iwav=nwindows-1;
	 }
	 	
	
       /* first subtract -dx0*s  */	
	 if (first == 0) {
	   for (it = 0; it < lfilt; it++) {	     
	 tmp_recgather[itmp+it]+=-ampspk*s[iwav][it];
	    }	    	    	 
	 } 
	 	   
       /* now add dx0pc   */
           for (ishot = 0; ishot < nshots; ishot++) {	     
	     ix=ishot*n1;	    
	     istartpc=ishot*n2*n1+itmppc;
	     	     
	     for (isample = sampspk; isample < n1; isample++) {
	     tmp_recgather[ix+isample]+=ampspk*pc[istartpc+isample];
	     }	      	 
	   }
	 }        
       }
      
  	for (i = 0; i < nshots; i++) {
	  for (j = 0; j < n1; j++) {
	   dpminds[i*n2+ireceiver][j]+=tmp_recgather[i*n1+j];
	   tmp_recgather[i*n1+j]=0;	   
	  }
	 }              
     }  
       
   t2 = cputime();
   //samess("CPU-time to determine dx0pc - dx0s = %.3f", t2-t3);	
   
        
/*  calculate alfa to update x0=x0+alfa*dx0 */
       
       tmp1=0;
       tmp2=0;
       for (ishot = 0; ishot < nshots; ishot++) { 
         for (ireceiver = 0; ireceiver < n2; ireceiver++) {
          for (j = 0; j < n1; j++) {
           tmp1+=mask_alfa[ishot][ireceiver]*(dpminds[ishot*n2+ireceiver][j])*(B[ishot*n2+ireceiver][j]-x0s[ishot*n2+ireceiver][j]);
	   tmp2+=mask_alfa[ishot][ireceiver]*(dpminds[ishot*n2+ireceiver][j])*(dpminds[ishot*n2+ireceiver][j]);	 
	  }  
         } 
       }   
       
       alfa=-tmp1/tmp2;
       /*if (verbose ) {
       samess("scaling factor alfa for updating x0 = %.15f", alfa);   
       } */
                     
   t3 = cputime();
          
/*  having alfa we can update x0 = x0 + alfa*dx0, also update B = B + alfa*dx0*pc, s will be estimated by minimizing B - x0*s,   and finally also D = ppcH+x0*pcpcH                 dx0 = residu * sH + x0s*pcH           */       
 /* to save memory space dpminds is first added (dx0s should not be included) to B and in a second step will be compensated for the extra dx0s (this way we don't have to keep dx0s in the memory) */
  
	for (i = 0; i < nshots*n2; i++) {
	  for (j = 0; j < n1; j++) {
	   B[i][j]+=alfa*dpminds[i][j];
	  }
	 }   
          
       tmp1=0;   
       for (ispike = 0; ispike < nspikes; ispike++) {
               
         ampspk=alfa*dx0_amp[ispike];
	 sampspk=dx0_sample[ispike];
	 shotspk=dx0_shot[ispike];
	 recspk=dx0_receiver[ispike];
	 itmp=shotspk*n2+recspk;
	 
	/* if (nwindows > 1) {
	   if (sampspk < 150) iwav=0;
	   if (sampspk >= 150 && sampspk < 300) iwav=1;
	   if (sampspk >= 300 && sampspk < 450) iwav=2;
	   if (sampspk >= 450 && sampspk < 600) iwav=3;
	   if (sampspk >= 600) iwav=4;
	 }        */      
        /* within this loop compensate B for the extra dx0s term   */
	
	 //iwav=0;	 
         if (nwindows > 1) {
	   for (iwindow = 0; iwindow < nwindows; iwindow++) {
	    if (sampspk >= iwindow*nsamp_window && sampspk < (iwindow+1)*nsamp_window) iwav=iwindow;  
	   }
	   if (sampspk >= (nwindows-1)*nsamp_window) iwav=nwindows-1;
          }
	 
	 
	 if (first == 0) {
	   for (it = 0; it < lfilt; it++) {	     
	     B[itmp][sampspk-nrot+it]+=ampspk*s[iwav][it];	     
	    }	    	    	 
	 }  	
          x0[itmp][sampspk]+=ampspk;                 
      }  

       t2 = cputime();
	
	
/* make the basis wavelet zero before estimating the updated basis wavelet   */	

	for (i = 0; i < lfilt; i++) {
	 s[0][i]=0;
	}
	
/*  Here s will be estimated by minimizing (B - x0*s)2 : matching x0 to B shot by shot and averaging the resulting s over the nuber of shots  */	 

tmp2_x0 = calloc2float(n1,n2);				
if(tmp2_x0 == NULL) saerr("memory allocation error tmp2_x0");  
	 
      for (ishot = 0; ishot <nshots; ishot++) {
             
         for (ireceiver = 0; ireceiver < n2; ireceiver++) {
	   if (mode==1 && ireceiver < ishot+Ngapl && ireceiver > ishot-Ngapl) {

              for (j = 0; j < n1; j++) {
	        tmp2_x0[ireceiver][j]=x0[ishot*n2+ireceiver][j];   
	        tmp_B[ireceiver][j]=B[ishot*n2+ireceiver][j];
	        mask_s[ireceiver][j]=0;	 
	       }
	   }	      
	   else {
	   	   
	    for (j = 0; j < n1; j++) {
	     tmp2_x0[ireceiver][j]=x0[ishot*n2+ireceiver][j];
	     tmp_B[ireceiver][j]=B[ishot*n2+ireceiver][j];	 
	      }
	    if (nwindows > 1) { 
	     for (j = 0; j < nsamp_window-10; j++) {	    
	      mask_s[ireceiver][j]=1;	      		       
	     } 
	    for (j = nsamp_window-10; j < nsamp_window+30; j++) {	
	      mask_s[ireceiver][j]=1-(j-(nsamp_window-10))*0.025;      	      	 
	     }  
	    for (j = nsamp_window+30; j < n1; j++) {
	      mask_s[ireceiver][j]=0;		       
	     }
	    }
	    else {
	    for (j = 0; j < n1; j++) {
	      mask_s[ireceiver][j]=1;		       
	     }
	    } 
	      	      	      	     	           
	   }
	   
	     
       }
                	 
	least_filterc(tmp2_x0, tmp_B, mask_s, n1, n2, lfilt, nrot, 0.01, 0.0000001, 1000, s_tmp, tmp_x0s); 
 
	 for (it=0;it<lfilt;it++) {
		s[0][it]+=s_tmp[it]/nshots;
	   } 	  
       }
       
       
       /***************** There is a bug when executing this line  ******************/
       free(tmp2_x0);
       /****************************************************************************/ 							
       
   
   /*bring wavelet to frequency domain and apply absorbtion correction to the basis wavelet s[0][it]  */ 
              if (nwindows > 1) {    
		sign = -1;
		scale = 1.0;
			
		for(j = 0; j < lfilt; j++) rdata_wav[j] = s[0][j]*scale;				
		for(j = lfilt; j < optn; j++) rdata_wav[j] = 0.0;
              
                rcmfft(&rdata_wav[0], (complex *)&tmpdata_wav[0], optn, 1, optn, nfreq, sign);

                scale = 1.0/optn;      

                for(iwav = 1; iwav < nwindows; iwav++) { 

		for(j = 0; j < nfreq; j++) {	
		
		 if (j > 0) { 
		  param=j/(2*d1*(nfreq-1)*freq_base);
		  tmp100 = 1-(220*log(param)/Qfactor);
		 }
		 else {
		  tmp100=1;
		 }

		 //cdata_wav[j].r = tmpdata_wav[j*2]*scale*(1-(220*tmp100/(480-iwav*90)));
		 cdata_wav[j].r = tmpdata_wav[j*2]*scale*(1+iwav*(tmp100-1)/(nwindows-1));
		 cdata_wav[j].i = tmpdata_wav[j*2+1]*scale;
		
		}
		
                sign = 1;

		crmfft(&cdata_wav[0], &datout_wav[0], optn, 1, nfreq, optn, sign);  
		 
		for(it = 0; it < lfilt; it++) s[iwav][it]=datout_wav[it];
		}
			 
	      }	
		       	
	
   /*    end absorption correction       */	
	 
	for (iwav=0;iwav<nwindows;iwav++)  {
          for (it=0;it<lfilt;it++) {
           s_rev[iwav][lfilt-1-it]=s[iwav][it];
	  }  
	}

   /*apply the averaged wavelet s to x0  */	 
	  	
	iwav=0;	
	for (ix=0;ix<nshots*n2;ix++) {		
	  for (itout=lfilt-1-nrot; itout < n1-nrot; itout++) {
	    x0s[ix][itout] = 0 ; 
	    itin=itout-lfilt+1+nrot;		   
	    for (it=0;it<lfilt;it++) {
	      if (nwindows > 1) {
	        for (iwindow = 0; iwindow < nwindows; iwindow++) {
	         if (itin >= iwindow*nsamp_window && itin < (iwindow+1)*nsamp_window) iwav=iwindow;  
	        }	      
                 if (itin >= (nwindows-1)*nsamp_window) iwav=nwindows-1;
	       }    	    
	       x0s[ix][itout] += (x0[ix][itin])*s_rev[iwav][it];
	       itin++;
	    }
	  }
	}
	  
       t3 = cputime();
       //samess("CPU-time to determine s and x0s = %.3f", t3-t2); 
 	 
  /* determine dpin in time-domain = -(I+x0c)H*residu; here first part: x0H*residu  */
  
  if (mode==1) {	 

	for (ishot = 0; ishot < nshots; ishot++) {
         v_shot_start=MAX(ishot-Ngapl,0);
	 v_shot_end=MIN(ishot+Ngapr,nshots);	  
          for (ireceiver = 0; ireceiver < nreceivers; ireceiver++) {	   
	   for (isample = 0; isample < n1; isample++) {
	     if (fabs(x0[ishot*nreceivers+ireceiver][isample]) > 0) {	      
	      ampspk=x0[ishot*nreceivers+ireceiver][isample];
	      for (ivshot = v_shot_start; ivshot < v_shot_end; ivshot++) {
	        for (it = isample; it < MIN(nsamp_gap+isample,n1); it++) {
	        tmp_dpin[(ivshot-ishot+Ngapl)*nsamp_gap-isample+it]+=ampspk*(B[ivshot*nreceivers+ireceiver][it]-x0s[ivshot*nreceivers+ireceiver][it]);	 		 		 		 
	        }		
	       }	       
	      }
	     }	     
	    }
	    
 /*  one limited receiver gather of dpin is now determined and can add it to dpin with - sign  */	    
	    for (ipshot = v_shot_start; ipshot < v_shot_end; ipshot++) {	    
	     for (ipsample = 0; ipsample < nsamp_gap; ipsample++) {
	       dpin[ipshot*Ngap*nsamp_gap+(ishot-ipshot+Ngapl)*nsamp_gap+ipsample]=-tmp_dpin[(ipshot-ishot+Ngapl)*nsamp_gap+ipsample];
	              	       
	tmp_dpin[(ipshot-ishot+Ngapl)*nsamp_gap+ipsample]=0;       
	      }
	    }
	   }
       t2 = cputime();
       //if (verbose) samess("CPU-time to determine dpin=x0' * V = %.3f", t2-t3); 	   
	   	   
/*               apply 2D correction on -x0'*V                 */
	 
  	ishot=0;	
					
	while (ishot<nshots) {	

		sign = -1;
		scale = 1.0;
			
		for(ireceiver = 0; ireceiver < Ngap; ireceiver++) {
			for(j = 0; j < nsamp_gap; j++) 
				rdata_pin[ireceiver*optn+j] = dpin[ishot*Ngap*nsamp_gap+ireceiver*nsamp_gap+j]*scale;				
			for(j = nsamp_gap; j < optn; j++) 
				rdata_pin[ireceiver*optn+j] = 0.0;
		}
              
                rcmfft(&rdata_pin[0], (complex *)&tmpdata_pin[0], optn, Ngap, optn, nfreq, sign);

             scale = 1.0/optn;
		
		for(ireceiver = 0; ireceiver < Ngap; ireceiver++) {
			for(j = 0; j < nfreq; j++) {
			  cdata_pinc[ireceiver*nfreq+j].r = corr2D[j]*(tmpdata_pin[ireceiver*nfreq*2+j*2]+tmpdata_pin[ireceiver*nfreq*2+j*2+1])*scale;
			  cdata_pinc[ireceiver*nfreq+j].i = corr2D[j]*(-tmpdata_pin[ireceiver*nfreq*2+j*2]+tmpdata_pin[ireceiver*nfreq*2+j*2+1])*scale;
			}
		}		  
		
                 sign = 1;

		 crmfft(&cdata_pinc[0], &datout_pinc[0], optn, Ngap, nfreq, optn, sign);
		 		 			              
                 if (nsamp_gap>0 && nsamp_gap<optn) {				
		     for (ireceiver = 0; ireceiver < Ngap; ireceiver++) {
		       itmp=ishot*n2+(ireceiver+ishot-Ngapl);
		       if (itmp < 0 || itmp > n2*n2-1 ) {
		         for (j = 0; j < nsamp_gap; j++) {
			   dpin[ishot*Ngap*nsamp_gap+ireceiver*nsamp_gap+j]=0;
		           }		       
		         }
		       else {
		        nmask=mask[itmp];
			for (j = 0; j < nmask; j++) {
			   dpin[ishot*Ngap*nsamp_gap+ireceiver*nsamp_gap+j]=0;
		         }			
			for (j = nmask; j < nsamp_gap; j++) {
			   dpin[ishot*Ngap*nsamp_gap+ireceiver*nsamp_gap+j]=datout_pinc[ireceiver*optn+j];
		         }
		       }	 	 				
		     }
		 }				
			  			            			
	        ishot+=1;							
	}       
 
       t3 = cputime();	  	      
/* add -residu  to complete dpin */
  		   
  	for (ishot = 0; ishot < nshots; ishot++) { 
	  v_shot_start=MAX(ishot-Ngapl,0);
	  v_shot_end=MIN(ishot+Ngapr,nshots);
	  for (ipshot = v_shot_start; ipshot < v_shot_end; ipshot++) {
	    itmp=ishot*n2+ipshot;
	    nmask=mask[itmp];	
	       
	    for (ipsample = nmask; ipsample < nsamp_gap; ipsample++) {
	       
dpin[ishot*Ngap*nsamp_gap+(ipshot-ishot+Ngapl)*nsamp_gap+ipsample]+=-(B[itmp][ipsample]-x0s[itmp][ipsample]);	       	       	         
	      }
	    }	   
	   }  	   
/* apply 2D correction on dpin in freq. domain */		   
	   
   	ishot=0;	
					
	while (ishot<nshots) {	

	  sign = -1;
	  scale = 1.0;
			
	  for(ireceiver = 0; ireceiver < Ngap; ireceiver++) {
	    for(j = 0; j < nsamp_gap; j++) rdata_pin[ireceiver*optn+j] = dpin[ishot*Ngap*nsamp_gap+ireceiver*nsamp_gap+j]*scale;
	    for(j = nsamp_gap; j < optn; j++) rdata_pin[ireceiver*optn+j] = 0.0;
	    }

          rcmfft(&rdata_pin[0], (complex *)&tmpdata_pin[0], optn, Ngap, optn, nfreq, sign);

          scale = 1.0/optn;
		
	  for(ireceiver = 0; ireceiver < Ngap; ireceiver++) {
	    for(j = 0; j < nfreq; j++) {
	      cdata_pinc[ireceiver*nfreq+j].r =	corr2D[j]*(tmpdata_pin[ireceiver*nfreq*2+j*2]-tmpdata_pin[ireceiver*nfreq*2+j*2+1])*scale;
	      cdata_pinc[ireceiver*nfreq+j].i =	corr2D[j]*(tmpdata_pin[ireceiver*nfreq*2+j*2]+tmpdata_pin[ireceiver*nfreq*2+j*2+1])*scale;
	     }
	   }		  
		  
          sign = 1;
          crmfft(&cdata_pinc[0], &datout_pinc[0], optn, Ngap, nfreq, optn, sign);
	  			              
          if (nsamp_gap>0 && nsamp_gap<optn) {				
	    for (ireceiver = 0; ireceiver < Ngap; ireceiver++) {		       
	      nmask=mask[MIN(MAX(ishot*n2+(ireceiver+ishot-Ngapl),0),n2*n2-1)];
	      for (j = 0; j < nmask; j++) {
	        dpinc[ishot*Ngap*nsamp_gap+ireceiver*nsamp_gap+j]=0;
               }
	      for (j = nmask; j < nsamp_gap; j++) {
		dpinc[ishot*Ngap*nsamp_gap+ireceiver*nsamp_gap+j]=datout_pinc[ireceiver*optn+j];
	        }	 				
	      }
            }	    
	    		 	  			        	        
	  ishot+=1;							
	} 	   
	    
        t2 = cputime();
	//samess("CPU-time to determine dpin and dpinc = %.3f", t2-t3);
        
  /*  scale dpin by calculating "dpin + x0c*dpin" and match that to the residual. Instead of using x0c we will use dpinc to appply the 2D correction. Use dpminds to save "dpin + x0*dpin". No need to define an extra variable */ 
    
        for (i = 0; i < nshots*n2; i++) {
	  for (j = 0; j < n1; j++) {
	   dpminds[i][j]=0;
	   }
	 }
    
        for (ireceiver = 0; ireceiver < n2; ireceiver++) { 
         recspk=ireceiver; 
          for (ishot = 0; ishot < nshots; ishot++) {
            shot_start=MAX(ishot-Ngapl,0);
	    shot_end=MIN(ishot+Ngapr,nshots);
            shotspk=ishot;
            for (isample = 0; isample < n1; isample++) {
              if (fabs(x0[ishot*nreceivers+ireceiver][isample]) > 0 ) {      
	        sampspk=isample;	  
	        ampspk=x0[shotspk*n2+recspk][sampspk];	 
	        itmppc=shotspk*n1-sampspk;	   	        
                for (ipshot = shot_start; ipshot < shot_end; ipshot++) {	     
	          ix=ipshot*n1;	    
	          istartpc=ipshot*Ngap*nsamp_gap+(shotspk-ipshot+Ngapl)*nsamp_gap-sampspk;	     
	          for (ipsample = sampspk; ipsample < sampspk+nsamp_gap; ipsample++) {	  
	           tmp_recgather[ix+ipsample]+=ampspk*dpinc[istartpc+ipsample];	            
	            }	      	 
	          }
	       }
	    }         
          }
// samess("7.4 tmp_x0[0][0]=%f receiver=%d  n2=%d",tmp_x0[0][0],ireceiver,n2);      
  	  for (i = 0; i < nshots; i++) {
	    ix=i*n1;
	    it=i*n2;
	    for (j = 0; j < n1; j++) {
	     dpminds[it+ireceiver][j]=tmp_recgather[ix+j];
	     tmp_recgather[ix+j]=0;	   
	     }
	   }              
         } 
//  samess("7.5 tmp_x0[0][0]=%f",tmp_x0[0][0]);   
 /*     add dpin to Iplusxodpin,    */
       for (ishot = 0; ishot < nshots; ishot++) { 	 
	 for (ireceiver = MAX(ishot-Ngapl,0); ireceiver < MIN(ishot+Ngapr,nreceivers); ireceiver++) {
          for (isample = 0; isample < nsamp_gap; isample++) {	    
	     dpminds[ishot*nreceivers+ireceiver][isample]+=dpin[ishot*Ngap*nsamp_gap+(ireceiver-ishot+Ngapl)*nsamp_gap+isample];	     
	     	     	     	 
	  } 	   
         } 
       }
       
       tmp1=0;
       tmp2=0;
       for (ishot = 0; ishot < nshots; ishot++) { 
        for (ireceiver = 0; ireceiver < n2; ireceiver++) {
         for (j = 0; j < n1; j++) {
           tmp1+=mask_beta[ishot][ireceiver]*(dpminds[ishot*n2+ireceiver][j])*(B[ishot*n2+ireceiver][j]-x0s[ishot*n2+ireceiver][j]);
	   tmp2+=mask_beta[ishot][ireceiver]*(dpminds[ishot*n2+ireceiver][j])*(dpminds[ishot*n2+ireceiver][j]);	 
	  }  
         } 
       }   

       beta=-tmp1/tmp2;
       /*if (verbose) {
       samess("beta = %.5f", beta);   
       } */
            
  	for (i = 0; i < nshots*n2; i++) {
	  for (j = 0; j < n1; j++) {
	   B[i][j]+=beta*dpminds[i][j];
	  }
	 }   
	 	 
         for (i = 0; i < nshots*Ngap*nsamp_gap; i++) {	    
	     pin[i]+=beta*dpin[i];	     
	     pinc[i]+=beta*dpinc[i];	     	     	 
	  } 
	  
	 t3 = cputime();
         //samess("CPU-time to determine dpin + x0*dpin, scale and update p,pc and residual = %.3f", t3-t2);      
 	 
     } 	 /* end if mode ==1 */ 
/*  Bring residu and dpinc to frequency domain to determine dx0. Use Preal and Pim to store residu.     */ 
	
	ishot=0;	
					
	while (ishot<nshots) {	
          sign = -1;
          scale = 1.0;
	  for(ireceiver = 0; ireceiver < n2; ireceiver++) {
	    nmask=mask[ishot*n2 + ireceiver];
	    for(j = 0; j < nmask; j++) rdata[ireceiver*optn+j] = 0.0;
	    for(j = nmask; j < n1; j++)
	      rdata[ireceiver*optn+j] = (B[ishot*n2+ireceiver][j]-x0s[ishot*n2+ireceiver][j])*scale;
	    for(j = n1; j < optn; j++) rdata[ireceiver*optn+j] = 0.0;
           }
	  rcmfft(&rdata[0], (complex *)&tmpdata[0], optn, n2, optn, nfreq, sign);
								
	  for(i = 0; i < nreceivers; i++) {
	    for(j = 0; j < nfreq; j++) {
	      Preal[j*nreceivers*nshots+i*nshots+ishot] = tmpdata[i*nfreq*2+j*2];
              Pim[j*nreceivers*nshots+i*nshots+ishot] = tmpdata[i*nfreq*2+j*2+1];
	      }
	    }
/*  get dpinc in frequency domain, (apply 2D correction is allready done!) and add to near offsets of Pc   */
	  if (mode==1) {	 	
	    for(ireceiver = 0; ireceiver < Ngap; ireceiver++) {
	      for(j = 0; j < nsamp_gap; j++) 
	        rdata_pin[ireceiver*optn+j] = beta*dpinc[ishot*Ngap*nsamp_gap+ireceiver*nsamp_gap+j]*scale;
	      for(j = nsamp_gap; j < optn; j++) rdata_pin[ireceiver*optn+j] = 0.0;
	     }
	     
            rcmfft(&rdata_pin[0], (complex *)&tmpdata_pin[0], optn, Ngap, optn, nfreq, sign);
										
	    for (ireceiver = MAX(ishot-Ngapl,0); ireceiver < MIN(ishot+Ngapr,nreceivers); ireceiver++) {
	      for(j = 0; j < nfreq; j++) {		 
	        Pcreal[j*nreceivers*nshots+ireceiver*nshots+ishot] += tmpdata_pin[(ireceiver-ishot+Ngapl)*nfreq*2+j*2];
		Pcim[j*nreceivers*nshots+ireceiver*nshots+ishot] += tmpdata_pin[(ireceiver-ishot+Ngapl)*nfreq*2+j*2+1];
		}
	      }
	   }	
	 ishot+=1;							
	}	
	
          t2 = cputime();
	 //samess("CPU-time to get residu and pin in frequency domain and ifft pinc = %.3f", t2-t3);
	  		
/*  correlation  */
 		
	tmp1=0;
	tmp2=0;	
		
        for(ifreq = 0; ifreq < nfreq; ifreq++) {
	    index_freq_start=ifreq*nreceivers*nshots;
	  for (ishot = 0; ishot < nshots; ishot++) {
	    index_out_shot=2*ishot*nreceivers*nfreq+2*ifreq;
	    index_shot_start=index_freq_start+ishot*nshots;
	    for (ireceiver = 0; ireceiver < nreceivers; ireceiver++) {
	      index_receiver_start=index_freq_start+ireceiver*nshots;
	      index_out=index_out_shot+2*ireceiver*nfreq;
	     for (ix = 0; ix < nshots; ix++) {	                
	      tmp1+=Preal[index_receiver_start+ix]*Pcreal[index_shot_start+ix]+Pim[index_receiver_start+ix]*Pcim[index_shot_start+ix];	      	    
	      tmp2+=Pim[index_receiver_start+ix]*Pcreal[index_shot_start+ix]-Preal[index_receiver_start+ix]*Pcim[index_shot_start+ix];							              
	      }
	      PPcH[index_out]=tmp1;
	      PPcH[index_out+1]=tmp2;
	      tmp1=0;
	      tmp2=0;  
	    }
	  }
        }   
	

 	t3 = cputime();
	//samess("CPU-time to apply correlation in frequency domain = %.3f", t3-t2);
		  	
      ishot=0;
			
	while (ishot<nshots) {		
		
	  scale = 1.0/optn;
	  for(i = 0; i < n2; i++) {
	    for(j = 0; j < nfreq; j++) {
	      cdata[i*nfreq+j].r = PPcH[2*ishot*n2*nfreq+2*i*nfreq+2*j]*scale;
	      cdata[i*nfreq+j].i = PPcH[2*ishot*n2*nfreq+2*i*nfreq+2*j+1]*scale;
	      }
	    }		  
		  
          sign = 1;

	  crmfft(&cdata[0], &datout[0], optn, n2, nfreq, optn, sign);
				             
           if (n1_orig>0 && n1_orig<optn) {
	    for (i = 0; i < n2; i++) {
	     for (j = 0; j < n1_orig; j++) {
	       dx0[ishot*n2*n1+i*n1+j]=-datout[i*optn+j];			    
	       }	 
              }  
	     }
	   ishot+=1;	 
	}	 	 

	t2 = cputime();
	//samess("CPU-time to ifft dx0 = %.3f", t2-t3); 
		  
/* update pc and p by replacing the near offsets */	  	 
  if (mode==1) {	
       for (ishot = 0; ishot < nshots; ishot++) { 	 
	 for (ireceiver = MAX(ishot-Ngapl,0); ireceiver < MIN(ishot+Ngapr,nreceivers); ireceiver++) {
          for (isample = 0; isample < nsamp_gap; isample++) {	    	     pc[ishot*nreceivers*n1+ireceiver*n1+isample]=pinc[ishot*Ngap*nsamp_gap+(ireceiver-ishot+Ngapl)*nsamp_gap+isample];
p[ishot*nreceivers*n1+ireceiver*n1+isample]=pin[ishot*Ngap*nsamp_gap+(ireceiver-ishot+Ngapl)*nsamp_gap+isample];	  
	  
	    } 	   
         } 
       }
     }  	 

   /*add residu*sH to dxo  */	
		
	for (ix=0;ix<nshots*n2;ix++) {		
	  //for (itout=lfilt-1-nrot;itout<n1-nrot;itout++) {
	   nmask=mask[ix];
	  for (itout=MAX(nrot,nmask);itout<n1-lfilt+nrot+1;itout++) {
	    itin=itout-nrot;		   
	    for (it=0;it<lfilt;it++) {
	      if (nwindows > 1) {	        	      
	       for (iwindow = 0; iwindow < nwindows; iwindow++) {
	         if (itin >= iwindow*nsamp_window && itin < (iwindow+1)*nsamp_window) iwav=iwindow;  
	        }	      
                 if (itin >= (nwindows-1)*nsamp_window) iwav=nwindows-1;
	      } 		    
	     dx0[ix*n1+itout] += (B[ix][itin]-x0s[ix][itin])*s[iwav][it];
	     itin++;
	      }
	    }
	  }       	 
  /*   calculate objective function    */	
       tmp1=0;
       tmp2=0;
       for (ishot = 0; ishot < nshots; ishot++) { 
         for (ireceiver = 0; ireceiver < n2; ireceiver++) {
          for (j = 0; j < n1; j++) {
	    tmp1+=(B[ishot*n2+ireceiver][j]-x0s[ishot*n2+ireceiver][j])*(B[ishot*n2+ireceiver][j]-x0s[ishot*n2+ireceiver][j]);
	  }  
         } 
       }	 
	 
       if (verbose) samess("objective function value after iteration %d = %.3f",iter+1, tmp1);	 		 
       t3 = cputime();	
       
/*  write wavelet(s) after each iteration  */       
//samess("10.tmp_x0[0][0]=%f",tmp_x0[0][0]);         
	if (file_wav != NULL) {		   
	  if (QC==0 && file_QC_dx0==NULL) {
	    dom1 = SA_AXIS_TIME;
	    type = SA_TYPE_REAL;
	    trid = (short) axis2trid(dom1, dom2, type);			
	  }	  
	  hdrs_in[0].duse = iter+1;
	for (iwav = 0; iwav < nwindows; iwav++) {
	  hdrs_in[0].year = iwav+1;
	 for  (it = 0; it < lfilt; it++) s_tmp[it]=s[iwav][it];
	  ret = write_data(file_wav,s_tmp,lfilt,1,f1,f2,d1,d2,type,hdrs_in);
	  if (ret < 0 ) saerr("error on writing output file_wav.");	
	  }       
        }
       
       	 
       //if (verbose) samess("Total CPU-time after iteration %d = %.3f",iter+1, t3-tstart);
       t2 = cputime();
	
       first = 0;
       iter+=1;
       if (verbose) samess("iteration %d finished", iter);
		
     }
	
/*  write shot records to output file  */	
   
       
      ishot=0; 	
	while (ishot<nshots) {	

	  for (i = 0; i < n2; i++) {
	    for (j = 0; j < n1; j++) {
	      datout[i*n1_orig+j]=x0s[ishot*n2+i][j];
	      }	 				
	    }
	    
	  dom1 = SA_AXIS_TIME;
	  type = SA_TYPE_REAL;
	  trid = (short) axis2trid(dom1, dom2, type);
	  /*n1 = optn;*/
	  /*d1 = d1/2;*/
	  f1 = 0.0;	
	  f2=(float) hdrs_tot[ishot*n2].offset;	
				 
  /* update headers for output file */
	  for(i = 0; i < n2; i++){
	    hdrs_in[i]=hdrs_tot[ishot*n2+i];
	  }			

/* Write extra data to output */

	  ret = set_keys(keys, seqnr, nkeys);
	  if (ret < 0 ) sawarn("error on writing keys.");
	  ret = set_axis(dom1, dom2);
	  if (ret < 0 ) saerr("error on writing axis.");

	  ret = write_data(file_prim,datout,n1,n2,f1,f2,d1,d2,type,hdrs_in);
	  if (ret < 0 ) saerr("error on writing output primaries");
		
	  /* write data residuals */
	  if (file_res != NULL) {
	    for (i = 0; i < n2; i++) {
	      for (j = 0; j < n1; j++) {
               datout[i*n1_orig+j]=B[ishot*n2+i][j]-x0s[ishot*n2+i][j];
	        }	 				
	      }
	    ret = write_data(file_res,datout,n1,n2,f1,f2,d1,d2,type,hdrs_in);
	    if (ret < 0 ) saerr("error on writing output residual");
	   }
		
	  /* write input data with reconstructed near offsets */
	  if (file_p_recon != NULL) {
	   for (i = 0; i < n2; i++) {
	    for (j = 0; j < n1; j++) {	
	      datout[i*n1_orig+j]=p[ishot*n2*n1+i*n1+j];
	      }	 				
	     } 
	   ret = write_data(file_p_recon,datout,n1,n2,f1,f2,d1,d2,type,hdrs_in);
	     if (ret < 0 ) saerr("error on writing output reconstructed data");
	   }

	  /* write estimated multiples */
	  if (file_mult != NULL) {
	   for (i = 0; i < n2; i++) {
	     for (j = 0; j < n1; j++) {
	       datout[i*n1_orig+j]=-B[ishot*n2+i][j]+p[ishot*n2*n1+i*n1+j];
               }
	     }       
	   ret = write_data(file_mult,datout,n1,n2,f1,f2,d1,d2,type,hdrs_in);
	   if (ret < 0 ) saerr("error on writing output multiples");	
	   }

	  /* write estimated impulsed response spikes x0 */
	  if (file_x0 != NULL) {
	   for (i = 0; i < n2; i++) {
	     for (j = 0; j < n1; j++) {
	       datout[i*n1_orig+j]=x0[ishot*n2+i][j];
               }
	     }       
	   ret = write_data(file_x0,datout,n1,n2,f1,f2,d1,d2,type,hdrs_in);
	   if (ret < 0 ) saerr("error on writing output x0");	
	   }
	   
	 first = 0;
         ishot+=1;	
						
	}
	
	/* if (file_wav != NULL) {
	ret = write_data(file_wav,s,lfilt,1,f1,f2,d1,d2,type,hdrs_in);
		if (ret < 0 ) saerr("error on writing output file_wav.");	
	}	*/		
        ret = close_file(file_prim);
	if (ret < 0) sawarn("err %d on closing output file",ret);
	if (file_res != NULL) {
	ret = close_file(file_res);
	if (ret < 0) sawarn("err %d on closing output file",ret);
	}
	if (file_mult != NULL) {
	ret = close_file(file_mult);
	if (ret < 0) sawarn("err %d on closing output file",ret);
	}
	if (file_p_recon != NULL) {
	ret = close_file(file_p_recon);
	if (ret < 0) sawarn("err %d on closing output file",ret);
	}	
	if (file_QC != NULL && QC > 0) {
	ret = close_file(file_QC);
	if (ret < 0) sawarn("err %d on closing output file",ret);
	}	
	if (file_wav != NULL) {
	ret = close_file(file_wav);
	if (ret < 0) sawarn("err %d on closing output file",ret);
	}
	if (file_QC_dx0 != NULL) {
	ret = close_file(file_QC_dx0);
	if (ret < 0) sawarn("err %d on closing output file",ret);
	}	
			
	free(hdrs_in);
	free(hdrs_tot);

	free(pc);
	free(tmpdata);
	free(datout);
	free(mask);
	free(corr2D);
	free(offsettaper);
	free(xtaper);
	free(ttaper);
	free(wortelT);
	free(dx0_shot);
	free(dx0_receiver);
	free(dx0_sample);
	free(dx0_amp);
	free(dx0);
	free(dpminds);
	free(x0);
	free(B);
	free(mask_s);
	free(mask_alfa);
	free(s);
	free(QC_shots);
	free(rdata);
	free(cdata);
	free(x0s);	
	free(s_rev);
	free(tmp_recgather);

	if (mode==1) {
	free(tmpdata_pin);
	free(datout_pinc);
	free(no_orig);	
	free(pin);
	free(pinc);
	free(dpin);
	free(dpinc);
	free(rdata_pin);
	free(cdata_pinc); 	
	free(mask_beta);
	}
	 
	     
	for (ikey=0; ikey<MAX_KEYS; ikey++) free(keys[ikey]);
	t2 = cputime();
	if (verbose) samess("Total CPU-time = %.3f", t2-tstart);

	exit(0);
}
