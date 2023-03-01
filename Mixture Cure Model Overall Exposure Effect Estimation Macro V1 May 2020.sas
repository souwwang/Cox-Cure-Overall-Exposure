*************************************************************************************************************************************************************************** ;
*	ESTIMATING OVERALL EXPOSURE EFFECTS FOR COX PROPORTIONAL HAZARDS CURE MODELS WITH INTERVAL-CENSORED SURVIVAL DATA                                                       ;
*   Version 1.1, Uploaded 02/28/2023 												                                                                                        ;
*																					               										                                    ;
*	Code written by Wei Wang (2023)													               										                                    ;
* 																					                										                                ;
*	Reference:																		               										                                    ;
* 	Wang W, Cong N, Ye A, Zhang H, Zhang B. Exposure assessment for Cox proportional hazards cure models with interval-censored survival data.                              ;
* 	Biom J. 2022 Jan;64(1):91-104. 														                                                                                    ;
* 																					                                                                                        ;
*************************************************************************************************************************************************************************** ;
* 																					                										                                ;
*	PROGRAM DESCRIPTION:															                                                                                        ;
*   ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- ;
* 	This program applied a post-estimation approach ('average-predicted value' method) by using model-predicted estimates difference to assess the overall exposure 	    ;
* 	effect on the restricted mean survival time scale for Cox proportional hazards cure model to analyze interval-censored survival data in the presence of a cured         ;
* 	fraction. For the model parameter estimation, simple parametric models as fractional polynomials or restricted cubic splines are utilized to approximate the baseline   ;
* 	logarithm cumulative hazard function, or alternatively, the full likelihood is specified through a piecewise linear approximation for the cumulative baseline hazard    ;
* 	function. This SAS macro can be applied for one binary exposure, one interval-censored outcome and covariates to predict the susceptible probability and the            ;
* 	non-susceptible survival time.																				                                                            ;
* 																					                                                                                        ;
*	DATASET REQUIREMENTS:															                                                                                        ;
*   ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- ;
*	The macro requires a SAS data set with the following variables:                                                                                                         ;
*   (1) Survival outcome variable 1: left endpoint (SURVLtime)                                                                                                              ;
*   (2) Survival outcome variable 2: right endpoint (SURVRtime)                                                                                                             ;
*   (3) Survival event indicator (SURVevent)                                                                                                                                ;
*   (4) List of covariates for the susceptible probability (CURECov)                                                                                                        ;
*   (5) List of covariates for the non-suscpetible survival outcome (SURVCov)                                                                                               ;
*   (6) Binary exposure (EXP).                                                                                                                                              ;
*                                                                                                                                                                           ;
*   Of note, the covariates for the susceptible probability and the non-susceptible survival outcome can be same, or it can be different.                                   ;
* 																					                										                                ;
*   Format of data set exmaple:                                                                                                                                             ;
* 																					                										                                ;
*	ID	SURVLTIME  SURVRTIME    SURVevent    Cov1     Cov2     Cov3        EXP                                                                                			    ;
*	1	3.2	       5.5          1            13.2	    1		 21		     0			                                                                                    ;
*	2	10.5	   .            0            16.9		0		 23		     1				                                                                                ;
*	3	7.1	       10.6         1            12.0		0		 23		 	 0					                                                                            ;
*	4	8.2	       .            0            23.5		1		 12		 	 1					                                                                            ;
*	5	4.5	       .            0            21.0		1		 13		 	 1				                                                                                ;
*	6	12.0	   15.6         1            12.7		1		 45		 	 1						                                                                        ;
*																					                                                                                        ;
*	MODEL SPECIFICATION																                                                                                        ;
*   ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- ;
*	X: binary exposure                                                                                                                                                      ;
*   W1: covariates for the susceptible probability                                                         			                                                        ;
*   W2: baseline covariate for the non-susceptible survival outcome                                                			                                                ;
*   Y: latent binary indicator with Y = 1 indicating that the subject will eventuallly experience the failure event (uncured), otherwise, if Y = 0 indicating that the      ;
*   individual will never experience such event (cured)                                                                    			                                        ;
*   T: time-to-event outcome                                                                                       			                                                ;
*																					                                                                                        ;
*	Uncured probability, logit(y) = alpha0 + alpha1 * x + alpha2 * w1    	   				                                                                                ;
*																					                                                                                        ;
*   Proportional hazards model for time-to-event outcome defined through the hazard function, h(t) = h0(t)exp(beta1 * x + beta2 * w2)                                       ;
*																					                                                                                        ;
*	MACRO VARIABLES:																                                                                                        ;
*   ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- ;
*	DATA: 		    SAS data set to contain all analysis variables                                                                              	                        ;
*																					                                                                                        ;
*	subjid:		    subject id vairable 	                                                                                                                                ;
*																					                                                                                        ;
*	CURECov:	    List of covariates for the survival susceptibility                                                                                                      ;
*																					                                                                                        ;
*	SURVCov:	    List of covariates for the interval-censored survival outcome                                                                                           ;
*																					                                                                                        ;
*	EXP: 		    binary exposure                                      	                                                                                                ;
*																					                                                                                        ;
*	Nboot:          number of bootstrap samples to calculate the 95% CI for the FP/RCS parametric approximation approach , default Nboot = 200 (suggested)	                ;
*																					                                                                                        ;
*	SURVLtime:		interval censored time-to-event variable left endpoint (should not missing)	                                                                            ;
*																					                                                                                        ;
*	SURVRtime:		interval censored time-to-event variable right endpoint (can be missing depending on whether the survival time is interval-censored or right-censored)  ;
*																					                                                                                        ;
*	SURVevent:		event indicator, event = 1, censor = 0. Of note, if event = 1 (interval censored), SURVRtime should be non-missing and SURVRtime is greater             ;
*					than SURVLtime. If event = 0 (right censored), SURVRtime should be missing.																                ;
*																					                                                                                        ;
*	T_ACCESS:		t*, time point at which the overall exposure effects are defined, default (999999) is the last observed event time, the macro users may replace 999999  ;
*					with other arbitrarily designated t* (t* > 0, strongly recommended to be less than the last observed event time, see discussion section of published    ;
*					reference by Wang et al. SMMR)								                                                                                            ;
*																					                                                                                        ;
*	REFERENCE:		can be 0, 1, 2 or 3. The estimated overall exposure effects are determined by the baseline covariate (predictors for the susceptibility probability and ;
*					survival outcome, heterogeneous overall exposure effect estimation for different baseline covariates) and the estimated 			    				;
*					and presented average overall exposure effects are calculated as the average of overall exposure effects over the empirical covariates distribution from;
*					the designated reference population (see reference for more details). The reference population can be exposed group, unexposed group, total             ;
*					population or any designated population with arbitrarily determined baseline covariate distribution, default vaule for reference is 1.                  ;
*					0: the unexposed group as the reference population, and the average overall exposure effects are calculated over the empirical covariates distribution  ;
*					   of the unexposed group															                                                                    ;
*					1: the exposed group as the reference population, and the average overall exposure effects are calculated over the empirical covariates distribution    ;
*					   of the exposed group															                                                                        ;
*					2: the total population as the reference group, and the average overall exposure effects aare calculated over the empirical covariates distribution     ;
*					   of the total population   															                                                                ;
*					3: users provided population with pre-specified baseline covariate distribution as the reference group. If reference is 3, the Macro user should        ;
*					   provide a temporary dataset named as REF to show the baseline covariate distribution in which the users want to estimate the overall exposure effects;
*					   The dataset REF should follow the following format:															                                        ;
*																					                                                                                        ;
*	                          Cov1         Cov2         Cov3         NSUB                                                                                      			    ;
*                             0		       0		    0			 10                                                                                                     ;
*                             0		       0		    1			 10                                                                                                     ;
*                             0		       1		    0			 20                                                                                                     ;
*                             0		       1		    1			 20                                                                                                     ;
*                             1		       0		    0			 10                                                                                                     ;
*                             1		       0		    1			 10                                                                                                     ;
*                             1		       1		    0			 10                                                                                                     ;
*                             1		       1		    1			 10                                                                                                     ;
*																					                                                                                        ;
*					   If the dataset REF is defined as above, the reference population is defined as following, in this reference population,                              ;
*					   number of subjects with Cov1 = 0, Cov2 = 0 and Cov3 = 0 is 10, number of subjects with Cov1 = 0, Cov2 = 0 and                                        ;
*					   Cov3 = 1 is 10,	number of subjects with Cov1 = 0, Cov2 = 1 and Cov3 = 0 is 20 etc.                                                                  ;
*																					                                                                                        ;
*	SEEDBOOT:		the random-number seed to generate bootstrap samples 	                                                                                                ;
*																					                                                                                        ;
*	PRINTALL:		Outpout data set print out or not	                                                                                                                    ;
*					F: output data set is not printed out								                                                                                    ;
*					T: output data set is printed out, default												                                                                ;
*																					                                                                                        ;
*	OUT:	        Output dataset                                                                                                                                          ;
*																					                                                                                        ;
*																					                                                                                        ;
*	EXAMPLE CODE:																	                                                                                        ;
*																					                                                                                        ;
*   %include 'cureexp.sas'														                                                                                            ;
*																					                                                                                        ;
*   %CUREEXP   (data             = aa11,                                                                                                                                    ;
*               subjid           = ObsLHS,                                                                                                                                  ;
*               CUREcov          = sexF Duration F10Cigs,                                                                                                                   ;  
*               SURVcov          = sexF Duration F10Cigs,                                                                                                                   ;
*               EXP              = SI_UC,                                                                                                                                   ;
*               nboot            = 200,                                                                                                                                     ;
*               SURVLtime        = timept1,                                                                                                                                 ;
*               SURVRtime        = timept2,                                                                                                                                 ;
*               SURVevent        = relapse,                                                                                                                                 ;
*			    T_ACCESS         = 999999,                                                                                                                                  ;
*			    REFERENCE        = 1,                                                                                                                                       ;
*               seedboot         = 1617893,                                                                                                                                 ;
*               Printall         = T,                                                                                                                                       ;
*               Out             = out01                                                                                                                                     ;
*               );                                                                                                                                                          ;
* 																					                                                                                        ;
*************************************************************************************************************************************************************************** ;
* 																					                                                                                        ;
*************************************************************************************************************************************************************************** ;
***********************************************************************************Final Macro******************************************************************************;
*************************************************************************************************************************************************************************** ;

dm  'log;clear;out;clear;';

OPTIONS ls=130  ps=57 NOCENTER DATE PAGENO=1; 

** Define SubMacro;

** RCS is sub Macro to fit the survival outcome with the restricted cubic splines with 1 or 2 or 3 interior knots to approximate the baseline log cumulative hazard function. 
Parameter estimates and fit criteria (AIC and BIC) are provided in the output data set;

%macro rcs();

/*%put &linpredCURE &linpredSURV;*/

  proc nlmixed data=aa53 cov;

    parms / bydata data=Init1;
    linpredSURV  = &linpredSURV; 
    linpredcrs = gam10 + gam11 * lnt + gam12 * nu11;
    linpred = linpredsurv + linpredcrs;
    linpredcrsR = gam10 + gam11 * lnRt + gam12 * nuR11;
    linpredR = linpredsurv + linpredcrsR;
	cureprop = exp(min(&linpredCURE, 600))/(1+exp(min(&linpredCURE, 600)));

	STL = exp(-exp(min(linpred, 600)));
	STR = exp(-exp(min(linpredR, 600)));
	logevent = log(cureprop) + log(STL - STR);
	logcensor = log((1 - cureprop) + cureprop * STL);

    if &SURVevent = 1 then loglike = logevent;
	if &SURVevent = 0 then loglike = logcensor;

	model &SURVevent ~ general(loglike);
 
    ods output ParameterEstimates=temp15;
    ods output FitStatistics = temp16; 
    ods output ConvergenceStatus = temp17 (keep = replicate status rename = (status = value));
	by replicate;
  run;

  %if %sysfunc(exist(temp16)) %then %do;

    data temp18;
      set temp16 temp17 (in = a);
      if a then Descr = 'Con';
    run;

    proc sort data = temp18;
      by replicate;
    run;

    proc transpose data = temp18 out = temp19 (rename = (AIC__smaller_is_better_ = AIC AICC__smaller_is_better_ = AICC BIC__smaller_is_better_ = BIC));
      id descr;
      var value;
      by replicate;
    run;

    data gg45;
      set temp19;
      drop _name_;
      if con = 0.0 then converge = 'Yes';
      else converge = 'No';
      drop con;
      method = 'RCS';
      nknot = 1;
      ID = 45;
    run;

    data ff45;
      set temp15;
      method = 'RCS';
      nknot = 1;
      ID = 45;
    run;

    proc datasets lib = work;
      delete temp15 temp16 temp17 temp18 temp19;
    run;

  %end;

  %else %do;

    data gg45;
      format converge $3.;
      converge = 'No';
      method = 'RCS';
      nknot = 1;
      ID = 45;
    run;

    data ff45;
      format converge $3.;
      converge = 'No';
      method = 'RCS';
      nknot = 1;
      ID = 45;
    run;

  %end;

  proc nlmixed data=aa53 cov;

    parms / bydata data=Init2;
    linpredSURV  = &linpredSURV; 
    linpredcrs = gam20 + gam21 * lnt + gam22 * nu21 + gam23 * nu22;
    linpred = linpredsurv + linpredcrs;
    linpredcrsR = gam20 + gam21 * lnRt + gam22 * nuR21 + gam23 * nuR22;
    linpredR = linpredsurv + linpredcrsR;
	cureprop = exp(min(&linpredCURE, 600))/(1+exp(min(&linpredCURE, 600)));

	STL = exp(-exp(min(linpred, 600)));
	STR = exp(-exp(min(linpredR, 600)));
	logevent = log(cureprop) + log(STL - STR);
	logcensor = log((1 - cureprop) + cureprop * STL);

    if &SURVevent = 1 then loglike = logevent;
	if &SURVevent = 0 then loglike = logcensor;

	model &SURVevent ~ general(loglike);

    ods output ParameterEstimates=temp15;
    ods output FitStatistics = temp16; 
    ods output ConvergenceStatus = temp17 (keep = replicate status rename = (status = value));
    by replicate; 
  run;

  %if %sysfunc(exist(temp16)) %then %do;

    data temp18;
      set temp16 temp17 (in = a);
      if a then Descr = 'Con';
    run;

    proc sort data = temp18;
      by replicate;
    run;

    proc transpose data = temp18 out = temp19 (rename = (AIC__smaller_is_better_ = AIC AICC__smaller_is_better_ = AICC BIC__smaller_is_better_ = BIC));
      id descr;
      var value;
      by replicate;
    run;

    data gg46;
      set temp19;
      drop _name_;
      if con = 0.0 then converge = 'Yes';
      else converge = 'No';
      drop con;
      method = 'RCS';
      nknot = 2;
      ID = 46;
    run;

    data ff46;
      set temp15;
      method = 'RCS';
      nknot = 2;
      ID = 46;
    run;

    proc datasets lib = work;
      delete temp15 temp16 temp17 temp18 temp19;
    run;

  %end;

  %else %do;

    data gg46;
      format converge $3.;
      converge = 'No';
      method = 'RCS';
      nknot = 2;
      ID = 46;
    run;

    data ff46;
      format converge $3.;
      converge = 'No';
      method = 'RCS';
      nknot = 2;
      ID = 46;
    run;

  %end;

  proc nlmixed data=aa53 cov;

    parms / bydata data=Init3;
    linpredSURV  = &linpredSURV; 
    linpredcrs = gam30 + gam31 * lnt + gam32 * nu31 + gam33 * nu32 + gam34 * nu33;
    linpred = linpredsurv + linpredcrs;
    linpredcrsR = gam30 + gam31 * lnRt + gam32 * nuR31 + gam33 * nuR32 + gam34 * nuR33;
    linpredR = linpredsurv + linpredcrsR;
	cureprop = exp(min(&linpredCURE, 600))/(1+exp(min(&linpredCURE, 600)));

	STL = exp(-exp(min(linpred, 600)));
	STR = exp(-exp(min(linpredR, 600)));
	logevent = log(cureprop) + log(STL - STR);
	logcensor = log((1 - cureprop) + cureprop * STL);

    if &SURVevent = 1 then loglike = logevent;
	if &SURVevent = 0 then loglike = logcensor;

	model &SURVevent ~ general(loglike);

    ods output ParameterEstimates=temp15;
    ods output FitStatistics = temp16; 
    ods output ConvergenceStatus = temp17 (keep = replicate status rename = (status = value));
	by replicate;
  run;

  %if %sysfunc(exist(temp16)) %then %do;

    data temp18;
      set temp16 temp17 (in = a);
      if a then Descr = 'Con';
    run;

    proc sort data = temp18;
      by replicate;
    run;

    proc transpose data = temp18 out = temp19 (rename = (AIC__smaller_is_better_ = AIC AICC__smaller_is_better_ = AICC BIC__smaller_is_better_ = BIC));
      id descr;
      var value;
      by replicate;
    run;

    data gg47;
      set temp19;
      drop _name_;
      if con = 0.0 then converge = 'Yes';
      else converge = 'No';
      drop con;
      method = 'RCS';
      nknot = 3;
      ID = 47;
    run;

    data ff47;
      set temp15;
      method = 'RCS';
      nknot = 3;
      ID = 47;
    run;

    proc datasets lib = work;
      delete temp15 temp16 temp17 temp18 temp19;
    run;

  %end;

  %else %do;

    data gg47;
      format converge $3.;
      converge = 'No';
      method = 'RCS';
      nknot = 3;
      ID = 47;
    run;

    data ff47;
      format converge $3.;
      converge = 'No';
      method = 'RCS';
      nknot = 3;
      ID = 47;
    run;

  %end;

%mend;

** SINGLE is sub Macro to fit the survival outcome with the fractional polynomials with degree 1 to approximate the baseline log cumulative hazard function. 
Parameter estimates and fit criteria (AIC and BIC) are provided in the output data set;

%macro single(in1, in2, out1, out3, out2, n, id);

/*  dm  'log;clear;out;clear;';*/

  data temp21;
    set &in2;
  run;

  proc glm data = temp21;
    model logbasecumhaz = t&n;
    by replicate;
    ods output ParameterEstimates  = temp211;
  run;

  data temp212;
    set temp211;
    if parameter = 'Intercept' then parameter1 = 'lamd0';
    if parameter = "t&n" then parameter1 = 'lamd1';
	drop parameter;
	rename parameter1 = parameter;
    keep replicate parameter1 estimate;
  run;

  data temp213;
    set SURVests2 (keep = replicate model parameter estimate effect) temp212 (in = a);
    if a then model = 'FP';
  run;
  
  data temp214;
    set  beta3 temp213;
  run;

  proc sort data = temp214;
    by replicate;
  run;

  proc nlmixed data=&in1 cov;
    parms /bydata data=temp214;

    linpredSURV  = &linpredSURV; 
    linpredcrs = lamd0 + lamd1 * t&n;
    linpred = linpredsurv + linpredcrs;
    linpredcrsR = lamd0 + lamd1 * tR&n;
    linpredR = linpredsurv + linpredcrsR;
	cureprop = exp(min(&linpredCURE, 600))/(1+exp(min(&linpredCURE, 600)));

	STL = exp(-exp(min(linpred, 600)));
	STR = exp(-exp(min(linpredR, 600)));
	logevent = log(cureprop) + log(STL - STR);
	logcensor = log((1 - cureprop) + cureprop * STL);

    if &SURVevent = 1 then loglike = logevent;
	if &SURVevent = 0 then loglike = logcensor;

	model &SURVevent ~ general(loglike);
 
    ods output ParameterEstimates=temp25;
    ods output FitStatistics = temp26; 
    ods output ConvergenceStatus = temp27 (keep = replicate status rename = (status = value));
	by replicate;
  run;

  %if %sysfunc(exist(temp26)) %then %do;

    data temp28;
      set temp26 temp27 (in = a);
      if a then Descr = 'Con';
    run;

    proc sort data = temp28;
      by replicate;
    run;

    proc transpose data = temp28 out = temp29 (rename = (AIC__smaller_is_better_ = AIC AICC__smaller_is_better_ = AICC BIC__smaller_is_better_ = BIC));
      id descr;
      var value;
      by replicate;
    run;

    data &out2;
      set temp29;
      drop _name_;
      if con = 0.0 then converge = 'Yes';
      else converge = 'No';
      drop con;
      method = 'FP';
      P1 = &n;
      ID = &ID;
    run;

    data &out1;
      set temp25;
      method = 'FP';
      P1 = &n;
      ID = &ID;
    run;

    proc datasets lib = work;
      delete temp25 temp26 temp27 temp28 temp29;
    run;

  %end;

  %else %do;

    data &out2;
      format converge $3.;
      converge = 'No';
      method = 'FP';
      P1 = &n;
      ID = &ID;
    run;

    data &out1;
      method = 'FP';
      P1 = &n;
      ID = &ID;
    run;

  %end;

  proc datasets lib = work;
    delete temp21 temp211 temp212 temp213 temp214;
  run;

%mend;

** DOUBLE is sub Macro to fit the survival outcome with the fractional polynomials with degree 2 to approximate the baseline log cumulative hazard function. 
Parameter estimates and fit criteria (AIC and BIC) are provided in the output data set;

%macro double(in1, in2, out1, out3, out2, m, n, id);

/*  dm  'log;clear;out;clear;';*/

  data temp31;
    set &in2;
  run;

  proc glm data = temp31;
    model logbasecumhaz = t&m t&n;
    by replicate;
    ods output ParameterEstimates  = temp311;
  run;

  data temp312;
    set temp311;
    if parameter = 'Intercept' then parameter1 = 'lamd0';
    if parameter = "t&m" then parameter1 = 'lamd1';
    if parameter = "t&n" then parameter1 = 'lamd2';
	drop parameter;
	rename parameter1 = parameter;
    keep replicate parameter1 estimate;
  run;

  data temp313;
    set SURVests2 (keep = replicate model parameter estimate effect) temp312 (in = a);
    if a then model = 'FP';
  run;

  proc sort data = temp313;
    by replicate;
  run;

  data temp314;
set  beta3 temp313;
run;

  proc sort data = temp314;
    by replicate;
  run;

  proc nlmixed data=&in1 cov;

    parms /bydata data=temp314;

	linpredSURV  = &linpredSURV; 
    linpredcrs = lamd0 + lamd1 * t&m + lamd2 * t&n;
    linpred = linpredsurv + linpredcrs;
    linpredcrsR = lamd0 + lamd1 * tR&m + lamd2 * tR&n;
    linpredR = linpredsurv + linpredcrsR;
	cureprop = exp(min(&linpredCURE, 600))/(1+exp(min(&linpredCURE, 600)));

	STL = exp(-exp(min(linpred, 600)));
	STR = exp(-exp(min(linpredR, 600)));
	logevent = log(cureprop) + log(STL - STR);
	logcensor = log((1 - cureprop) + cureprop * STL);

    if &SURVevent = 1 then loglike = logevent;
	if &SURVevent = 0 then loglike = logcensor;

	model &SURVevent ~ general(loglike);

    ods output ParameterEstimates=temp35;
    ods output FitStatistics = temp36; 
    ods output ConvergenceStatus = temp37 (keep = replicate status rename = (status = value));
    by replicate;
  run;

  %if %sysfunc(exist(temp36)) %then %do;

    data temp38;
      set temp36 temp37 (in = a);
      if a then Descr = 'Con';
    run;

    proc sort data = temp38;
      by replicate;
    run;

    proc transpose data = temp38 out = temp39 (rename = (AIC__smaller_is_better_ = AIC AICC__smaller_is_better_ = AICC BIC__smaller_is_better_ = BIC));
      id descr;
      var value;
      by replicate;
    run;

    data &out2;
      set temp39;
      drop _name_;
      if con = 0.0 then converge = 'Yes';
      else converge = 'No';
      drop con;
      method = 'FP';
      P1 = &m;
      P2 = &n;
      ID = &ID;
    run;

    data &out1;
      set temp35;
      method = 'FP';
      P1 = &m;
      P2 = &n;
      ID = &ID;
    run;

    proc datasets lib = work;
      delete temp35 temp36 temp37 temp38 temp39;
    run;

  %end;

  %else %do;

    data &out2;
      format converge $3.;
      converge = 'No';
      method = 'FP';
      P1 = &m;
      P2 = &n;
      ID = &ID;
    run;

    data &out1;
      method = 'FP';
      P1 = &m;
      P2 = &n;
      ID = &ID;
    run;

  %end;

  proc datasets lib = work;
    delete temp31 temp311 temp312 temp313 temp314;
  run;

%mend;

** EXP is sub Macro to calculate the mean potential outcomes used for the exposure effect estimation under the potential outcome framework. In the EXP sub Macro, it also called 16 
sub sub Macros (INTE11, INTE211, INTE212, INTE213, INTE214, INTE311, INTE312, INTE313, INTE314, INTE315, INTE411, INTE412, INTE413, INTE414, INTE415 and INTE416) to calculate the 
mean potential outcomes for FP baseline hazard function (1 submcro), RCS with degree of freedom 2 (4 submacros), RCS with degree of freedom 3 (5 submacros) and RCS with degree of 
freedom 4 (6 submacros) baseline hazard function respectively.;

/* FP Potential Outcome Mean Calculation */

%macro inte11(indata, outdata);

  data temp62;
    set &indata;
  run;

  proc iml;
    use temp60r;
    read all into surv;
    use temp61r;
    read all into baseline;

	%if &stat = 1 | &stat = 2 %then %do;
	  %put 1;
      use temp62;
      read all var{id} into subjid;
      replicate = subjid[,1];
      obs = nrow(replicate);
      read all var{&SURVCov} into cov;

	  surv1 = j(obs,1,1)|| cov;
      surv0 = j(obs,1,0)|| cov;
    %end;

    %else %if &stat = 3 | &stat = 4 %then %do;
	  %put 2;
      use temp62;
      read all var{id} into subjid;
      replicate = subjid[,1];
      obs = nrow(replicate);

	  surv1 = j(obs,1,1);
      surv0 = j(obs,1,0);
    %end;

    alpha = surv[1, 1:ncol(surv)];

    if &id >= 1 & &id <= 8 then do;
      lamd0 = baseline[1, 1];
      lamd1 = baseline[1, 2];
    end;

    if &id > 8 & &id <= 44 then do;
      lamd0 = baseline[1, 1];
      lamd1 = baseline[1, 2];
      lamd2 = baseline[1, 3];
    end;

    out0=j(obs,1,.);
    out1=j(obs,1,.);
    outt0=j(obs,1,.);
    outt1=j(obs,1,.);

/* Restricted Mean Survival Time(tmax) potential outcome (1)*/

    start fun1(s) global(alpha, sur1, lamd0, lamd1, lamd2); 
      n = min(((&sx) + sur1*t(alpha)), 600);
      v = exp((-1)*exp(n));
      return(v); 
    finish; 

/* Restricted Mean Survival Time(tmax) potential outcome (0)*/

    start fun0(s) global(alpha, sur0, lamd0, lamd1, lamd2); 
      n = min(((&sx) + sur0*t(alpha)), 600);
      v = exp((-1)*exp(n));
      return(v); 
    finish; 

/* Survival Probability (t*) potential outcome (1)*/

    start funn1(s) global(alpha, sur1, lamd0, lamd1, lamd2); 
      n = min(((&sx) + sur1*t(alpha)), 600);
      v = exp((-1)*exp(n));
      return(v); 
    finish; 

/* Survival Probability (t*) potential outcome (0)*/

    start funn0(s) global(alpha, sur0, lamd0, lamd1, lamd2); 
      n = min(((&sx) + sur0*t(alpha)), 600);
      v = exp((-1)*exp(n));
      return(v); 
    finish; 

    do sub = 1 to obs;
      sur1 = surv1[sub,];
      sur0 = surv0[sub,];
      a = { 0  &rmst };
	  outt0[sub] = funn0(&rmst);
	  outt1[sub] = funn1(&rmst);
      call quad(z0, "fun0", a);
      out0[sub] = z0;
      call quad(z1, "fun1", a);
      out1[sub] = z1;
    end;

    postprobs=subjid||out1||out0||outt1||outt0;
    cname = {"id" "RMST1" "RMST0" "S1" "S0"};

    create &outdata from postprobs  [ colname=cname ];
    append from postprobs; 

    quit;

%mend;

/* RCS Mean Potential Outcome Calculation */

/* RCS DF = 2*/

/* Assessed t* <= smin */

%macro inte211(indata, outdata);

  data temp72;
    set &indata;
  run;

  proc iml;
    use temp60r;
    read all into surv;
    use temp61r;
    read all into baseline;
    use temp13r;
    read all into rcs;

	%if &stat = 1 | &stat = 2 %then %do;
	  %put 1;
	  use temp72;
      read all var{id} into subjid;
      replicate = subjid[,1];
      obs = nrow(replicate);
      read all var{&SURVCov} into cov;

      surv1 = j(obs,1,1)|| cov;
      surv0 = j(obs,1,0)|| cov;
	%end;

    %else %if &stat = 3 | &stat = 4 %then %do;
	  %put 2;
	  use temp72;
      read all var{id} into subjid;
      replicate = subjid[,1];
      obs = nrow(replicate);

      surv1 = j(obs,1,1);
      surv0 = j(obs,1,0);
	%end;

    alpha = surv[1, 1:ncol(surv)];

    gam10 = baseline[1, 1];
    gam11 = baseline[1, 2];
    gam12 = baseline[1, 3];

    k_min = rcs[1, 1];
    k_max = rcs[1, 2];
    lamda11 = rcs[1, 3];
    k11 = rcs[1, 9];

    outt0=j(obs,1,.);
    outt1=j(obs,1,.);
    outcome10=j(obs,1,.);
    outcome20=j(obs,1,.);
    outcome30=j(obs,1,.);
    outcome11=j(obs,1,.);
    outcome21=j(obs,1,.);
    outcome31=j(obs,1,.);

    start fun10(s) global(alpha, sur0, gam10, gam11, gam12, k_min, k_max, lamda11, k11); 
      n = min(((&sx1) + sur0*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun11(s) global(alpha, sur1, gam10, gam11, gam12, k_min, k_max, lamda11, k11); 
      n = min(((&sx1) + sur1*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

/* Survival Probability (t*) potential outcome (1)*/

    start funn11(s) global(alpha, sur1, gam10, gam11, gam12, k_min, k_max, lamda11, k11); 
      n = min(((&sx1) + sur1*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

/* Survival Probability (t*) potential outcome (0)*/

    start funn10(s) global(alpha, sur0, gam10, gam11, gam12, k_min, k_max, lamda11, k11); 
      n = min(((&sx1) + sur0*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    do sub = 1 to obs;
      sur1 = surv1[sub,];
      sur0 = surv0[sub,];
      a1 = { 0  &rmst };

	  outt0[sub] = funn10(&rmst);
	  outt1[sub] = funn11(&rmst);
      call quad(z10,"fun10",a1); 
      outcome10[sub] = z10;
	  call quad(z11,"fun11",a1); 
      outcome11[sub] = z11;

    end;

    out0 = outcome10;
    out1 = outcome11;

    postprobs=subjid||out1||out0||outt1||outt0;
    cname = {"id" "RMST1" "RMST0" "S1" "S0"};

    create &outdata from postprobs  [ colname=cname ];
    append from postprobs; 
  quit;

%mend;

/* Assessed smin < t* <= s11 */

%macro inte212(indata, outdata);

  data temp72;
    set &indata;
  run;

  proc iml;
    use temp60r;
    read all into surv;
    use temp61r;
    read all into baseline;
    use temp13r;
    read all into rcs;

	%if &stat = 1 | &stat = 2 %then %do;
	  %put 1;
	  use temp72;
      read all var{id} into subjid;
      replicate = subjid[,1];
      obs = nrow(replicate);
      read all var{&SURVCov} into cov;

      surv1 = j(obs,1,1)|| cov;
      surv0 = j(obs,1,0)|| cov;
	%end;

    %else %if &stat = 3 | &stat = 4 %then %do;
	  %put 2;
	  use temp72;
      read all var{id} into subjid;
      replicate = subjid[,1];
      obs = nrow(replicate);

      surv1 = j(obs,1,1);
      surv0 = j(obs,1,0);
	%end;

    alpha = surv[1, 1:ncol(surv)];

    gam10 = baseline[1, 1];
    gam11 = baseline[1, 2];
    gam12 = baseline[1, 3];

    k_min = rcs[1, 1];
    k_max = rcs[1, 2];
    lamda11 = rcs[1, 3];
    k11 = rcs[1, 9];

    outt0=j(obs,1,.);
    outt1=j(obs,1,.);
    outcome10=j(obs,1,.);
    outcome20=j(obs,1,.);
    outcome30=j(obs,1,.);
    outcome11=j(obs,1,.);
    outcome21=j(obs,1,.);
    outcome31=j(obs,1,.);

    start fun10(s) global(alpha, sur0, gam10, gam11, gam12, k_min, k_max, lamda11, k11); 
      n = min(((&sx1) + sur0*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun20(s) global(alpha, sur0, gam10, gam11, gam12, k_min, k_max, lamda11, k11); 
      n = min(((&sx2) + sur0*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun11(s) global(alpha, sur1, gam10, gam11, gam12, k_min, k_max, lamda11, k11); 
      n = min(((&sx1) + sur1*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun21(s) global(alpha, sur1, gam10, gam11, gam12, k_min, k_max, lamda11, k11); 
      n = min(((&sx2) + sur1*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

/* Survival Probability (t*) potential outcome (1)*/

    start funn21(s) global(alpha, sur1, gam10, gam11, gam12, k_min, k_max, lamda11, k11); 
      n = min(((&sx2) + sur1*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

/* Survival Probability (t*) potential outcome (0)*/

    start funn20(s) global(alpha, sur0, gam10, gam11, gam12, k_min, k_max, lamda11, k11); 
      n = min(((&sx2) + sur0*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    do sub = 1 to obs;
      sur1 = surv1[sub,];
      sur0 = surv0[sub,];
      a1 = { 0  &smin };
      a2 = { &smin &rmst};

	  outt0[sub] = funn20(&rmst);
	  outt1[sub] = funn21(&rmst);
      call quad(z10,"fun10",a1); 
      outcome10[sub] = z10;
      call quad(z20,"fun20",a2); 
      outcome20[sub] = z20;
	  call quad(z11,"fun11",a1); 
      outcome11[sub] = z11;
      call quad(z21,"fun21",a2); 
      outcome21[sub] = z21;

    end;

    out0 = outcome10 + outcome20;
    out1 = outcome11 + outcome21;

    postprobs=subjid||out1||out0||outt1||outt0;
    cname = {"id" "RMST1" "RMST0" "S1" "S0"};

    create &outdata from postprobs  [ colname=cname ];
    append from postprobs; 
  quit;

%mend;

/* Assessed s11 < t* <= smax */

%macro inte213(indata, outdata);

  data temp72;
    set &indata;
  run;

  proc iml;
    use temp60r;
    read all into surv;
    use temp61r;
    read all into baseline;
    use temp13r;
    read all into rcs;

	%if &stat = 1 | &stat = 2 %then %do;
	  %put 1;
	  use temp72;
      read all var{id} into subjid;
      replicate = subjid[,1];
      obs = nrow(replicate);
      read all var{&SURVCov} into cov;

      surv1 = j(obs,1,1)|| cov;
      surv0 = j(obs,1,0)|| cov;
	%end;

    %else %if &stat = 3 | &stat = 4 %then %do;
	  %put 2;
	  use temp72;
      read all var{id} into subjid;
      replicate = subjid[,1];
      obs = nrow(replicate);

      surv1 = j(obs,1,1);
      surv0 = j(obs,1,0);
	%end;

    alpha = surv[1, 1:ncol(surv)];

    gam10 = baseline[1, 1];
    gam11 = baseline[1, 2];
    gam12 = baseline[1, 3];

    k_min = rcs[1, 1];
    k_max = rcs[1, 2];
    lamda11 = rcs[1, 3];
    k11 = rcs[1, 9];

    outt0=j(obs,1,.);
    outt1=j(obs,1,.);
    outcome10=j(obs,1,.);
    outcome20=j(obs,1,.);
    outcome30=j(obs,1,.);
    outcome11=j(obs,1,.);
    outcome21=j(obs,1,.);
    outcome31=j(obs,1,.);

    start fun10(s) global(alpha, sur0, gam10, gam11, gam12, k_min, k_max, lamda11, k11); 
      n = min(((&sx1) + sur0*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun20(s) global(alpha, sur0, gam10, gam11, gam12, k_min, k_max, lamda11, k11); 
      n = min(((&sx2) + sur0*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun30(s) global(alpha, sur0, gam10, gam11, gam12, k_min, k_max, lamda11, k11); 
      n = min(((&sx3) + sur0*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    start fun11(s) global(alpha, sur1, gam10, gam11, gam12, k_min, k_max, lamda11, k11); 
      n = min(((&sx1) + sur1*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun21(s) global(alpha, sur1, gam10, gam11, gam12, k_min, k_max, lamda11, k11); 
      n = min(((&sx2) + sur1*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun31(s) global(alpha, sur1, gam10, gam11, gam12, k_min, k_max, lamda11, k11); 
      n = min(((&sx3) + sur1*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

/* Survival Probability (t*) potential outcome (1)*/

    start funn31(s) global(alpha, sur1, gam10, gam11, gam12, k_min, k_max, lamda11, k11); 
      n = min(((&sx3) + sur1*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

/* Survival Probability (t*) potential outcome (0)*/

	start funn30(s) global(alpha, sur0, gam10, gam11, gam12, k_min, k_max, lamda11, k11); 
      n = min(((&sx3) + sur0*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    do sub = 1 to obs;
      sur1 = surv1[sub,];
      sur0 = surv0[sub,];
      a1 = { 0  &smin };
      a2 = { &smin &s11};
      a3 = { &s11 &rmst };

	  outt0[sub] = funn30(&rmst);
	  outt1[sub] = funn31(&rmst);
      call quad(z10,"fun10",a1); 
      outcome10[sub] = z10;
      call quad(z20,"fun20",a2); 
      outcome20[sub] = z20;
      call quad(z30,"fun30",a3); 
      outcome30[sub] = z30;
	  call quad(z11,"fun11",a1); 
      outcome11[sub] = z11;
      call quad(z21,"fun21",a2); 
      outcome21[sub] = z21;
      call quad(z31,"fun31",a3); 
      outcome31[sub] = z31;

    end;

    out0 = outcome10 + outcome20 + outcome30;
    out1 = outcome11 + outcome21 + outcome31;

    postprobs=subjid||out1||out0||outt1||outt0;
    cname = {"id" "RMST1" "RMST0" "S1" "S0"};

    create &outdata from postprobs  [ colname=cname ];
    append from postprobs; 
  quit;

%mend;

/* Assessed t* > smax */

%macro inte214(indata, outdata);

  data temp72;
    set &indata;
  run;

  proc iml;
    use temp60r;
    read all into surv;
    use temp61r;
    read all into baseline;
    use temp13r;
    read all into rcs;

	%if &stat = 1 | &stat = 2 %then %do;
	  %put 1;
	  use temp72;
      read all var{id} into subjid;
      replicate = subjid[,1];
      obs = nrow(replicate);
      read all var{&SURVCov} into cov;

      surv1 = j(obs,1,1)|| cov;
      surv0 = j(obs,1,0)|| cov;
	%end;

    %else %if &stat = 3 | &stat = 4 %then %do;
	  %put 2;
	  use temp72;
      read all var{id} into subjid;
      replicate = subjid[,1];
      obs = nrow(replicate);

      surv1 = j(obs,1,1);
      surv0 = j(obs,1,0);
	%end;

    alpha = surv[1, 1:ncol(surv)];

    gam10 = baseline[1, 1];
    gam11 = baseline[1, 2];
    gam12 = baseline[1, 3];

    k_min = rcs[1, 1];
    k_max = rcs[1, 2];
    lamda11 = rcs[1, 3];
    k11 = rcs[1, 9];

    outt0=j(obs,1,.);
    outt1=j(obs,1,.);
    outcome10=j(obs,1,.);
    outcome20=j(obs,1,.);
    outcome30=j(obs,1,.);
    outcome40=j(obs,1,.);
    outcome11=j(obs,1,.);
    outcome21=j(obs,1,.);
    outcome31=j(obs,1,.);
    outcome41=j(obs,1,.);

    start fun10(s) global(alpha, sur0, gam10, gam11, gam12, k_min, k_max, lamda11, k11); 
      n = min(((&sx1) + sur0*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun20(s) global(alpha, sur0, gam10, gam11, gam12, k_min, k_max, lamda11, k11); 
      n = min(((&sx2) + sur0*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun30(s) global(alpha, sur0, gam10, gam11, gam12, k_min, k_max, lamda11, k11); 
      n = min(((&sx3) + sur0*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

	start fun40(s) global(alpha, sur0, gam10, gam11, gam12, k_min, k_max, lamda11, k11); 
      n = min(((&sx4) + sur0*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    start fun11(s) global(alpha, sur1, gam10, gam11, gam12, k_min, k_max, lamda11, k11); 
      n = min(((&sx1) + sur1*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun21(s) global(alpha, sur1, gam10, gam11, gam12, k_min, k_max, lamda11, k11); 
      n = min(((&sx2) + sur1*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun31(s) global(alpha, sur1, gam10, gam11, gam12, k_min, k_max, lamda11, k11); 
      n = min(((&sx3) + sur1*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    start fun41(s) global(alpha, sur1, gam10, gam11, gam12, k_min, k_max, lamda11, k11); 
      n = min(((&sx4) + sur1*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

/* Survival Probability (t*) potential outcome (1)*/

    start funn41(s) global(alpha, sur1, gam10, gam11, gam12, k_min, k_max, lamda11, k11); 
      n = min(((&sx4) + sur1*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

/* Survival Probability (t*) potential outcome (0)*/

	start funn40(s) global(alpha, sur0, gam10, gam11, gam12, k_min, k_max, lamda11, k11); 
      n = min(((&sx4) + sur0*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    do sub = 1 to obs;
      sur1 = surv1[sub,];
      sur0 = surv0[sub,];
      a1 = { 0  &smin };
      a2 = { &smin &s11};
      a3 = { &s11 &smax };
      a4 = { &smax &rmst };

	  outt0[sub] = funn40(&rmst);
	  outt1[sub] = funn41(&rmst);
      call quad(z10,"fun10",a1); 
      outcome10[sub] = z10;
      call quad(z20,"fun20",a2); 
      outcome20[sub] = z20;
      call quad(z30,"fun30",a3); 
      outcome30[sub] = z30;
      call quad(z40,"fun40",a4); 
      outcome40[sub] = z40;
	  call quad(z11,"fun11",a1); 
      outcome11[sub] = z11;
      call quad(z21,"fun21",a2); 
      outcome21[sub] = z21;
      call quad(z31,"fun31",a3); 
      outcome31[sub] = z31;
      call quad(z41,"fun41",a4); 
      outcome41[sub] = z41;

    end;

    out0 = outcome10 + outcome20 + outcome30 + outcome40;
    out1 = outcome11 + outcome21 + outcome31 + outcome41;

    postprobs=subjid||out1||out0||outt1||outt0;
    cname = {"id" "RMST1" "RMST0" "S1" "S0"};

    create &outdata from postprobs  [ colname=cname ];
    append from postprobs; 
  quit;

%mend;

/* RCS DF = 3*/

/* Assessed t* <= smin */

%macro inte311(indata, outdata);

  data temp82;
    set &indata;
  run;

  proc iml;
    use temp60r;
    read all into surv;
    use temp61r;
    read all into baseline;
    use temp13r;
    read all into rcs;

	%if &stat = 1 | &stat = 2 %then %do;
	  %put 1;
	  use temp82;
      read all var{id} into subjid;
      replicate = subjid[,1];
      obs = nrow(replicate);
      read all var{&SURVCov} into cov;

      surv1 = j(obs,1,1)|| cov;
      surv0 = j(obs,1,0)|| cov;
	%end;

    %else %if &stat = 3 | &stat = 4 %then %do;
	  %put 2;
	  use temp82;
      read all var{id} into subjid;
      replicate = subjid[,1];
      obs = nrow(replicate);

      surv1 = j(obs,1,1);
      surv0 = j(obs,1,0);
	%end;

    alpha = surv[1, 1:ncol(surv)];

    gam20 = baseline[1, 1];
    gam21 = baseline[1, 2];
    gam22 = baseline[1, 3];
    gam23 = baseline[1, 4];

    k_min = rcs[1, 1];
    k_max = rcs[1, 2];
    lamda21 = rcs[1, 4];
    lamda22 = rcs[1, 5];
	k21 = rcs[1, 10];
    k22 = rcs[1, 11];

    outt0=j(obs,1,.);
    outt1=j(obs,1,.);
    outcome10=j(obs,1,.);
    outcome20=j(obs,1,.);
    outcome30=j(obs,1,.);
    outcome40=j(obs,1,.);
    outcome11=j(obs,1,.);
    outcome21=j(obs,1,.);
    outcome31=j(obs,1,.);
    outcome41=j(obs,1,.);

    start fun10(s) global(alpha, sur0, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx1) + sur0*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun11(s) global(alpha, sur1, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx1) + sur1*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

/* Survival Probability (t*) potential outcome (1)*/

    start funn11(s) global(alpha, sur1, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx1) + sur1*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

/* Survival Probability (t*) potential outcome (0)*/

    start funn10(s) global(alpha, sur0, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx1) + sur0*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    do sub = 1 to obs;
      sur1 = surv1[sub,];
      sur0 = surv0[sub,];
      a1 = { 0  &rmst }; 

  	  outt0[sub] = funn10(&rmst);
	  outt1[sub] = funn11(&rmst);
      call quad(z10,"fun10",a1); 
      outcome10[sub] = z10;
	  call quad(z11,"fun11",a1); 
      outcome11[sub] = z11;
 
    end;

    out0 = outcome10;
    out1 = outcome11;

    postprobs=subjid||out1||out0||outt1||outt0;
    cname = {"id" "RMST1" "RMST0" "S1" "S0"};

    create &outdata from postprobs  [ colname=cname ];
    append from postprobs; 
  quit;

%mend;

/* Assessed smin < t* <= s21 */

%macro inte312(indata, outdata);

  data temp82;
    set &indata;
  run;

 proc iml;
    use temp60r;
    read all into surv;
    use temp61r;
    read all into baseline;
    use temp13r;
    read all into rcs;

	%if &stat = 1 | &stat = 2 %then %do;
	  %put 1;
	  use temp82;
      read all var{id} into subjid;
      replicate = subjid[,1];
      obs = nrow(replicate);
      read all var{&SURVCov} into cov;

      surv1 = j(obs,1,1)|| cov;
      surv0 = j(obs,1,0)|| cov;
	%end;

    %else %if &stat = 3 | &stat = 4 %then %do;
	  %put 2;
	  use temp82;
      read all var{id} into subjid;
      replicate = subjid[,1];
      obs = nrow(replicate);

      surv1 = j(obs,1,1);
      surv0 = j(obs,1,0);
	%end;

    alpha = surv[1, 1:ncol(surv)];

    gam20 = baseline[1, 1];
    gam21 = baseline[1, 2];
    gam22 = baseline[1, 3];
    gam23 = baseline[1, 4];

    k_min = rcs[1, 1];
    k_max = rcs[1, 2];
    lamda21 = rcs[1, 4];
    lamda22 = rcs[1, 5];
	k21 = rcs[1, 10];
    k22 = rcs[1, 11];

    outt0=j(obs,1,.);
    outt1=j(obs,1,.);
    outcome10=j(obs,1,.);
    outcome20=j(obs,1,.);
    outcome30=j(obs,1,.);
    outcome40=j(obs,1,.);
    outcome11=j(obs,1,.);
    outcome21=j(obs,1,.);
    outcome31=j(obs,1,.);
    outcome41=j(obs,1,.);

    start fun10(s) global(alpha, sur0, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx1) + sur0*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun20(s) global(alpha, sur0, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx2) + sur0*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun11(s) global(alpha, sur1, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx1) + sur1*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun21(s) global(alpha, sur1, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx2) + sur1*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

/* Survival Probability (t*) potential outcome (1)*/

    start funn21(s) global(alpha, sur1, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx2) + sur1*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

/* Survival Probability (t*) potential outcome (0)*/

    start funn20(s) global(alpha, sur0, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx2) + sur0*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    do sub = 1 to obs;
      sur1 = surv1[sub,];
      sur0 = surv0[sub,];
      a1 = { 0  &smin }; 
      a2 = { &smin &rmst }; 

  	  outt0[sub] = funn20(&rmst);
	  outt1[sub] = funn21(&rmst);
      call quad(z10,"fun10",a1); 
      outcome10[sub] = z10;
      call quad(z20,"fun20",a2); 
      outcome20[sub] = z20;
	  call quad(z11,"fun11",a1); 
      outcome11[sub] = z11;
      call quad(z21,"fun21",a2); 
      outcome21[sub] = z21;

    end;

    out0 = outcome10 + outcome20;
    out1 = outcome11 + outcome21;

    postprobs=subjid||out1||out0||outt1||outt0;
    cname = {"id" "RMST1" "RMST0" "S1" "S0"};

    create &outdata from postprobs  [ colname=cname ];
    append from postprobs; 
  quit;

%mend;

/* Assessed s21 < t* <= s22 */

%macro inte313(indata, outdata);

  data temp82;
    set &indata;
  run;

 proc iml;
    use temp60r;
    read all into surv;
    use temp61r;
    read all into baseline;
    use temp13r;
    read all into rcs;

	%if &stat = 1 | &stat = 2 %then %do;
	  %put 1;
	  use temp82;
      read all var{id} into subjid;
      replicate = subjid[,1];
      obs = nrow(replicate);
      read all var{&SURVCov} into cov;

      surv1 = j(obs,1,1)|| cov;
      surv0 = j(obs,1,0)|| cov;
	%end;

    %else %if &stat = 3 | &stat = 4 %then %do;
	  %put 2;
	  use temp82;
      read all var{id} into subjid;
      replicate = subjid[,1];
      obs = nrow(replicate);

      surv1 = j(obs,1,1);
      surv0 = j(obs,1,0);
	%end;

    alpha = surv[1, 1:ncol(surv)];

    gam20 = baseline[1, 1];
    gam21 = baseline[1, 2];
    gam22 = baseline[1, 3];
    gam23 = baseline[1, 4];

    k_min = rcs[1, 1];
    k_max = rcs[1, 2];
    lamda21 = rcs[1, 4];
    lamda22 = rcs[1, 5];
	k21 = rcs[1, 10];
    k22 = rcs[1, 11];

    outt0=j(obs,1,.);
    outt1=j(obs,1,.);
    outcome10=j(obs,1,.);
    outcome20=j(obs,1,.);
    outcome30=j(obs,1,.);
    outcome40=j(obs,1,.);
    outcome11=j(obs,1,.);
    outcome21=j(obs,1,.);
    outcome31=j(obs,1,.);
    outcome41=j(obs,1,.);

    start fun10(s) global(alpha, sur0, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx1) + sur0*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun20(s) global(alpha, sur0, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx2) + sur0*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun30(s) global(alpha, sur0, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx3) + sur0*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    start fun11(s) global(alpha, sur1, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx1) + sur1*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun21(s) global(alpha, sur1, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx2) + sur1*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun31(s) global(alpha, sur1, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx3) + sur1*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

/* Survival Probability (t*) potential outcome (1)*/

    start funn31(s) global(alpha, sur1, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx3) + sur1*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

/* Survival Probability (t*) potential outcome (0)*/

    start funn30(s) global(alpha, sur0, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx3) + sur0*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    do sub = 1 to obs;
      sur1 = surv1[sub,];
      sur0 = surv0[sub,];
      a1 = { 0  &smin }; 
      a2 = { &smin &s21 }; 
      a3 = { &s21 &rmst }; 

  	  outt0[sub] = funn30(&rmst);
	  outt1[sub] = funn31(&rmst);
      call quad(z10,"fun10",a1); 
      outcome10[sub] = z10;
      call quad(z20,"fun20",a2); 
      outcome20[sub] = z20;
      call quad(z30,"fun30",a3); 
      outcome30[sub] = z30;
	  call quad(z11,"fun11",a1); 
      outcome11[sub] = z11;
      call quad(z21,"fun21",a2); 
      outcome21[sub] = z21;
      call quad(z31,"fun31",a3); 
      outcome31[sub] = z31;

    end;

    out0 = outcome10 + outcome20 + outcome30;
    out1 = outcome11 + outcome21 + outcome31;

    postprobs=subjid||out1||out0||outt1||outt0;
    cname = {"id" "RMST1" "RMST0" "S1" "S0"};

    create &outdata from postprobs  [ colname=cname ];
    append from postprobs; 
  quit;

%mend;

/* Assessed s22 < t* <= smax */

%macro inte314(indata, outdata);

  data temp82;
    set &indata;
  run;

 proc iml;
    use temp60r;
    read all into surv;
    use temp61r;
    read all into baseline;
    use temp13r;
    read all into rcs;

	%if &stat = 1 | &stat = 2 %then %do;
	  %put 1;
	  use temp82;
      read all var{id} into subjid;
      replicate = subjid[,1];
      obs = nrow(replicate);
      read all var{&SURVCov} into cov;

      surv1 = j(obs,1,1)|| cov;
      surv0 = j(obs,1,0)|| cov;
	%end;

    %else %if &stat = 3 | &stat = 4 %then %do;
	  %put 2;
	  use temp82;
      read all var{id} into subjid;
      replicate = subjid[,1];
      obs = nrow(replicate);

      surv1 = j(obs,1,1);
      surv0 = j(obs,1,0);
	%end;

    alpha = surv[1, 1:ncol(surv)];

    gam20 = baseline[1, 1];
    gam21 = baseline[1, 2];
    gam22 = baseline[1, 3];
    gam23 = baseline[1, 4];

    k_min = rcs[1, 1];
    k_max = rcs[1, 2];
    lamda21 = rcs[1, 4];
    lamda22 = rcs[1, 5];
	k21 = rcs[1, 10];
    k22 = rcs[1, 11];

    outt0=j(obs,1,.);
    outt1=j(obs,1,.);
    outcome10=j(obs,1,.);
    outcome20=j(obs,1,.);
    outcome30=j(obs,1,.);
    outcome40=j(obs,1,.);
    outcome11=j(obs,1,.);
    outcome21=j(obs,1,.);
    outcome31=j(obs,1,.);
    outcome41=j(obs,1,.);

    start fun10(s) global(alpha, sur0, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx1) + sur0*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun20(s) global(alpha, sur0, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx2) + sur0*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun30(s) global(alpha, sur0, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx3) + sur0*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    start fun40(s) global(alpha, sur0, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx4) + sur0*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    start fun11(s) global(alpha, sur1, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx1) + sur1*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun21(s) global(alpha, sur1, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx2) + sur1*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun31(s) global(alpha, sur1, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx3) + sur1*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    start fun41(s) global(alpha, sur1, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx4) + sur1*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

/* Survival Probability (t*) potential outcome (1)*/

    start funn41(s) global(alpha, sur1, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx4) + sur1*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

/* Survival Probability (t*) potential outcome (0)*/

    start funn40(s) global(alpha, sur0, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx4) + sur0*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    do sub = 1 to obs;
      sur1 = surv1[sub,];
      sur0 = surv0[sub,];
      a1 = { 0  &smin }; 
      a2 = { &smin &s21 }; 
      a3 = { &s21 &s22 }; 
      a4 = { &s22  &rmst }; 

  	  outt0[sub] = funn40(&rmst);
	  outt1[sub] = funn41(&rmst);
      call quad(z10,"fun10",a1); 
      outcome10[sub] = z10;
      call quad(z20,"fun20",a2); 
      outcome20[sub] = z20;
      call quad(z30,"fun30",a3); 
      outcome30[sub] = z30;
      call quad(z40,"fun40",a4); 
      outcome40[sub] = z40;
	  call quad(z11,"fun11",a1); 
      outcome11[sub] = z11;
      call quad(z21,"fun21",a2); 
      outcome21[sub] = z21;
      call quad(z31,"fun31",a3); 
      outcome31[sub] = z31;
      call quad(z41,"fun41",a4); 
      outcome41[sub] = z41;

    end;

    out0 = outcome10 + outcome20 + outcome30 + outcome40;
    out1 = outcome11 + outcome21 + outcome31 + outcome41;

    postprobs=subjid||out1||out0||outt1||outt0;
    cname = {"id" "RMST1" "RMST0" "S1" "S0"};

    create &outdata from postprobs  [ colname=cname ];
    append from postprobs; 
  quit;

%mend;

/* Assessed t* > smax */

%macro inte315(indata, outdata);

  data temp82;
    set &indata;
  run;

 proc iml;
    use temp60r;
    read all into surv;
    use temp61r;
    read all into baseline;
    use temp13r;
    read all into rcs;

	%if &stat = 1 | &stat = 2 %then %do;
	  %put 1;
	  use temp82;
      read all var{id} into subjid;
      replicate = subjid[,1];
      obs = nrow(replicate);
      read all var{&SURVCov} into cov;

      surv1 = j(obs,1,1)|| cov;
      surv0 = j(obs,1,0)|| cov;
	%end;

    %else %if &stat = 3 | &stat = 4 %then %do;
	  %put 2;
	  use temp82;
      read all var{id} into subjid;
      replicate = subjid[,1];
      obs = nrow(replicate);

      surv1 = j(obs,1,1);
      surv0 = j(obs,1,0);
	%end;

    alpha = surv[1, 1:ncol(surv)];

    gam20 = baseline[1, 1];
    gam21 = baseline[1, 2];
    gam22 = baseline[1, 3];
    gam23 = baseline[1, 4];

    k_min = rcs[1, 1];
    k_max = rcs[1, 2];
    lamda21 = rcs[1, 4];
    lamda22 = rcs[1, 5];
	k21 = rcs[1, 10];
    k22 = rcs[1, 11];

    outt0=j(obs,1,.);
    outt1=j(obs,1,.);
    outcome10=j(obs,1,.);
    outcome20=j(obs,1,.);
    outcome30=j(obs,1,.);
    outcome40=j(obs,1,.);
    outcome50=j(obs,1,.);
    outcome11=j(obs,1,.);
    outcome21=j(obs,1,.);
    outcome31=j(obs,1,.);
    outcome41=j(obs,1,.);
    outcome51=j(obs,1,.);

    start fun10(s) global(alpha, sur0, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx1) + sur0*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun20(s) global(alpha, sur0, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx2) + sur0*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun30(s) global(alpha, sur0, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx3) + sur0*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    start fun40(s) global(alpha, sur0, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx4) + sur0*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    start fun50(s) global(alpha, sur0, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx5) + sur0*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    start fun11(s) global(alpha, sur1, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx1) + sur1*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun21(s) global(alpha, sur1, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx2) + sur1*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun31(s) global(alpha, sur1, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx3) + sur1*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    start fun41(s) global(alpha, sur1, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx4) + sur1*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    start fun51(s) global(alpha, sur1, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx5) + sur1*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

/* Survival Probability (t*) potential outcome (1)*/

    start funn51(s) global(alpha, sur1, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx5) + sur1*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

/* Survival Probability (t*) potential outcome (0)*/

    start funn50(s) global(alpha, sur0, gam20, gam21, gam22, gam23, k_min, k_max, lamda21, lamda22, k21, k22); 
      n = min(((&sx5) + sur0*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    do sub = 1 to obs;
      sur1 = surv1[sub,];
      sur0 = surv0[sub,];
      a1 = { 0  &smin }; 
      a2 = { &smin &s21 }; 
      a3 = { &s21 &s22 }; 
      a4 = { &s22  &smax }; 
      a5 = { &smax &rmst };

  	  outt0[sub] = funn50(&rmst);
	  outt1[sub] = funn51(&rmst);
      call quad(z10,"fun10",a1); 
      outcome10[sub] = z10;
      call quad(z20,"fun20",a2); 
      outcome20[sub] = z20;
      call quad(z30,"fun30",a3); 
      outcome30[sub] = z30;
      call quad(z40,"fun40",a4); 
      outcome40[sub] = z40;
      call quad(z50,"fun50",a5); 
      outcome50[sub] = z50;
	  call quad(z11,"fun11",a1); 
      outcome11[sub] = z11;
      call quad(z21,"fun21",a2); 
      outcome21[sub] = z21;
      call quad(z31,"fun31",a3); 
      outcome31[sub] = z31;
      call quad(z41,"fun41",a4); 
      outcome41[sub] = z41;
      call quad(z51,"fun51",a5); 
      outcome51[sub] = z51;

    end;

    out0 = outcome10 + outcome20 + outcome30 + outcome40 + outcome50;
    out1 = outcome11 + outcome21 + outcome31 + outcome41 + outcome51;

    postprobs=subjid||out1||out0||outt1||outt0;
    cname = {"id" "RMST1" "RMST0" "S1" "S0"};

    create &outdata from postprobs  [ colname=cname ];
    append from postprobs; 
  quit;

%mend;

/* RCS DF = 4*/

/* Assessed t* <= smin */

%macro inte411(indata, outdata);

  data temp92;
    set &indata;
  run;

 proc iml;
    use temp60r;
    read all into surv;
    use temp61r;
    read all into baseline;
    use temp13r;
    read all into rcs;

	%if &stat = 1 | &stat = 2 %then %do;
	  %put 1;
	  use temp92;
      read all var{id} into subjid;
      replicate = subjid[,1];
      obs = nrow(replicate);
      read all var{&SURVCov} into cov;

      surv1 = j(obs,1,1)|| cov;
      surv0 = j(obs,1,0)|| cov;
	%end;

    %else %if &stat = 3 | &stat = 4 %then %do;
	  %put 2;
	  use temp92;
      read all var{id} into subjid;
      replicate = subjid[,1];
      obs = nrow(replicate);

      surv1 = j(obs,1,1);
      surv0 = j(obs,1,0);
	%end;

    alpha = surv[1, 1:ncol(surv)];

    gam30 = baseline[1, 1];
    gam31 = baseline[1, 2];
    gam32 = baseline[1, 3];
    gam33 = baseline[1, 4];
    gam34 = baseline[1, 5];

    k_min = rcs[1, 1];
    k_max = rcs[1, 2];
    lamda31 = rcs[1, 6];
    lamda32 = rcs[1, 7];
    lamda33 = rcs[1, 8];

    k31 = rcs[1, 12];
    k32 = rcs[1, 13];
    k33 = rcs[1, 14];

    outt0=j(obs,1,.);
    outt1=j(obs,1,.);
    outcome10=j(obs,1,.);
    outcome20=j(obs,1,.);
    outcome30=j(obs,1,.);
    outcome40=j(obs,1,.);
    outcome50=j(obs,1,.);
    outcome11=j(obs,1,.);
    outcome21=j(obs,1,.);
    outcome31=j(obs,1,.);
    outcome41=j(obs,1,.);
    outcome51=j(obs,1,.);

    start fun10(s) global(alpha, sur0, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx1) + sur0*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun11(s) global(alpha, sur1, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx1) + sur1*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

/* Survival Probability (t*) potential outcome (1)*/

    start funn11(s) global(alpha, sur1, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx1) + sur1*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

/* Survival Probability (t*) potential outcome (0)*/

    start funn10(s) global(alpha, sur0, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx1) + sur0*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    do sub = 1 to obs;
      sur1 = surv1[sub,];
      sur0 = surv0[sub,];
      a1   = { 0  &rmst }; 

  	  outt0[sub] = funn10(&rmst);
	  outt1[sub] = funn11(&rmst);
      call quad(z10,"fun10",a1); 
      outcome10[sub] = z10;
	  call quad(z11,"fun11",a1); 
      outcome11[sub] = z11;

    end;

    out0 = outcome10;
    out1 = outcome11;

    postprobs=subjid||out1||out0||outt1||outt0;
    cname = {"id" "RMST1" "RMST0" "S1" "S0"};

    create &outdata from postprobs  [ colname=cname ];
    append from postprobs; 
  quit;

%mend;

/* Assessed smin < t* <= s31 */

%macro inte412(indata, outdata);

  data temp92;
    set &indata;
  run;

 proc iml;
    use temp60r;
    read all into surv;
    use temp61r;
    read all into baseline;
    use temp13r;
    read all into rcs;

	%if &stat = 1 | &stat = 2 %then %do;
	  %put 1;
	  use temp92;
      read all var{id} into subjid;
      replicate = subjid[,1];
      obs = nrow(replicate);
      read all var{&SURVCov} into cov;

      surv1 = j(obs,1,1)|| cov;
      surv0 = j(obs,1,0)|| cov;
	%end;

    %else %if &stat = 3 | &stat = 4 %then %do;
	  %put 2;
	  use temp92;
      read all var{id} into subjid;
      replicate = subjid[,1];
      obs = nrow(replicate);

      surv1 = j(obs,1,1);
      surv0 = j(obs,1,0);
	%end;

    alpha = surv[1, 1:ncol(surv)];

    gam30 = baseline[1, 1];
    gam31 = baseline[1, 2];
    gam32 = baseline[1, 3];
    gam33 = baseline[1, 4];
    gam34 = baseline[1, 5];

    k_min = rcs[1, 1];
    k_max = rcs[1, 2];
    lamda31 = rcs[1, 6];
    lamda32 = rcs[1, 7];
    lamda33 = rcs[1, 8];

    k31 = rcs[1, 12];
    k32 = rcs[1, 13];
    k33 = rcs[1, 14];

    outt0=j(obs,1,.);
    outt1=j(obs,1,.);
    outcome10=j(obs,1,.);
    outcome20=j(obs,1,.);
    outcome30=j(obs,1,.);
    outcome40=j(obs,1,.);
    outcome50=j(obs,1,.);
    outcome11=j(obs,1,.);
    outcome21=j(obs,1,.);
    outcome31=j(obs,1,.);
    outcome41=j(obs,1,.);
    outcome51=j(obs,1,.);

    start fun10(s) global(alpha, sur0, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx1) + sur0*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun20(s) global(alpha, sur0, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx2) + sur0*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun11(s) global(alpha, sur1, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx1) + sur1*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun21(s) global(alpha, sur1, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx2) + sur1*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

/* Survival Probability (t*) potential outcome (1)*/

    start funn21(s) global(alpha, sur1, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx2) + sur1*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

/* Survival Probability (t*) potential outcome (0)*/

    start funn20(s) global(alpha, sur0, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx2) + sur0*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    do sub = 1 to obs;
      sur1 = surv1[sub,];
      sur0 = surv0[sub,];
      a1   = { 0  &smin }; 
      a2   = { &smin &rmst }; 

  	  outt0[sub] = funn20(&rmst);
	  outt1[sub] = funn21(&rmst);
      call quad(z10,"fun10",a1); 
      outcome10[sub] = z10;
      call quad(z20,"fun20",a2); 
      outcome20[sub] = z20;
	  call quad(z11,"fun11",a1); 
      outcome11[sub] = z11;
      call quad(z21,"fun21",a2); 
      outcome21[sub] = z21;

    end;

    out0 = outcome10 + outcome20;
    out1 = outcome11 + outcome21;

    postprobs=subjid||out1||out0||outt1||outt0;
    cname = {"id" "RMST1" "RMST0" "S1" "S0"};

    create &outdata from postprobs  [ colname=cname ];
    append from postprobs; 
  quit;

%mend;

/* Assessed s31 < t* <= s32 */

%macro inte413(indata, outdata);

  data temp92;
    set &indata;
  run;

 proc iml;
    use temp60r;
    read all into surv;
    use temp61r;
    read all into baseline;
    use temp13r;
    read all into rcs;

	%if &stat = 1 | &stat = 2 %then %do;
	  %put 1;
	  use temp92;
      read all var{id} into subjid;
      replicate = subjid[,1];
      obs = nrow(replicate);
      read all var{&SURVCov} into cov;

      surv1 = j(obs,1,1)|| cov;
      surv0 = j(obs,1,0)|| cov;
	%end;

    %else %if &stat = 3 | &stat = 4 %then %do;
	  %put 2;
	  use temp92;
      read all var{id} into subjid;
      replicate = subjid[,1];
      obs = nrow(replicate);

      surv1 = j(obs,1,1);
      surv0 = j(obs,1,0);
	%end;

    alpha = surv[1, 1:ncol(surv)];

    gam30 = baseline[1, 1];
    gam31 = baseline[1, 2];
    gam32 = baseline[1, 3];
    gam33 = baseline[1, 4];
    gam34 = baseline[1, 5];

    k_min = rcs[1, 1];
    k_max = rcs[1, 2];
    lamda31 = rcs[1, 6];
    lamda32 = rcs[1, 7];
    lamda33 = rcs[1, 8];

    k31 = rcs[1, 12];
    k32 = rcs[1, 13];
    k33 = rcs[1, 14];

    outt0=j(obs,1,.);
    outt1=j(obs,1,.);
    outcome10=j(obs,1,.);
    outcome20=j(obs,1,.);
    outcome30=j(obs,1,.);
    outcome40=j(obs,1,.);
    outcome50=j(obs,1,.);
    outcome11=j(obs,1,.);
    outcome21=j(obs,1,.);
    outcome31=j(obs,1,.);
    outcome41=j(obs,1,.);
    outcome51=j(obs,1,.);

    start fun10(s) global(alpha, sur0, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx1) + sur0*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun20(s) global(alpha, sur0, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx2) + sur0*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun30(s) global(alpha, sur0, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx3) + sur0*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    start fun11(s) global(alpha, sur1, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx1) + sur1*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun21(s) global(alpha, sur1, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx2) + sur1*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun31(s) global(alpha, sur1, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx3) + sur1*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

/* Survival Probability (t*) potential outcome (1)*/

    start funn31(s) global(alpha, sur1, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx3) + sur1*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

/* Survival Probability (t*) potential outcome (0)*/

    start funn30(s) global(alpha, sur0, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx3) + sur0*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    do sub = 1 to obs;
      sur1 = surv1[sub,];
      sur0 = surv0[sub,];
      a1   = { 0  &smin }; 
      a2   = { &smin &s31 }; 
      a3   = { &s31 &rmst }; 

  	  outt0[sub] = funn30(&rmst);
	  outt1[sub] = funn31(&rmst);
      call quad(z10,"fun10",a1); 
      outcome10[sub] = z10;
      call quad(z20,"fun20",a2); 
      outcome20[sub] = z20;
      call quad(z30,"fun30",a3); 
      outcome30[sub] = z30;
	  call quad(z11,"fun11",a1); 
      outcome11[sub] = z11;
      call quad(z21,"fun21",a2); 
      outcome21[sub] = z21;
      call quad(z31,"fun31",a3); 
      outcome31[sub] = z31;

    end;

    out0 = outcome10 + outcome20 + outcome30;
    out1 = outcome11 + outcome21 + outcome31;

    postprobs=subjid||out1||out0||outt1||outt0;
    cname = {"id" "RMST1" "RMST0" "S1" "S0"};

    create &outdata from postprobs  [ colname=cname ];
    append from postprobs; 
  quit;

%mend;

/* Assessed s32 < t* <= s33 */

%macro inte414(indata, outdata);

  data temp92;
    set &indata;
  run;

 proc iml;
    use temp60r;
    read all into surv;
    use temp61r;
    read all into baseline;
    use temp13r;
    read all into rcs;

	%if &stat = 1 | &stat = 2 %then %do;
	  %put 1;
	  use temp92;
      read all var{id} into subjid;
      replicate = subjid[,1];
      obs = nrow(replicate);
      read all var{&SURVCov} into cov;

      surv1 = j(obs,1,1)|| cov;
      surv0 = j(obs,1,0)|| cov;
	%end;

    %else %if &stat = 3 | &stat = 4 %then %do;
	  %put 2;
	  use temp92;
      read all var{id} into subjid;
      replicate = subjid[,1];
      obs = nrow(replicate);

      surv1 = j(obs,1,1);
      surv0 = j(obs,1,0);
	%end;

    alpha = surv[1, 1:ncol(surv)];

    gam30 = baseline[1, 1];
    gam31 = baseline[1, 2];
    gam32 = baseline[1, 3];
    gam33 = baseline[1, 4];
    gam34 = baseline[1, 5];

    k_min = rcs[1, 1];
    k_max = rcs[1, 2];
    lamda31 = rcs[1, 6];
    lamda32 = rcs[1, 7];
    lamda33 = rcs[1, 8];

    k31 = rcs[1, 12];
    k32 = rcs[1, 13];
    k33 = rcs[1, 14];

    outt0=j(obs,1,.);
    outt1=j(obs,1,.);
    outcome10=j(obs,1,.);
    outcome20=j(obs,1,.);
    outcome30=j(obs,1,.);
    outcome40=j(obs,1,.);
    outcome50=j(obs,1,.);
    outcome11=j(obs,1,.);
    outcome21=j(obs,1,.);
    outcome31=j(obs,1,.);
    outcome41=j(obs,1,.);
    outcome51=j(obs,1,.);

    start fun10(s) global(alpha, sur0, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx1) + sur0*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun20(s) global(alpha, sur0, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx2) + sur0*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun30(s) global(alpha, sur0, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx3) + sur0*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    start fun40(s) global(alpha, sur0, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx4) + sur0*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    start fun11(s) global(alpha, sur1, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx1) + sur1*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun21(s) global(alpha, sur1, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx2) + sur1*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun31(s) global(alpha, sur1, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx3) + sur1*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    start fun41(s) global(alpha, sur1, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx4) + sur1*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

/* Survival Probability (t*) potential outcome (1)*/

    start funn41(s) global(alpha, sur1, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx4) + sur1*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

/* Survival Probability (t*) potential outcome (0)*/

    start funn40(s) global(alpha, sur0, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx4) + sur0*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    do sub = 1 to obs;
      sur1 = surv1[sub,];
      sur0 = surv0[sub,];
      a1   = { 0  &smin }; 
      a2   = { &smin &s31 }; 
      a3   = { &s31 &s32 }; 
      a4   = { &s32 &rmst }; 

  	  outt0[sub] = funn40(&rmst);
	  outt1[sub] = funn41(&rmst);
      call quad(z10,"fun10",a1); 
      outcome10[sub] = z10;
      call quad(z20,"fun20",a2); 
      outcome20[sub] = z20;
      call quad(z30,"fun30",a3); 
      outcome30[sub] = z30;
      call quad(z40,"fun40",a4); 
      outcome40[sub] = z40;
	  call quad(z11,"fun11",a1); 
      outcome11[sub] = z11;
      call quad(z21,"fun21",a2); 
      outcome21[sub] = z21;
      call quad(z31,"fun31",a3); 
      outcome31[sub] = z31;
      call quad(z41,"fun41",a4); 
      outcome41[sub] = z41;

    end;

    out0 = outcome10 + outcome20 + outcome30 + outcome40;
    out1 = outcome11 + outcome21 + outcome31 + outcome41;

    postprobs=subjid||out1||out0||outt1||outt0;
    cname = {"id" "RMST1" "RMST0" "S1" "S0"};

    create &outdata from postprobs  [ colname=cname ];
    append from postprobs; 
  quit;

%mend;

/* Assessed s33 < t* <= smax */

%macro inte415(indata, outdata);

  data temp92;
    set &indata;
  run;

 proc iml;
    use temp60r;
    read all into surv;
    use temp61r;
    read all into baseline;
    use temp13r;
    read all into rcs;

	%if &stat = 1 | &stat = 2 %then %do;
	  %put 1;
	  use temp92;
      read all var{id} into subjid;
      replicate = subjid[,1];
      obs = nrow(replicate);
      read all var{&SURVCov} into cov;

      surv1 = j(obs,1,1)|| cov;
      surv0 = j(obs,1,0)|| cov;
	%end;

    %else %if &stat = 3 | &stat = 4 %then %do;
	  %put 2;
	  use temp92;
      read all var{id} into subjid;
      replicate = subjid[,1];
      obs = nrow(replicate);

      surv1 = j(obs,1,1);
      surv0 = j(obs,1,0);
	%end;

    alpha = surv[1, 1:ncol(surv)];

    gam30 = baseline[1, 1];
    gam31 = baseline[1, 2];
    gam32 = baseline[1, 3];
    gam33 = baseline[1, 4];
    gam34 = baseline[1, 5];

    k_min = rcs[1, 1];
    k_max = rcs[1, 2];
    lamda31 = rcs[1, 6];
    lamda32 = rcs[1, 7];
    lamda33 = rcs[1, 8];

    k31 = rcs[1, 12];
    k32 = rcs[1, 13];
    k33 = rcs[1, 14];

    outt0=j(obs,1,.);
    outt1=j(obs,1,.);
    outcome10=j(obs,1,.);
    outcome20=j(obs,1,.);
    outcome30=j(obs,1,.);
    outcome40=j(obs,1,.);
    outcome50=j(obs,1,.);
    outcome11=j(obs,1,.);
    outcome21=j(obs,1,.);
    outcome31=j(obs,1,.);
    outcome41=j(obs,1,.);
    outcome51=j(obs,1,.);

    start fun10(s) global(alpha, sur0, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx1) + sur0*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun20(s) global(alpha, sur0, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx2) + sur0*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun30(s) global(alpha, sur0, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx3) + sur0*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    start fun40(s) global(alpha, sur0, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx4) + sur0*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    start fun50(s) global(alpha, sur0, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx5) + sur0*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    start fun11(s) global(alpha, sur1, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx1) + sur1*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun21(s) global(alpha, sur1, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx2) + sur1*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun31(s) global(alpha, sur1, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx3) + sur1*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    start fun41(s) global(alpha, sur1, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx4) + sur1*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    start fun51(s) global(alpha, sur1, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx5) + sur1*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

/* Survival Probability (t*) potential outcome (1)*/

    start funn51(s) global(alpha, sur1, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx5) + sur1*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

/* Survival Probability (t*) potential outcome (0)*/

    start funn50(s) global(alpha, sur0, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx5) + sur0*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    do sub = 1 to obs;
      sur1 = surv1[sub,];
      sur0 = surv0[sub,];
      a1   = { 0  &smin }; 
      a2   = { &smin &s31 }; 
      a3   = { &s31 &s32 }; 
      a4   = { &s32 &s33 }; 
      a5   = { &s33  &rmst }; 

  	  outt0[sub] = funn50(&rmst);
	  outt1[sub] = funn51(&rmst);
      call quad(z10,"fun10",a1); 
      outcome10[sub] = z10;
      call quad(z20,"fun20",a2); 
      outcome20[sub] = z20;
      call quad(z30,"fun30",a3); 
      outcome30[sub] = z30;
      call quad(z40,"fun40",a4); 
      outcome40[sub] = z40;
      call quad(z50,"fun50",a5); 
      outcome50[sub] = z50;
	  call quad(z11,"fun11",a1); 
      outcome11[sub] = z11;
      call quad(z21,"fun21",a2); 
      outcome21[sub] = z21;
      call quad(z31,"fun31",a3); 
      outcome31[sub] = z31;
      call quad(z41,"fun41",a4); 
      outcome41[sub] = z41;
      call quad(z51,"fun51",a5); 
      outcome51[sub] = z51;

    end;

    out0 = outcome10 + outcome20 + outcome30 + outcome40 + outcome50;
    out1 = outcome11 + outcome21 + outcome31 + outcome41 + outcome51;

    postprobs=subjid||out1||out0||outt1||outt0;
    cname = {"id" "RMST1" "RMST0" "S1" "S0"};

    create &outdata from postprobs  [ colname=cname ];
    append from postprobs; 
  quit;

%mend;

/* Assessed t* > smax */

%macro inte416(indata, outdata);

  data temp92;
    set &indata;
  run;

 proc iml;
    use temp60r;
    read all into surv;
    use temp61r;
    read all into baseline;
    use temp13r;
    read all into rcs;

	%if &stat = 1 | &stat = 2 %then %do;
	  %put 1;
	  use temp92;
      read all var{id} into subjid;
      replicate = subjid[,1];
      obs = nrow(replicate);
      read all var{&SURVCov} into cov;

      surv1 = j(obs,1,1)|| cov;
      surv0 = j(obs,1,0)|| cov;
	%end;

    %else %if &stat = 3 | &stat = 4 %then %do;
	  %put 2;
	  use temp92;
      read all var{id} into subjid;
      replicate = subjid[,1];
      obs = nrow(replicate);

      surv1 = j(obs,1,1);
      surv0 = j(obs,1,0);
	%end;

    alpha = surv[1, 1:ncol(surv)];

    gam30 = baseline[1, 1];
    gam31 = baseline[1, 2];
    gam32 = baseline[1, 3];
    gam33 = baseline[1, 4];
    gam34 = baseline[1, 5];

    k_min = rcs[1, 1];
    k_max = rcs[1, 2];
    lamda31 = rcs[1, 6];
    lamda32 = rcs[1, 7];
    lamda33 = rcs[1, 8];

    k31 = rcs[1, 12];
    k32 = rcs[1, 13];
    k33 = rcs[1, 14];

    outt0=j(obs,1,.);
    outt1=j(obs,1,.);
    outcome10=j(obs,1,.);
    outcome20=j(obs,1,.);
    outcome30=j(obs,1,.);
    outcome40=j(obs,1,.);
    outcome50=j(obs,1,.);
    outcome60=j(obs,1,.);
    outcome11=j(obs,1,.);
    outcome21=j(obs,1,.);
    outcome31=j(obs,1,.);
    outcome41=j(obs,1,.);
    outcome51=j(obs,1,.);
    outcome61=j(obs,1,.);

    start fun10(s) global(alpha, sur0, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx1) + sur0*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun20(s) global(alpha, sur0, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx2) + sur0*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun30(s) global(alpha, sur0, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx3) + sur0*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    start fun40(s) global(alpha, sur0, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx4) + sur0*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    start fun50(s) global(alpha, sur0, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx5) + sur0*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    start fun60(s) global(alpha, sur0, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx6) + sur0*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    start fun11(s) global(alpha, sur1, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx1) + sur1*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun21(s) global(alpha, sur1, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx2) + sur1*t(alpha)), 600);
      v =  exp((-1)*exp(n));
      return(v); 
    finish; 

    start fun31(s) global(alpha, sur1, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx3) + sur1*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    start fun41(s) global(alpha, sur1, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx4) + sur1*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    start fun51(s) global(alpha, sur1, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx5) + sur1*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    start fun61(s) global(alpha, sur1, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx6) + sur1*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

/* Survival Probability (t*) potential outcome (1)*/

    start funn61(s) global(alpha, sur1, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx6) + sur1*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

/* Survival Probability (t*) potential outcome (0)*/

    start funn60(s) global(alpha, sur0, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32, lamda33, k31, k32, k33); 
      n = min(((&sx6) + sur0*t(alpha)), 600);
      v = exp((-1)*exp(n)) ;
      return(v); 
    finish; 

    do sub = 1 to obs;
      sur1 = surv1[sub,];
      sur0 = surv0[sub,];
      a1   = { 0  &smin }; 
      a2   = { &smin &s31 }; 
      a3   = { &s31 &s32 }; 
      a4   = { &s32 &s33 }; 
      a5   = { &s33  &smax}; 
      a6 = { &smax &rmst };

  	  outt0[sub] = funn60(&rmst);
	  outt1[sub] = funn61(&rmst);
      call quad(z10,"fun10",a1); 
      outcome10[sub] = z10;
      call quad(z20,"fun20",a2); 
      outcome20[sub] = z20;
      call quad(z30,"fun30",a3); 
      outcome30[sub] = z30;
      call quad(z40,"fun40",a4); 
      outcome40[sub] = z40;
      call quad(z50,"fun50",a5); 
      outcome50[sub] = z50;
      call quad(z60,"fun60",a6); 
      outcome60[sub] = z60;
	  call quad(z11,"fun11",a1); 
      outcome11[sub] = z11;
      call quad(z21,"fun21",a2); 
      outcome21[sub] = z21;
      call quad(z31,"fun31",a3); 
      outcome31[sub] = z31;
      call quad(z41,"fun41",a4); 
      outcome41[sub] = z41;
      call quad(z51,"fun51",a5); 
      outcome51[sub] = z51;
      call quad(z61,"fun61",a6); 
      outcome61[sub] = z61;

    end;

    out0 = outcome10 + outcome20 + outcome30 + outcome40 + outcome50 + outcome60;
    out1 = outcome11 + outcome21 + outcome31 + outcome41 + outcome51 + outcome61;

    postprobs=subjid||out1||out0||outt1||outt0;
    cname = {"id" "RMST1" "RMST0" "S1" "S0"};

    create &outdata from postprobs  [ colname=cname ];
    append from postprobs; 
  quit;

%mend;

%macro exp();

  %do i = 1 %to (&nboot + 1);

/*    dm  'log;clear;out;clear;';*/

    data temp13r;
      set temp13;
      if replicate = (&i - 1);
      drop replicate;
    run;

    data _dat32;
      set _dat05;
      if replicate = (&i - 1);
    run;

    data temp46r;
      set temp46;
      if replicate = (&i - 1);
      drop replicate;
    run;

    data temp60r;
      set temp60;
      if replicate = (&i - 1);
      drop replicate;
    run;

    data temp63r;
      set temp63;
      if replicate = (&i - 1);
      drop replicate;
    run;

    data _null_;
      set temp13r;
      call symput('smin', exp(k_min));
      call symput('smax', exp(k_max));
      call symput('s11', exp(k11));
      call symput('s21', exp(k21));
      call symput('s22', exp(k22));
      call symput('s31', exp(k31));
      call symput('s32', exp(k32));
      call symput('s33', exp(k33));
    run;

    data _null_;
      set temp46r;
      call symput('ID', ID);
    run;

	data _null_;
      set temp63r;
      call symput('rmst', max_time);
    run;

	%put &rmst &id &smin &smax &s11 &s21 &s22 &s31 &s32 &s33;

    %if &id ge 1 and &id le 44 %then %do;

      %if &id ge 1 and &id le 8 %then %do;

        data temp61r;
          set temp61;
          if replicate = (&i - 1);
          keep lamd0 lamd1;
        run;

      %end;

      %if &id ge 9 and &id le 44 %then %do;

        data temp61r;
          set temp61;
          if replicate = (&i - 1);
          keep lamd0 lamd1 lamd2;
        run;

      %end;

      data _null_;
        format parameter parameter1 $250.;
        if &id = 1 then do;
          parameter = 'lamd0 + lamd1 * s ** (-2)';
          call symput('sx',parameter);
          parameter1 = '(-2) * lamd1 * s ** (-3)';
          call symput('dsx',parameter1);
        end;
        if &id = 2 then do;
          parameter = 'lamd0 + lamd1 * s ** (-1)';
          call symput('sx',parameter);
          parameter1 = '(-1) * lamd1 * s ** (-2)';
          call symput('dsx',parameter1);
        end;
        if &id = 3 then do;
          parameter = 'lamd0 + lamd1 * s ** (-0.5)';
          call symput('sx',parameter);
          parameter1 = '(-0.5) * lamd1 * s ** (-1.5)';
          call symput('dsx',parameter1);
        end;
        if &id = 4 then do;
          parameter = 'lamd0 + lamd1 * log(s)';
          call symput('sx',parameter);
          parameter1 = 'lamd1 / s';
          call symput('dsx',parameter1);
        end;
        if &id = 5 then do;
          parameter = 'lamd0 + lamd1 * s ** (0.5)';
          call symput('sx',parameter);
          parameter1 = '0.5 * lamd1 * s ** (-0.5)';
          call symput('dsx',parameter1);
        end;
        if &id = 6 then do;
          parameter = 'lamd0 + lamd1 * s ';
          call symput('sx',parameter);
          parameter1 = 'lamd1';
          call symput('dsx',parameter1);
        end;
        if &id = 7 then do;
          parameter = 'lamd0 + lamd1 * s **2';
          call symput('sx',parameter);
          parameter1 = '2 * lamd1 * s';
          call symput('dsx',parameter1);
        end;
        if &id = 8 then do;
          parameter = 'lamd0 + lamd1 * s ** 3';
          call symput('sx',parameter);
          parameter1 = '3 * lamd1 * s ** 2';
          call symput('dsx',parameter1);
        end;

        if &id = 9 then do;
          parameter = 'lamd0 + lamd1 * s ** (-2) + lamd2 * s ** (-2) * log(s)';
          call symput('sx',parameter);
          parameter1 = '(-2) * lamd1 * s ** (-3) + lamd2 * s** (-3) + (-2) * lamd2 * s ** (-3) * log(s)';
          call symput('dsx',parameter1);
        end;
        if &id = 10 then do;
          parameter = 'lamd0 + lamd1 * s ** (-2) + lamd2 * s ** (-1)';
          call symput('sx',parameter);
          parameter1 = '(-2) * lamd1 * s ** (-3) + (-1) * lamd2 * s ** (-2)';
          call symput('dsx',parameter1);
        end;
        if &id = 11 then do;
          parameter = 'lamd0 + lamd1 * s ** (-2) + lamd2 * s ** (-0.5)';
          call symput('sx',parameter);
          parameter1 = '(-2) * lamd1 * s ** (-3) + (-0.5) * lamd2 * s** (-1.5)';
          call symput('dsx',parameter1);
        end;
        if &id = 12 then do;
          parameter = 'lamd0 + lamd1 * s ** (-2) + lamd2 * log(s)';
          call symput('sx',parameter);
          parameter1 = '(-2) * lamd1 * s ** (-3) + lamd2 /s';
          call symput('dsx',parameter1);
        end;
        if &id = 13 then do;
          parameter = 'lamd0 + lamd1 * s ** (-2) + lamd2 * s ** 0.5';
          call symput('sx',parameter);
          parameter1 = '(-2) * lamd1 * s ** (-3) + (0.5) * lamd2 * s** (-0.5)';
          call symput('dsx',parameter1);
        end;
        if &id = 14 then do;
          parameter = 'lamd0 + lamd1 * s ** (-2) + lamd2 * s';
          call symput('sx',parameter);
          parameter1 = '(-2) * lamd1 * s ** (-3) + lamd2';
          call symput('dsx',parameter1);
        end;
        if &id = 15 then do;
          parameter = 'lamd0 + lamd1 * s ** (-2) + lamd2 * s ** 2';
          call symput('sx',parameter);
          parameter1 = '(-2) * lamd1 * s ** (-3) + 2*lamd2 * s';
          call symput('dsx',parameter1);
        end;
        if &id = 16 then do;
          parameter = 'lamd0 + lamd1 * s ** (-2) + lamd2 * s ** 3';
          call symput('sx',parameter);
          parameter1 = '(-2) * lamd1 * s ** (-3) + 3 * lamd2 * s** 2';
          call symput('dsx',parameter1);
        end;

        if &id = 17 then do;
          parameter = 'lamd0 + lamd1 * s ** (-1) + lamd2 * s ** (-1) * log(s)';
          call symput('sx',parameter);
          parameter1 = '(-1) * lamd1 * s ** (-2) + lamd2 * s** (-2) + (-1) * lamd2 * s ** (-2) * log(s)';
          call symput('dsx',parameter1);
        end;
        if &id = 18 then do;
          parameter = 'lamd0 + lamd1 * s ** (-1) + lamd2 * s ** (-0.5)';
          call symput('sx',parameter);
          parameter1 = '(-1) * lamd1 * s ** (-2) + (-0.5) * lamd2 * s ** (-1.5)';
          call symput('dsx',parameter1);
        end;
        if &id = 19 then do;
          parameter = 'lamd0 + lamd1 * s ** (-1) + lamd2 * log(s)';
          call symput('sx',parameter);
          parameter1 = '(-1) * lamd1 * s ** (-2) + lamd2 /s';
          call symput('dsx',parameter1);
        end;
        if &id = 20 then do;
          parameter = 'lamd0 + lamd1 * s ** (-1) + lamd2 * s ** 0.5';
          call symput('sx',parameter);
          parameter1 = '(-1) * lamd1 * s ** (-2) + 0.5 * lamd2 * s ** (-0.5)';
          call symput('dsx',parameter1);
        end;
        if &id = 21 then do;
          parameter = 'lamd0 + lamd1 * s ** (-1) + lamd2 * s';
          call symput('sx',parameter);
          parameter1 = '(-1) * lamd1 * s ** (-2) + lamd2';
          call symput('dsx',parameter1);
        end;
        if &id = 22 then do;
          parameter = 'lamd0 + lamd1 * s ** (-1) + lamd2 * s ** 2';
          call symput('sx',parameter);
          parameter1 = '(-1) * lamd1 * s ** (-2) + 2 * lamd2 * s';
          call symput('dsx',parameter1);
        end;
        if &id = 23 then do;
          parameter = 'lamd0 + lamd1 * s ** (-1) + lamd2 * s ** 3';
          call symput('sx',parameter);
          parameter1 = '(-1) * lamd1 * s ** (-2) + 3 * lamd2 * s** 2';
          call symput('dsx',parameter1);
        end;

        if &id = 24 then do;
          parameter = 'lamd0 + lamd1 * s ** (-0.5) + lamd2 * s ** (-0.5) * log(s)';
          call symput('sx',parameter);
          parameter1 = '(-0.5) * lamd1 * s ** (-1.5) + lamd2 * s** (-1.5) + (-0.5) * lamd2 * s ** (-1.5) * log(s)';
          call symput('dsx',parameter1);
        end;
        if &id = 25 then do;
          parameter = 'lamd0 + lamd1 * s ** (-0.5) + lamd2 * log(s)';
          call symput('sx',parameter);
          parameter1 = '(-0.5) * lamd1 * s ** (-1.5) + lamd2 / s';
          call symput('dsx',parameter1);
        end;
        if &id = 26 then do;
          parameter = 'lamd0 + lamd1 * s ** (-0.5) + lamd2 * s ** 0.5';
          call symput('sx',parameter);
          parameter1 = '(-0.5) * lamd1 * s ** (-1.5) + 0.5 * lamd2 * s ** (-0.5)';
          call symput('dsx',parameter1);
        end;
        if &id = 27 then do;
          parameter = 'lamd0 + lamd1 * s ** (-0.5) + lamd2 * s';
          call symput('sx',parameter);
          parameter1 = '(-0.5) * lamd1 * s ** (-1.5) + lamd2';
          call symput('dsx',parameter1);
        end;
        if &id = 28 then do;
          parameter = 'lamd0 + lamd1 * s ** (-0.5) + lamd2 * s ** 2';
          call symput('sx',parameter);
          parameter1 = '(-0.5) * lamd1 * s ** (-1.5) + 2 * lamd2 * s';
          call symput('dsx',parameter1);
        end;
        if &id = 29 then do;
          parameter = 'lamd0 + lamd1 * s ** (-0.5) + lamd2 * s ** 3';
          call symput('sx',parameter);
          parameter1 = '(-0.5) * lamd1 * s ** (-1.5) + 3 * lamd2 * s ** 2';
          call symput('dsx',parameter1);
        end;

        if &id = 30 then do;
          parameter = 'lamd0 + lamd1 * log(s) + lamd2 * log(s) * log(s)';
          call symput('sx',parameter);
          parameter1 = 'lamd1 / s  + 2 * lamd2 * log(s)/s';
          call symput('dsx',parameter1);
        end;
        if &id = 31 then do;
          parameter = 'lamd0 + lamd1 * log(s) + lamd2 * s ** (0.5)';
          call symput('sx',parameter);
          parameter1 = 'lamd1 / s  + 0.5 * lamd2 * s ** (-0.5)';
          call symput('dsx',parameter1);
        end;
        if &id = 32 then do;
          parameter = 'lamd0 + lamd1 * log(s) + lamd2 * s ';
          call symput('sx',parameter);
          parameter1 = 'lamd1 / s  + lamd2 ';
          call symput('dsx',parameter1);
        end;
        if &id = 33 then do;
          parameter = 'lamd0 + lamd1 * log(s) + lamd2 * s ** 2';
          call symput('sx',parameter);
          parameter1 = 'lamd1 / s  + 2 * lamd2 * s';
          call symput('dsx',parameter1);
        end;
        if &id = 34 then do;
          parameter = 'lamd0 + lamd1 * log(s) + lamd2 * s**3';
          call symput('sx',parameter);
          parameter1 = 'lamd1 / s + lamd2 * s**2 * 3';
          call symput('dsx',parameter1);
        end;

        if &id = 35 then do;
          parameter = 'lamd0 + lamd1 * s ** (0.5) + lamd2 * s ** (0.5) * log(s)';
          call symput('sx',parameter);
          parameter1 = '0.5 * lamd1 * s ** (-0.5) + lamd2 * s** (-0.5) + 0.5 * lamd2 * s ** (-0.5) * log(s)';
          call symput('dsx',parameter1);
        end;
        if &id = 36 then do;
          parameter = 'lamd0 + lamd1 * s ** (0.5) + lamd2 * s';
          call symput('sx',parameter);
          parameter1 = '0.5 * lamd1 * s ** (-0.5) + lamd2';
          call symput('dsx',parameter1);
        end;
        if &id = 37 then do;
          parameter = 'lamd0 + lamd1 * s ** (0.5) + lamd2 * s ** 2';
          call symput('sx',parameter);
          parameter1 = '0.5 * lamd1 * s ** (-0.5) + 2 * lamd2 * s';
          call symput('dsx',parameter1);
        end;
        if &id = 38 then do;
          parameter = 'lamd0 + lamd1 * s ** (0.5) + lamd2 * s ** 3';
          call symput('sx',parameter);
          parameter1 = '0.5 * lamd1 * s ** (-0.5) + 3 * lamd2 * s ** 2';
          call symput('dsx',parameter1);
        end;

        if &id = 39 then do;
          parameter = 'lamd0 + lamd1 * s + lamd2 * s * log(s)';
          call symput('sx',parameter);
          parameter1 = 'lamd1 + lamd2 + lamd2 * log(s)';
          call symput('dsx',parameter1);
        end;
        if &id = 40 then do;
          parameter = 'lamd0 + lamd1 * s + lamd2 * s ** 2';
          call symput('sx',parameter);
          parameter1 = 'lamd1 + 2 * lamd2 * s';
          call symput('dsx',parameter1);
        end;
        if &id = 41 then do;
          parameter = 'lamd0 + lamd1 * s + lamd2 * s ** 3';
          call symput('sx',parameter);
          parameter1 = 'lamd1 + 3 * lamd2 * s**2';
          call symput('dsx',parameter1);
        end;

        if &id = 42 then do;
          parameter = 'lamd0 + lamd1 * s**2 + lamd2 * s**2 * log(s)';
          call symput('sx',parameter);
          parameter1 = '2*lamd1*s + lamd2*s + 2*lamd2 * s* log(s)';
          call symput('dsx',parameter1);
        end;
        if &id = 43 then do;
          parameter = 'lamd0 + lamd1 * s**2 + lamd2 * s**3';
          call symput('sx',parameter);
          parameter1 = '2*lamd1*s + 3*lamd2*s**2';
          call symput('dsx',parameter1);
        end;

        if &id = 44 then do;
          parameter = 'lamd0 + lamd1 * s**3 + lamd2 * s**3 * log(s)';
          call symput('sx',parameter);
          parameter1 = '3*lamd1*s**2 + lamd2*s**2 + 3*lamd2 * s**2* log(s)';
          call symput('dsx',parameter1);
        end;
      run;

      %inte11(_dat32, final1);

    %end;

    %if &id = 45 %then %do;

        data temp61r;
          set temp61;
          if replicate = (&i - 1);
          keep gam10 gam11 gam12;
        run;

        data _null_;
          format parameter1 parameter2 parameter3 parameter4 parameter5 parameter6 parameter7 parameter8 $250.;
          parameter1 = 'gam10 + gam11 * log(s)';
          call symput('sx1',parameter1);
          parameter2 = 'gam11';
          call symput('dsx1',parameter2);
          parameter3 = 'gam10 + gam11 * log(s) + gam12 * (- lamda11 * (log(s) - k_min)**3)';
          call symput('sx2',parameter3);
          parameter4 = 'gam11 + gam12 * (- lamda11 * (log(s) - k_min)**2*3)';
          call symput('dsx2',parameter4);
          parameter5 = 'gam10 + gam11 * log(s) + gam12 * ((log(s) - k11)**3 - lamda11 * (log(s) - k_min)**3)';
          call symput('sx3',parameter5);
          parameter6 = 'gam11 + gam12 * ((log(s) - k11)**2*3 - lamda11 * (log(s) - k_min)**2*3)';
          call symput('dsx3',parameter6);
          parameter7 = 'gam10 + gam11 * log(s) + gam12 * ((log(s) - k11)**3 - lamda11 * (log(s) - k_min)**3 - (1 - lamda11) * (log(s) - k_max)**3)';
          call symput('sx4',parameter7);
          parameter8 = 'gam11 + gam12 * ((log(s) - k11)**2*3 - lamda11 * (log(s) - k_min)**2*3 - (1 - lamda11) * (log(s) - k_max)**2*3)';
          call symput('dsx4',parameter8);
        run;

        %if &rmst <= &smin %then %do;

          %inte211(_dat32, final1);

        %end;

        %else %if &rmst <= &s11 %then %do;

          %inte212(_dat32, final1);

        %end;

		%else %if &rmst <= &smax %then %do;

          %inte213(_dat32, final1);

        %end;

    	%else %do;

          %inte214(_dat32, final1);

        %end;

      %end;

      %if &id = 46 %then %do;

        data temp61r;
          set temp61;
          if replicate = (&i - 1);
          keep gam20 gam21 gam22 gam23;
        run;

        data _null_;
          format parameter1 parameter2 parameter3 parameter4 parameter5 parameter6 parameter7 parameter8 parameter9 parameter10 $500.;
          parameter1 = 'gam20 + gam21 * log(s)';
          call symput('sx1',parameter1);
          parameter2 = 'gam21';
          call symput('dsx1',parameter2);

          parameter3 = 'gam20 + gam21 * log(s) - gam22 * (lamda21 * (log(s) - k_min)**3) - gam23 * (lamda22 * (log(s) - k_min)**3)';
          call symput('sx2',parameter3);
          parameter4 = 'gam21 + gam22 * (- lamda21 * (log(s) - k_min)**2*3) + gam23 * (- lamda22 * (log(s) - k_min)**2*3)';
          call symput('dsx2',parameter4);

          parameter5 = 'gam20 + gam21 * log(s) + gam22 * ((log(s) - k21)**3 - lamda21 * (log(s) - k_min)**3) - gam23 * (lamda22 * (log(s) - k_min)**3)';
          call symput('sx3',parameter5);
          parameter6 = 'gam21 + gam22 * ((log(s) - k21)**2*3 - lamda21 * (log(s) - k_min)**2*3) - gam23 * (3 * lamda22 * (log(s) - k_min)**2)';
          call symput('dsx3',parameter6);

          parameter7 = 'gam20 + gam21 * log(s) + gam22 * ((log(s) - k21)**3 - lamda21 * (log(s) - k_min)**3) + gam23 * ((log(s) - k22)**3 - lamda22 * (log(s) - k_min)**3)';
          call symput('sx4',parameter7);
          parameter8 = 'gam21 + gam22 * ((log(s) - k21)**2*3 - lamda21 * (log(s) - k_min)**2*3) + gam23 * (3 * (log(s) - k22)**2 - 3 * lamda22 * (log(s) - k_min)**2)';
          call symput('dsx4',parameter8);

          parameter9 = 'gam20 + gam21 * log(s) + gam22 * ((log(s) - k21)**3 - lamda21 * (log(s) - k_min)**3 - (1 - lamda21) * (log(s) - k_max)**3)+gam23 * ((log(s) - k22)**3 - lamda22 * (log(s) - k_min)**3 - (1 - lamda22) * (log(s) - k_max)**3)';
          call symput('sx5',parameter9);
          parameter10 = 'gam21 + gam22 * ((log(s) - k21)**2*3 - lamda21 * (log(s) - k_min)**2*3 - (1 - lamda21) * (log(s) - k_max)**2*3)+gam23 * ((log(s) - k22)**2*3 - lamda22 * (log(s) - k_min)**2*3 - (1 - lamda22) * (log(s) - k_max)**2*3)';
          call symput('dsx5',parameter10);

        run;

        %if &rmst <= &smin %then %do;

          %inte311(_dat32, final1);

        %end;

        %else %if &rmst <= &s21 %then %do;

          %inte312(_dat32, final1);

        %end;

		%else %if &rmst <= &s22 %then %do;

          %inte313(_dat32, final1);

        %end;

		%else %if &rmst <= &smax %then %do;

          %inte314(_dat32, final1);

        %end;

		%else %do;

          %inte315(_dat32, final1);

        %end;

      %end;

      %if &id = 47 %then %do;

        data temp61r;
          set temp61;
          if replicate = (&i - 1);
          keep gam30 gam31 gam32 gam33 gam34;
        run;

        data _null_;
          format parameter1 parameter2 parameter3 parameter4 parameter5 parameter6 parameter7 parameter8 parameter9 parameter10 parameter11 $500.;
          parameter1 = 'gam30 + gam31 * log(s)';
          call symput('sx1',parameter1);
          parameter2 = 'gam31';
          call symput('dsx1',parameter2);

          parameter3 = 'gam30 + gam31 * log(s) + gam32 * (- lamda31 * (log(s) - k_min)**3) + gam33 * (-lamda32 * (log(s) - k_min)**3) + gam34 * (-lamda33 * (log(s) - k_min)**3)';
          call symput('sx2',parameter3);
          parameter4 = 'gam31 + gam32 * (- lamda31 * (log(s) - k_min)**2*3) + gam33 * (- lamda32 * (log(s) - k_min)**2*3)+ gam34 * (- lamda33 * (log(s) - k_min)**2*3)';
          call symput('dsx2',parameter4);

          parameter5 = 'gam30 + gam31 * log(s) + gam32 * ((log(s) - k31)**3 - lamda31 * (log(s) - k_min)**3) + gam33 * (-lamda32 * (log(s) - k_min)**3)+ gam34 * (-lamda33 * (log(s) - k_min)**3)';
          call symput('sx3',parameter5);
          parameter6 = 'gam31 + gam32 * ((log(s) - k31)**2*3 - lamda31 * (log(s) - k_min)**2*3) + gam33 * (-lamda32 * (log(s) - k_min)**2*3) + gam34 * (-lamda33 * (log(s) - k_min)**2*3)';
          call symput('dsx3',parameter6);

          parameter7 = 'gam30 + gam31 * log(s) + gam32 * ((log(s) - k31)**3 - lamda31 * (log(s) - k_min)**3) + gam33 * ((log(s) - k32)**3 - lamda32 * (log(s) - k_min)**3) + gam34 * (-lamda33 * (log(s) - k_min)**3)';
          call symput('sx4',parameter7);
          parameter8 = 'gam31 + gam32 * ((log(s) - k31)**2*3 - lamda31 * (log(s) - k_min)**2*3) + gam33 * (3 * (log(s) - k32)**2 - 3 * lamda32 * (log(s) - k_min)**2) + gam34 * (-lamda33 * (log(s) - k_min)**2*3)';
          call symput('dsx4',parameter8);

          parameter9 = 'gam30 + gam31 * log(s) + gam32 * ((log(s) - k31)**3-lamda31*(log(s)-k_min)**3)+gam33 *((log(s)-k32)**3-lamda32 *(log(s)-k_min)**3)+gam34 *((log(s) - k33)**3-lamda33*(log(s)-k_min)**3)';
          call symput('sx5',parameter9);
          parameter10 = 'gam31+gam32*((log(s)-k31)**2*3-lamda31*(log(s)-k_min)**2*3)+gam33*(3 * (log(s) - k32)**2 - 3 * lamda32 * (log(s) - k_min)**2) + gam34 * (3 * (log(s) - k33)**2 -lamda33 * (log(s) - k_min)**2*3)';
          call symput('dsx5',parameter10);

		  parameter11 = 'gam30+gam31*log(s)+gam32*((log(s)-k31)**3-lamda31*(log(s)-k_min)**3-(1-lamda31)*(log(s)-k_max)**3)+gam33*((log(s)-k32)**3-lamda32*(log(s)-k_min)**3-(1-lamda32)*(log(s)-k_max)**3)+gam34*((log(s)-k33)**3-lamda33*(log(s)-k_min)**3-(1-lamda33)*(log(s)-k_max)**3)';
          call symput('sx6',parameter11);
          *parameter12 = 'gam31+gam32*((log(s)-k31)**2*3-lamda31*(log(s)-k_min)**2*3-(1-lamda31)*(log(s)-k_max)**2*3)+gam33*((log(s)-k32)**2*3-lamda32*(log(s)-k_min)**2*3-(1-lamda32)*(log(s)-k_max)**2*3)+gam34*((log(s)-k33)**2*3-lamda33*(log(s)-k_min)**2*3-(1-lamda33)*(log(s)-k_max)**2*3)';
          *call symput('dsx6',parameter12);
        run;

		%put &sx6;

        %if &rmst <= &smin %then %do;

          %inte411(_dat32, final1);

        %end;

        %else %if &rmst <= &s31 %then %do;

          %inte412(_dat32, final1);

        %end;

		%else %if &rmst <= &s32 %then %do;

          %inte413(_dat32, final1);

        %end;

		%else %if &rmst <= &s33 %then %do;

          %inte414(_dat32, final1);

        %end;

		%else %if &rmst <= &smax %then %do;

          %inte415(_dat32, final1);

        %end;

		%else %do;

          %inte416(_dat32, final1);

        %end;

      %end;

      data final;
        set final final1 (in = a);
        if a then replicate = &i - 1;
      run;

      proc sort data = final;
        by replicate id;
      run;

      proc datasets lib = work;
        delete final1 temp61r temp63r temp60r temp46r _dat32 temp13r;
      run;

  %end;

%mend;

*******************************************************************************************************
                                             Main Macro;
*******************************************************************************************************;

%macro CUREEXP   (data             = , 
                  subjid           = ,
                  CUREcov          = ,
                  SURVcov          = ,
                  EXP              = , 
                  nboot            = 200, 
                  SURVLtime        = ,
                  SURVRtime        = ,
                  SURVevent        = ,
				  T_ACCESS         = 999999,
				  REFERENCE        = 1,
                  seedboot         = 1237893,
                  Printall         = T,
                  Out             = 
                  );

data zzz0;
  indata = "&data";
  subjid = "&subjid";
  survLtime = "&survLtime";
  survRtime = "&survRtime";
  survevent = "&survevent";
  curecov = "&curecov";
  survcov = "&survcov";
  exp = "&exp";
  reference = &reference;
  if indata = '' or SURVLtime=' ' or SURVRtime=' ' or SURVevent=' ' or EXP = ' ' 
  or compress(indata, '') = "''" or compress(SURVLtime, '') = "''" or compress(SURVRtime, '') = "''" or compress(SURVevent, '') = "''" or compress(EXP, '') = "''"  
  then sta = 0; 
  else if (curecov = '' or compress(curecov, '') = "''") and (survcov = '' or compress(survcov, '') = "''") then sta = 4; 
  else if (curecov = '' or compress(curecov, '') = "''") then sta = 2; 
  else if (survcov = '' or compress(survcov, '') = "''") then sta = 3; 
  else sta = 1;
run;

data _null_;
  set zzz0;
  call symput('stat', sta);
run;

/* 
Make sure input data set name, outcome (including survival outcome left boundary, survival outcome right boundary, and event indicator) 
and expsoure variables are not missing 
stat = 0: miss data set name, exposure or outcome variables
stat = 1: curvcov and survcov exist
stat = 2: curecov does not exist and survcov does exist
stat = 3: curecov does exist and survcov does not exist
stat = 4: curecov does not exist and survcov does not exist
*/

%if &stat = 0 %then %do;
  proc iml;
    print,"WARNING: NEED TO SPECIFY DATASET NAME, EXPOSURE AND OUTCOME VARIABLE","PROGRAM WILL TERMINATE";
  quit;
%end;

%else %do;

ods listing close;

*******************************************************************************************************
Determine Parametric Survival Model;
*******************************************************************************************************;

/* Simple Data Manipulation, to select knot of RCS, using non-censored observation middle point of left and right censoring */

data aa11;
  set &data;
run;

proc surveyselect data= aa11 out=boots
  method = urs
  samprate = 1 outhits rep = &nboot seed = &seedboot;
run;

data aa12;
  set aa11 (in = a) boots ;
  if a then replicate = 0;
  drop NumberHits;
run;

data aa1;
  set aa12;
  by replicate &subjid;
  retain subjid;
  if first.replicate then subjid = 0;
  subjid = subjid + 1;
run;

data aa31;
  set aa1;
  if &SURVRtime = . then time = &SURVLtime;
  else time = (&SURVRtime + &SURVLtime)/2;
run;

proc means data = aa31;
  var time;
  output out = aa42 (drop = _type_ _freq_) max(time) = max_time1;
  by replicate;
  where &SURVevent = 1;
run;

proc means data = aa31;
  var time;
  output out = aaa42 (drop = _type_ _freq_) max(&SURVLtime) = max_time2 min(&SURVLtime) = min_time2;
  by replicate;
  where &SURVevent = 0;
run;

data aa3;
  merge aa31 aa42 aaa42;
  by replicate;
  if max(max_time1, max_time2) = max_time1 then max_time = max_time2 - (max_time2 - min_time2)/5;
  if max(max_time1, max_time2) = max_time2 then max_time = max_time1;
  if time le max_time then sus = 1;
  else if time ne . then sus = 0;
  drop max_time max_time1 max_time2 min_time2;
run;

proc freq data = aa3;
  tables &SURVevent * &SURVRtime sus &SURVevent * sus/list missing;
  by replicate;
run;

/* Isolate all distinct survival time */

data temp1 (keep = replicate &SURVLtime rename = (&SURVLtime = SURVtime)) temp2 (keep = replicate &SURVRtime rename = (&SURVRtime = SURVtime));
  set aa3;
run;

data temp3;
  set temp1 temp2;
  event = 1;
  if SURVtime ne .;
run;

proc sort data = temp3 nodupkey;
  by replicate SURVtime;
run;

/* Logistic Initial Value */

proc logistic data=aa3 ;
  model sus = &exp &curecov; 
  ods output parameterestimates=beta
  ConvergenceStatus = con11;
  by replicate;
run;

data tt61;
  set beta;
  if replicate = 0;
  i = 1;
  np2 = _n_;
run;

data tt62;
  set tt61;
  by i;
  if last.i;
  keep np2;
run; 

data _null_;
  set tt62;
  call symput('NP2', np2);
run;

** create dataset of initial Cure binary model parameters for NLMIXED;
data beta1; length parameter $8; set beta END=lastobs;
  by replicate;
  format row 3.0 linpred $255.;
  row = (_n_ - replicate * &np2)-1;
  parameter = "a" || left(row);
  *create binary cure probability linear predictor for NLMIXED;
  if not first.replicate then linpred =  "+" || cats(parameter) || "*" || variable;
  else linpred =  cats(parameter);
run;

data beta2;
  set beta1;
  if replicate = 0;
  if _N_=1 then call symput('linpredCURE',linpred);
  else call symput('linpredCURE',trim(resolve('&linpredCURE'))||linpred);
run; 

data final6;
set beta2;
run;

%put &linpredCURE &np2;

data beta3;
  set beta1;
  format model $4.;
  model = 'Cure';
  rename variable = effect;
  keep replicate model variable parameter estimate;
run;

proc datasets lib = work;
  delete tt61 tt62 beta beta beta1 beta2;
run;

/* RCS */

data _dat; 
  set aa3; 
  _time = time;
  _event = &SURVevent;
  lnt = log(_time);
  if sus = 1;
run;

data _dat1;
  set _dat;
  if _event = 1;
run;

proc univariate data = _dat1 noprint;
  var lnt;
  output out=temp11 pctlpts=50 33 67 25 75 pctlpre = k pctlname = _50 _33 _67 _25 _75;
  by replicate;
run;

proc means data = _dat1 noprint;
  var lnt;
  output out = temp12 min(lnt) = k_min max(lnt) = k_max;
  by replicate;
run;

data temp13;
  merge temp11 temp12;
  by replicate;
  drop _type_ _freq_;
  lamda11 = (k_max - k_50)/(k_max - k_min);
  lamda21 = (k_max - k_33)/(k_max - k_min);
  lamda22 = (k_max - k_67)/(k_max - k_min);
  lamda31 = (k_max - k_25)/(k_max - k_min);
  lamda32 = (k_max - k_50)/(k_max - k_min);
  lamda33 = (k_max - k_75)/(k_max - k_min);
  k11 = k_50;
  k21 = k_33;
  k22 = k_67;
  k31 = k_25;
  k32 = k_50;
  k33 = k_75;
  drop k_50 k_33 k_67 k_25 k_75;
run;

data _dat11;
  merge _dat temp13;
  by replicate;

  if lnt le k_min then do;
    nu11 = 0;
    dnu11 = 0;
  end;
  else if lnt le k11 then do;
    nu11 = - lamda11 * (lnt - k_min)**3;
    dnu11 = - lamda11 * (lnt - k_min)**2*3;
  end;
  else if lnt le k_max then do;
    nu11 = (lnt - k11)**3 - lamda11 * (lnt - k_min)**3; 
    dnu11 = (lnt - k11)**2*3 - lamda11 * (lnt - k_min)**2*3; 
  end;
  else do;
    nu11 = (lnt - k11)**3 - lamda11 * (lnt - k_min)**3 - (1 - lamda11) * (lnt - k_max)**3; 
    dnu11 = (lnt - k11)**2*3 - lamda11 * (lnt - k_min)**2*3 - (1 - lamda11) * (lnt - k_max)**2*3; 
  end;

  if lnt le k_min then do;
    nu21 = 0;
    dnu21 = 0;
    nu22 = 0;
    dnu22 = 0;
  end;
  else if lnt le k21 then do;
    nu21 = - lamda21 * (lnt - k_min)**3;
    nu22 = - lamda22 * (lnt - k_min)**3;
    dnu21 = - lamda21 * (lnt - k_min)**2*3;
    dnu22 = - lamda22 * (lnt - k_min)**2*3;
  end;
  else if lnt le k22 then do;
    nu21 = (lnt - k21)**3 - lamda21 * (lnt - k_min)**3;
    nu22 =  - lamda22 * (lnt - k_min)**3;
    dnu21 = (lnt - k21)**2*3 - lamda21 * (lnt - k_min)**2*3;
    dnu22 =  - lamda22 * (lnt - k_min)**2*3;
  end;
  else if lnt le k_max then do;
    nu21 = (lnt - k21)**3 - lamda21 * (lnt - k_min)**3;
    nu22 = (lnt - k22)**3 - lamda22 * (lnt - k_min)**3;
    dnu21 = (lnt - k21)**2*3 - lamda21 * (lnt - k_min)**2*3;
    dnu22 = (lnt - k22)**2*3 - lamda22 * (lnt - k_min)**2*3;
  end;
  else do;
    nu21 = (lnt - k21)**3 - lamda21 * (lnt - k_min)**3 - (1 - lamda21) * (lnt - k_max)**3; 
    dnu21 = (lnt - k21)**2*3 - lamda21 * (lnt - k_min)**2*3 - (1 - lamda21) * (lnt - k_max)**2*3; 
    nu22 = (lnt - k22)**3 - lamda22 * (lnt - k_min)**3 - (1 - lamda22) * (lnt - k_max)**3; 
    dnu22 = (lnt - k22)**2*3 - lamda22 * (lnt - k_min)**2*3 - (1 - lamda22) * (lnt - k_max)**2*3; 
  end;

  if lnt le k_min then do;
    nu31 = 0;
    dnu31 = 0;
    nu32 = 0;
    dnu32 = 0;
    nu33 = 0;
    dnu33 = 0;
  end;
  else if lnt le k31 then do;
    nu31 = - lamda31 * (lnt - k_min)**3;
    nu32 = - lamda32 * (lnt - k_min)**3;
    nu33 = - lamda33 * (lnt - k_min)**3;
    dnu31 = - lamda31 * (lnt - k_min)**2*3;
    dnu32 = - lamda32 * (lnt - k_min)**2*3;
    dnu33 = - lamda33 * (lnt - k_min)**2*3;
  end;
  else if lnt le k32 then do;
    nu31 = (lnt - k31)**3 - lamda31 * (lnt - k_min)**3;
    nu32 = - lamda32 * (lnt - k_min)**3;
    nu33 = - lamda33 * (lnt - k_min)**3;
    dnu31 = (lnt - k31)**2*3 - lamda31 * (lnt - k_min)**2*3;
    dnu32 = - lamda32 * (lnt - k_min)**2*3;
    dnu33 = - lamda33 * (lnt - k_min)**2*3;
  end;
  else if lnt le k33 then do;
    nu31 = (lnt - k31)**3 - lamda31 * (lnt - k_min)**3;
    nu32 = (lnt - k32)**3 - lamda32 * (lnt - k_min)**3;
    nu33 = - lamda33 * (lnt - k_min)**3;
    dnu31 = (lnt - k31)**2*3 - lamda31 * (lnt - k_min)**2*3;
    dnu32 = (lnt - k32)**2*3 - lamda32 * (lnt - k_min)**2*3;
    dnu33 = - lamda33 * (lnt - k_min)**2*3;
  end;
  else if lnt le k_max then do;
    nu31 = (lnt - k31)**3 - lamda31 * (lnt - k_min)**3;
    nu32 = (lnt - k32)**3 - lamda32 * (lnt - k_min)**3;
    nu33 = (lnt - k33)**3 - lamda33 * (lnt - k_min)**3;
    dnu31 = (lnt - k31)**2*3 - lamda31 * (lnt - k_min)**2*3;
    dnu32 = (lnt - k32)**2*3 - lamda32 * (lnt - k_min)**2*3;
    dnu33 = (lnt - k33)**2*3 - lamda33 * (lnt - k_min)**2*3;
  end;
  else do;
    nu31 = (lnt - k31)**3 - lamda31 * (lnt - k_min)**3 - (1 - lamda31) * (lnt - k_max)**3; 
    dnu31 = (lnt - k31)**2*3 - lamda31 * (lnt - k_min)**2*3 - (1 - lamda31) * (lnt - k_max)**2*3; 
    nu32 = (lnt - k32)**3 - lamda32 * (lnt - k_min)**3 - (1 - lamda32) * (lnt - k_max)**3; 
    dnu32 = (lnt - k32)**2*3 - lamda32 * (lnt - k_min)**2*3 - (1 - lamda32) * (lnt - k_max)**2*3; 
    nu33 = (lnt - k33)**3 - lamda33 * (lnt - k_min)**3 - (1 - lamda33) * (lnt - k_max)**3; 
    dnu33 = (lnt - k33)**2*3 - lamda33 * (lnt - k_min)**2*3 - (1 - lamda33) * (lnt - k_max)**2*3; 
  end;
run;

data _dat21;
  set _dat11;
  if _event = 1;
  drop &survcov &curecov time &survLtime &survRtime &exp _event;
run;

proc phreg data=_dat;
  class &exp (ref="0" param=ref);
  model _time*_event(0)= &exp &SURVCov;
  ods output ParameterEstimates=SURVests1(rename=(Parameter=Effect));
  baseline out = bb5 survival = _all_ CUMHAZ=_ALL_;
  by replicate;
run;

data tt51;
  set survests1;
  if replicate = 0;
  i = 1;
  np1 = _n_;
run;

data tt52;
  set tt51;
  by i;
  if last.i;
  keep np1;
run; 

data _null_;
  set tt52;
  call symput('NP1', np1);
run;

*Survival model parameters;

data SURVests2; 
  length model $4 parameter $8; 
  set SURVests1;
  by replicate;
  format row 3.0;
  row = (_n_ - replicate * &np1)-1;
  parameter = "b" || left(row);
  if DF=0 then delete; *remove unestimated parameters;
  model="Surv";
  if not first.replicate then linpred =  "+" || cats(parameter) || "*" || effect;
  else linpred =  cats(parameter) || "*" || effect;
run; 

data final5;
set SURVests2;
run;

data SURVests3; 
  set SURVests2;
  if replicate = 0;
  if _N_=1 then call symput('linpredSURV',linpred);
  else call symput('linpredSURV',trim(resolve('&linpredSURV'))||linpred);
run;

%put &linpredSURV &np1;

proc datasets lib = work;
  delete temp11 temp12 SURVests3 SURVests1 tt51 tt52;
run;

proc sort data = _dat21;
  by replicate _time;
run;

data bb51;
  set bb5;
  if _time ne 0;
  keep replicate &exp &SURVcov _time CumHaz;
run;

proc sort data = bb51;
  by replicate _time;
run;

proc transpose data = SURVests2 out = survests4;
  var estimate;
  id parameter;
  by replicate;
run;

data bb53;
  merge survests4 bb51;
  by replicate;
  logbaseCumHaz = log(CumHaz) - (&linpredSURV);
  keep replicate _time logbaseCumHaz;
run;

data _dat22;
  merge _dat21 bb53;
  by replicate _time;
run;

proc glm data = _dat22;
  model logbasecumhaz = lnt nu11;
  by replicate;
  ods output ParameterEstimates  = tt11;
run;

proc transpose data = tt11 out = tt12 (drop = _name_ rename = (col1 = gam10 col2 = gam11 col3 = gam12));
  by replicate;
  var estimate;
run;

proc glm data = _dat22;
  model logbasecumhaz = lnt nu21 nu22;
  by replicate;
  ods output ParameterEstimates  = tt21;
run;

proc transpose data = tt21 out = tt22 (drop = _name_ rename = (col1 = gam20 col2 = gam21 col3 = gam22 col4 = gam23));
  by replicate;
  var estimate;
run;

proc glm data = _dat22;
  model logbasecumhaz = lnt nu31 nu32 nu33;
  by replicate;
  ods output ParameterEstimates  = tt31;
run;

proc transpose data = tt31 out = tt32 (drop = _name_ rename = (col1 = gam30 col2 = gam31 col3 = gam32 col4 = gam33 col5 = gam34));
  by replicate;
  var estimate;
run;

data bb54;
  merge tt12 tt22 tt32;
  by replicate;
run;

proc transpose data = bb54 out = bb55 (rename = (_name_ = parameter col1 = estimate));
  var gam10 -- gam34;
  by replicate;
run;

data SurvInit;
  set SURVests2 (keep = replicate model parameter estimate effect) bb55 (in = a);
  if a then model = 'RCS';
run;

data SurvInit0;
  set SurvInit;
  if parameter in ('gam10','gam20', 'gam30') then effect = 'Intercept';
  if parameter in ('gam11','gam21', 'gam31') then effect = 'lnt';
  if parameter in ('gam12') then effect = 'nu11';
  if parameter in ('gam22') then effect = 'nu21';
  if parameter in ('gam23') then effect = 'nu22';
  if parameter in ('gam32') then effect = 'nu31';
  if parameter in ('gam33') then effect = 'nu32';
  if parameter in ('gam34') then effect = 'nu33';
run;

data Init1;
  set beta3 SurvInit0;
  if substr(parameter, 1, 4) in ('gam2', 'gam3') then delete;
run; 

data Init2;
  set beta3 SurvInit0;
  if substr(parameter, 1, 4) in ('gam1', 'gam3') then delete;
run; 

data Init3;
  set beta3 SurvInit0;
  if substr(parameter, 1, 4) in ('gam2', 'gam1') then delete;
run; 

proc sort data = Init1;
  by replicate parameter;
run;

proc sort data = Init2;
  by replicate parameter;
run;

proc sort data = Init3;
  by replicate parameter;
run;

proc datasets lib = work;
  delete bb5 bb51 bb54 bb55 tt11 tt12 tt21 tt22 tt31 tt32 survinit0;
run;

data aa51;
  merge aa3 temp13;
  by replicate;
run;

data aa52;
  set aa51;
  lnt = log(&SURVLtime);
  if &SURVevent = 1 then lnRt = log(&SURVRtime);
run;

data aa53;
  set aa52;
  if lnt le k_min then do;
    nu11 = 0;
    dnu11 = 0;
  end;
  else if lnt le k11 then do;
    nu11 = - lamda11 * (lnt - k_min)**3;
    dnu11 = - lamda11 * (lnt - k_min)**2*3;
  end;
  else if lnt le k_max then do;
    nu11 = (lnt - k11)**3 - lamda11 * (lnt - k_min)**3; 
    dnu11 = (lnt - k11)**2*3 - lamda11 * (lnt - k_min)**2*3; 
  end;
  else do;
    nu11 = (lnt - k11)**3 - lamda11 * (lnt - k_min)**3 - (1 - lamda11) * (lnt - k_max)**3; 
    dnu11 = (lnt - k11)**2*3 - lamda11 * (lnt - k_min)**2*3 - (1 - lamda11) * (lnt - k_max)**2*3; 
  end;

  if lnt le k_min then do;
    nu21 = 0;
    dnu21 = 0;
    nu22 = 0;
    dnu22 = 0;
  end;
  else if lnt le k21 then do;
    nu21 = - lamda21 * (lnt - k_min)**3;
    nu22 = - lamda22 * (lnt - k_min)**3;
    dnu21 = - lamda21 * (lnt - k_min)**2*3;
    dnu22 = - lamda22 * (lnt - k_min)**2*3;
  end;
  else if lnt le k22 then do;
    nu21 = (lnt - k21)**3 - lamda21 * (lnt - k_min)**3;
    nu22 =  - lamda22 * (lnt - k_min)**3;
    dnu21 = (lnt - k21)**2*3 - lamda21 * (lnt - k_min)**2*3;
    dnu22 =  - lamda22 * (lnt - k_min)**2*3;
  end;
  else if lnt le k_max then do;
    nu21 = (lnt - k21)**3 - lamda21 * (lnt - k_min)**3;
    nu22 = (lnt - k22)**3 - lamda22 * (lnt - k_min)**3;
    dnu21 = (lnt - k21)**2*3 - lamda21 * (lnt - k_min)**2*3;
    dnu22 = (lnt - k22)**2*3 - lamda22 * (lnt - k_min)**2*3;
  end;
  else do;
    nu21 = (lnt - k21)**3 - lamda21 * (lnt - k_min)**3 - (1 - lamda21) * (lnt - k_max)**3; 
    dnu21 = (lnt - k21)**2*3 - lamda21 * (lnt - k_min)**2*3 - (1 - lamda21) * (lnt - k_max)**2*3; 
    nu22 = (lnt - k22)**3 - lamda22 * (lnt - k_min)**3 - (1 - lamda22) * (lnt - k_max)**3; 
    dnu22 = (lnt - k22)**2*3 - lamda22 * (lnt - k_min)**2*3 - (1 - lamda22) * (lnt - k_max)**2*3; 
  end;

  if lnt le k_min then do;
    nu31 = 0;
    dnu31 = 0;
    nu32 = 0;
    dnu32 = 0;
    nu33 = 0;
    dnu33 = 0;
  end;
  else if lnt le k31 then do;
    nu31 = - lamda31 * (lnt - k_min)**3;
    nu32 = - lamda32 * (lnt - k_min)**3;
    nu33 = - lamda33 * (lnt - k_min)**3;
    dnu31 = - lamda31 * (lnt - k_min)**2*3;
    dnu32 = - lamda32 * (lnt - k_min)**2*3;
    dnu33 = - lamda33 * (lnt - k_min)**2*3;
  end;
  else if lnt le k32 then do;
    nu31 = (lnt - k31)**3 - lamda31 * (lnt - k_min)**3;
    nu32 = - lamda32 * (lnt - k_min)**3;
    nu33 = - lamda33 * (lnt - k_min)**3;
    dnu31 = (lnt - k31)**2*3 - lamda31 * (lnt - k_min)**2*3;
    dnu32 = - lamda32 * (lnt - k_min)**2*3;
    dnu33 = - lamda33 * (lnt - k_min)**2*3;
  end;
  else if lnt le k33 then do;
    nu31 = (lnt - k31)**3 - lamda31 * (lnt - k_min)**3;
    nu32 = (lnt - k32)**3 - lamda32 * (lnt - k_min)**3;
    nu33 = - lamda33 * (lnt - k_min)**3;
    dnu31 = (lnt - k31)**2*3 - lamda31 * (lnt - k_min)**2*3;
    dnu32 = (lnt - k32)**2*3 - lamda32 * (lnt - k_min)**2*3;
    dnu33 = - lamda33 * (lnt - k_min)**2*3;
  end;
  else if lnt le k_max then do;
    nu31 = (lnt - k31)**3 - lamda31 * (lnt - k_min)**3;
    nu32 = (lnt - k32)**3 - lamda32 * (lnt - k_min)**3;
    nu33 = (lnt - k33)**3 - lamda33 * (lnt - k_min)**3;
    dnu31 = (lnt - k31)**2*3 - lamda31 * (lnt - k_min)**2*3;
    dnu32 = (lnt - k32)**2*3 - lamda32 * (lnt - k_min)**2*3;
    dnu33 = (lnt - k33)**2*3 - lamda33 * (lnt - k_min)**2*3;
  end;
  else do;
    nu31 = (lnt - k31)**3 - lamda31 * (lnt - k_min)**3 - (1 - lamda31) * (lnt - k_max)**3; 
    dnu31 = (lnt - k31)**2*3 - lamda31 * (lnt - k_min)**2*3 - (1 - lamda31) * (lnt - k_max)**2*3; 
    nu32 = (lnt - k32)**3 - lamda32 * (lnt - k_min)**3 - (1 - lamda32) * (lnt - k_max)**3; 
    dnu32 = (lnt - k32)**2*3 - lamda32 * (lnt - k_min)**2*3 - (1 - lamda32) * (lnt - k_max)**2*3; 
    nu33 = (lnt - k33)**3 - lamda33 * (lnt - k_min)**3 - (1 - lamda33) * (lnt - k_max)**3; 
    dnu33 = (lnt - k33)**2*3 - lamda33 * (lnt - k_min)**2*3 - (1 - lamda33) * (lnt - k_max)**2*3; 
  end;

  if &SURVevent = 1 then do;
    if lnRt le k_min then do;
      nuR11 = 0;
      dnuR11 = 0;
    end;
    else if lnRt le k11 then do;
      nuR11 = - lamda11 * (lnRt - k_min)**3;
      dnuR11 = - lamda11 * (lnRt - k_min)**2*3;
    end;
    else if lnRt le k_max then do;
      nuR11 = (lnRt - k11)**3 - lamda11 * (lnRt - k_min)**3; 
      dnuR11 = (lnRt - k11)**2*3 - lamda11 * (lnRt - k_min)**2*3; 
    end;
    else do;
      nuR11 = (lnRt - k11)**3 - lamda11 * (lnRt - k_min)**3 - (1 - lamda11) * (lnRt - k_max)**3; 
      dnuR11 = (lnRt - k11)**2*3 - lamda11 * (lnRt - k_min)**2*3 - (1 - lamda11) * (lnRt - k_max)**2*3; 
    end;

    if lnRt le k_min then do;
      nuR21 = 0;
      dnuR21 = 0;
      nuR22 = 0;
      dnuR22 = 0;
    end;
    else if lnRt le k21 then do;
      nuR21 = - lamda21 * (lnRt - k_min)**3;
      nuR22 = - lamda22 * (lnRt - k_min)**3;
      dnuR21 = - lamda21 * (lnRt - k_min)**2*3;
      dnuR22 = - lamda22 * (lnRt - k_min)**2*3;
    end;
    else if lnRt le k22 then do;
      nuR21 = (lnRt - k21)**3 - lamda21 * (lnRt - k_min)**3;
      nuR22 =  - lamda22 * (lnRt - k_min)**3;
      dnuR21 = (lnRt - k21)**2*3 - lamda21 * (lnRt - k_min)**2*3;
      dnuR22 =  - lamda22 * (lnRt - k_min)**2*3;
    end;
    else if lnRt le k_max then do;
      nuR21 = (lnRt - k21)**3 - lamda21 * (lnRt - k_min)**3;
      nuR22 = (lnRt - k22)**3 - lamda22 * (lnRt - k_min)**3;
      dnuR21 = (lnRt - k21)**2*3 - lamda21 * (lnRt - k_min)**2*3;
      dnuR22 = (lnRt - k22)**2*3 - lamda22 * (lnRt - k_min)**2*3;
    end;
    else do;
      nuR21 = (lnRt - k21)**3 - lamda21 * (lnRt - k_min)**3 - (1 - lamda21) * (lnRt - k_max)**3; 
      dnuR21 = (lnRt - k21)**2*3 - lamda21 * (lnRt - k_min)**2*3 - (1 - lamda21) * (lnRt - k_max)**2*3; 
      nuR22 = (lnRt - k22)**3 - lamda22 * (lnRt - k_min)**3 - (1 - lamda22) * (lnRt - k_max)**3; 
      dnuR22 = (lnRt - k22)**2*3 - lamda22 * (lnRt - k_min)**2*3 - (1 - lamda22) * (lnRt - k_max)**2*3; 
    end;

    if lnRt le k_min then do;
      nuR31 = 0;
      dnuR31 = 0;
      nuR32 = 0;
      dnuR32 = 0;
      nuR33 = 0;
      dnuR33 = 0;
    end;
    else if lnRt le k31 then do;
      nuR31 = - lamda31 * (lnRt - k_min)**3;
      nuR32 = - lamda32 * (lnRt - k_min)**3;
      nuR33 = - lamda33 * (lnRt - k_min)**3;
      dnuR31 = - lamda31 * (lnRt - k_min)**2*3;
      dnuR32 = - lamda32 * (lnRt - k_min)**2*3;
      dnuR33 = - lamda33 * (lnRt - k_min)**2*3;
    end;
    else if lnRt le k32 then do;
      nuR31 = (lnRt - k31)**3 - lamda31 * (lnRt - k_min)**3;
      nuR32 = - lamda32 * (lnRt - k_min)**3;
      nuR33 = - lamda33 * (lnRt - k_min)**3;
      dnuR31 = (lnRt - k31)**2*3 - lamda31 * (lnRt - k_min)**2*3;
      dnuR32 = - lamda32 * (lnRt - k_min)**2*3;
      dnuR33 = - lamda33 * (lnRt - k_min)**2*3;
    end;
    else if lnRt le k33 then do;
      nuR31 = (lnRt - k31)**3 - lamda31 * (lnRt - k_min)**3;
      nuR32 = (lnRt - k32)**3 - lamda32 * (lnRt - k_min)**3;
      nuR33 = - lamda33 * (lnRt - k_min)**3;
      dnuR31 = (lnRt - k31)**2*3 - lamda31 * (lnRt - k_min)**2*3;
      dnuR32 = (lnRt - k32)**2*3 - lamda32 * (lnRt - k_min)**2*3;
      dnuR33 = - lamda33 * (lnRt - k_min)**2*3;
    end;
    else if lnRt le k_max then do;
      nuR31 = (lnRt - k31)**3 - lamda31 * (lnRt - k_min)**3;
      nuR32 = (lnRt - k32)**3 - lamda32 * (lnRt - k_min)**3;
      nuR33 = (lnRt - k33)**3 - lamda33 * (lnRt - k_min)**3;
      dnuR31 = (lnRt - k31)**2*3 - lamda31 * (lnRt - k_min)**2*3;
      dnuR32 = (lnRt - k32)**2*3 - lamda32 * (lnRt - k_min)**2*3;
      dnuR33 = (lnRt - k33)**2*3 - lamda33 * (lnRt - k_min)**2*3;
    end;
    else do;
      nuR31 = (lnRt - k31)**3 - lamda31 * (lnRt - k_min)**3 - (1 - lamda31) * (lnRt - k_max)**3; 
      dnuR31 = (lnRt - k31)**2*3 - lamda31 * (lnRt - k_min)**2*3 - (1 - lamda31) * (lnRt - k_max)**2*3; 
      nuR32 = (lnRt - k32)**3 - lamda32 * (lnRt - k_min)**3 - (1 - lamda32) * (lnRt - k_max)**3; 
      dnuR32 = (lnRt - k32)**2*3 - lamda32 * (lnRt - k_min)**2*3 - (1 - lamda32) * (lnRt - k_max)**2*3; 
      nuR33 = (lnRt - k33)**3 - lamda33 * (lnRt - k_min)**3 - (1 - lamda33) * (lnRt - k_max)**3; 
      dnuR33 = (lnRt - k33)**2*3 - lamda33 * (lnRt - k_min)**2*3 - (1 - lamda33) * (lnRt - k_max)**2*3; 
    end;
  end;
run;

%rcs;

/* Fractional Polynomials */

data _dat31;
  set aa3;
  drop time sus;
  t1 = &SURVLtime ** (-2);
  t2 = &SURVLtime ** (-1);
  t3 = &SURVLtime ** (-0.5);
  t4 = log(&SURVLtime);
  t5 = &SURVLtime ** (0.5);
  t6 = &SURVLtime;
  t7 = &SURVLtime ** 2;
  t8 = &SURVLtime ** 3;

  t9 = &SURVLtime ** (-2) * log(&SURVLtime);
  t10 = &SURVLtime ** (-1) * log(&SURVLtime);
  t11 = &SURVLtime ** (-0.5) * log(&SURVLtime);
  t12 = log(&SURVLtime) * log(&SURVLtime);
  t13 = &SURVLtime ** (0.5) * log(&SURVLtime);
  t14 = &SURVLtime * log(&SURVLtime);
  t15 = &SURVLtime ** 2 * log(&SURVLtime);
  t16 = &SURVLtime ** 3 * log(&SURVLtime);

  if &SURVevent = 1 then do;
    tR1 = &SURVRtime ** (-2);
    tR2 = &SURVRtime ** (-1);
    tR3 = &SURVRtime ** (-0.5);
    tR4 = log(&SURVRtime);
    tR5 = &SURVRtime ** (0.5);
    tR6 = &SURVRtime;
    tR7 = &SURVRtime ** 2;
    tR8 = &SURVRtime ** 3;

    tR9 = &SURVRtime ** (-2) * log(&SURVRtime);
    tR10 = &SURVRtime ** (-1) * log(&SURVRtime);
    tR11 = &SURVRtime ** (-0.5) * log(&SURVRtime);
    tR12 = log(&SURVRtime) * log(&SURVRtime);
    tR13 = &SURVRtime ** (0.5) * log(&SURVRtime);
    tR14 = &SURVRtime * log(&SURVRtime);
    tR15 = &SURVRtime ** 2 * log(&SURVRtime);
    tR16 = &SURVRtime ** 3 * log(&SURVRtime);
  end;
run;

data bb56;
  set bb53;
  t1 = _time ** (-2);
  t2 = _time ** (-1);
  t3 = _time ** (-0.5);
  t4 = log(_time);
  t5 = _time ** (0.5);
  t6 = _time;
  t7 = _time ** 2;
  t8 = _time ** 3;

  t9 = _time ** (-2) * log(_time);
  t10 = _time ** (-1) * log(_time);
  t11 = _time ** (-0.5) * log(_time);
  t12 = log(_time) * log(_time);
  t13 = _time ** (0.5) * log(_time);
  t14 = _time * log(_time);
  t15 = _time ** 2 * log(_time);
  t16 = _time ** 3 * log(_time);
run;

%single(_dat31, bb56, ff01, hh01, gg01, 1, 1);
%single(_dat31, bb56, ff02, hh02, gg02, 2, 2);
%single(_dat31, bb56, ff03, hh03, gg03, 3, 3);
%single(_dat31, bb56, ff04, hh04, gg04, 4, 4);
%single(_dat31, bb56, ff05, hh05, gg05, 5, 5);
%single(_dat31, bb56, ff06, hh06, gg06, 6, 6);
%single(_dat31, bb56, ff07, hh07, gg07, 7, 7);
%single(_dat31, bb56, ff08, hh08, gg08, 8, 8);

%double(_dat31, bb56, ff09, hh09, gg09, 1, 9, 9);
%double(_dat31, bb56, ff10, hh10, gg10, 1, 2, 10);
%double(_dat31, bb56, ff11, hh11, gg11, 1, 3, 11);
%double(_dat31, bb56, ff12, hh12, gg12, 1, 4, 12);
%double(_dat31, bb56, ff13, hh13, gg13, 1, 5, 13);
%double(_dat31, bb56, ff14, hh14, gg14, 1, 6, 14);
%double(_dat31, bb56, ff15, hh15, gg15, 1, 7, 15);
%double(_dat31, bb56, ff16, hh16, gg16, 1, 8, 16);

%double(_dat31, bb56, ff17, hh17, gg17, 2, 10, 17);
%double(_dat31, bb56, ff18, hh18, gg18, 2, 3, 18);
%double(_dat31, bb56, ff19, hh19, gg19, 2, 4, 19);
%double(_dat31, bb56, ff20, hh20, gg20, 2, 5, 20);
%double(_dat31, bb56, ff21, hh21, gg21, 2, 6, 21);
%double(_dat31, bb56, ff22, hh22, gg22, 2, 7, 22);
%double(_dat31, bb56, ff23, hh23, gg23, 2, 8, 23);

%double(_dat31, bb56, ff24, hh24, gg24, 3, 11, 24);
%double(_dat31, bb56, ff25, hh25, gg25, 3, 4, 25);
%double(_dat31, bb56, ff26, hh26, gg26, 3, 5, 26);
%double(_dat31, bb56, ff27, hh27, gg27, 3, 6, 27);
%double(_dat31, bb56, ff28, hh28, gg28, 3, 7, 28);
%double(_dat31, bb56, ff29, hh29, gg29, 3, 8, 29);

%double(_dat31, bb56, ff30, hh30, gg30, 4, 12, 30);
%double(_dat31, bb56, ff31, hh31, gg31, 4, 5, 31);
%double(_dat31, bb56, ff32, hh32, gg32, 4, 6, 32);
%double(_dat31, bb56, ff33, hh33, gg33, 4, 7, 33);
%double(_dat31, bb56, ff34, hh34, gg34, 4, 8, 34);

%double(_dat31, bb56, ff35, hh35, gg35, 5, 13, 35);
%double(_dat31, bb56, ff36, hh36, gg36, 5, 6, 36);
%double(_dat31, bb56, ff37, hh37, gg37, 5, 7, 37);
%double(_dat31, bb56, ff38, hh38, gg38, 5, 8, 38);

%double(_dat31, bb56, ff39, hh39, gg39, 6, 14, 39);
%double(_dat31, bb56, ff40, hh40, gg40, 6, 7, 40);
%double(_dat31, bb56, ff41, hh41, gg41, 6, 8, 41);

%double(_dat31, bb56, ff42, hh42, gg42, 7, 15, 42);
%double(_dat31, bb56, ff43, hh43, gg43, 7, 8, 43);

%double(_dat31, bb56, ff44, hh44, gg44, 8, 16, 44);

/* Selection among 47 FP and RCS models */

data temp40;
  format method $3.0;
  set gg01 gg02 gg03 gg04 gg05 gg06 gg07 gg08 gg09 gg10
  gg11 gg12 gg13 gg14 gg15 gg16 gg17 gg18 gg19 gg20
  gg21 gg22 gg23 gg24 gg25 gg26 gg27 gg28 gg29 gg30
  gg31 gg32 gg33 gg34 gg35 gg36 gg37 gg38 gg39 gg40
  gg41 gg42 gg43 gg44 gg45 gg46 gg47;
run;

data temp50;
  format method $3.0;
  set ff01 ff02 ff03 ff04 ff05 ff06 ff07 ff08 ff09 ff10
  ff11 ff12 ff13 ff14 ff15 ff16 ff17 ff18 ff19 ff20
  ff21 ff22 ff23 ff24 ff25 ff26 ff27 ff28 ff29 ff30
  ff31 ff32 ff33 ff34 ff35 ff36 ff37 ff38 ff39 ff40
  ff41 ff42 ff43 ff44 ff45 ff46 ff47;
run;

data temp40;
  set temp40;
  where converge = 'Yes' and AIC ne .;
run;

proc sort data = temp50;
  by replicate id;
run;

proc sort data = temp40;
  by replicate aic;
run;

data temp46;
  set temp40;
  by replicate aic;
  if first.replicate;
run;

data temp43;
  merge temp46 (in = a keep = replicate id) temp50;
  by replicate id;
  if a;
run;

proc transpose data = temp43 out = temp60 (drop = _name_ id);
  by replicate id;
  id parameter;
  var estimate;
  where substr(parameter, 1, 1) in ('b');
run;

proc transpose data = temp43 out = temp64 (drop = _name_ id);
  by replicate id;
  id parameter;
  var estimate;
  where substr(parameter, 1, 1) in ('a');
run;

proc transpose data = temp43 out = temp61 (drop = _name_ id);
  by replicate id;
  id parameter;
  var estimate;
  where substr(parameter, 1, 1) in ('g', 'l');
run;

proc means data = temp3 noprint;
  by replicate;
  var SURVtime;
  output out = temp66 (drop = _type_ _freq_) max(SURVtime) = max_time;
run;

data final;
  replicate = -1;
run;

%if &reference = 0 %then %do;

  proc sort data = aa53 out = _dat00;
    by replicate &SURVCov &CURECov;
  run;

  proc means data = _dat00 noprint;
    var &exp;
    where &exp = 0;
    output out = _dat01(drop = _type_ _freq_) n(&exp) = nsub;
    by replicate &SURVCov &CURECov;
  run;

  proc means data = _dat00 noprint;
    var &exp;
    where &exp = 0;
    output out = _dat03(drop = _type_ _freq_) n(&exp) = total;
    by replicate;
  run;

%end;

%if &reference = 1 %then %do;

  proc sort data = aa53 out = _dat00;
    by replicate &SURVCov &CURECov;
  run;

  proc means data = _dat00 noprint;
    var &exp;
    where &exp = 1;
    output out = _dat01(drop = _type_ _freq_) n(&exp) = nsub;
    by replicate &SURVCov &CURECov;
  run;

  proc means data = _dat00 noprint;
    var &exp;
    where &exp = 1;
    output out = _dat03(drop = _type_ _freq_) n(&exp) = total;
    by replicate;
  run;

%end;

%if &reference = 2 %then %do;

  proc sort data = aa53 out = _dat00;
    by replicate &SURVCov &CURECov;
  run;

  proc means data = _dat00 noprint;
    var &exp;
    output out = _dat01(drop = _type_ _freq_) n(&exp) = nsub;
    by replicate &SURVCov &CURECov;
  run;

  proc means data = _dat00 noprint;
    var &exp;
    output out = _dat03(drop = _type_ _freq_) n(&exp) = total;
    by replicate;
  run;

%end;

%if &reference = 3 %then %do;

  proc sort data = ref out = _dat00;
    by &SURVCov &CURECov;
  run;

  data _dat01;
    set _dat00;
    do i = 1 to (&nboot + 1);
    replicate = i - 1;
    output;
    end;
    drop i;
  run;

  proc sort data = _dat01;
    by replicate &SURVCov &CURECov;
  run;

  proc means data = _dat01 noprint;
    var nsub;
    output out = _dat03(drop = _type_ _freq_) sum(nsub) = total;
    by replicate;
  run;

%end;

data _dat05;
  merge _dat01 _dat03;
  by replicate;
  percent = nsub/total * 100;
  drop nsub total;
  retain id;
  if first.replicate then id = 0;
  id = id + 1;
run;

%if &T_ACCESS ^= 999999 %then %do;
  data temp63;
	set temp66;
	max_time = &T_ACCESS;
  run;
%end;

%else %do;
  data temp63;
	set temp66;
  run;
%end;

%put &t_access;

%exp();

data final0;
  set final;
  if replicate ne (-1);
run;

data final11;
  merge final0 temp64 temp63;
  by replicate;
run;

data final12;
  merge final11 _dat05;
  by replicate id;
run;

data final131;
  set final12;
  &exp = 1;
  pi1 = exp(min(&linpredCURE, 600))/(1+exp(min(&linpredCURE, 600)));
  mu1 = (1-pi1)*max_time + pi1 * rmst1;
  muu1 = (1-pi1) + pi1 * s1;
run;

data final132;
  set final12;
  &exp = 0;
  pi0 = exp(min(&linpredCURE, 600))/(1+exp(min(&linpredCURE, 600)));
  mu0 = (1-pi0)*max_time + pi0 * rmst0;
  muu0 = (1-pi0) + pi0 * s0;
run;

data final13;
  merge final131 final132;
  by replicate id;
  exp = mu1 - mu0;
  expp = muu1 - muu0;
run;

data final14;
  set final13;
  exp = exp * percent/100;
  expp = expp * percent/100;
run;

proc means data = final14 noprint;
  var exp;
  by replicate;
  output out = final15 (drop = _type_ _freq_) sum(exp) = exp sum(expp) = expp;
run;

data ff0;
  set final15;
  if replicate = 0;
  drop replicate;
run;

data ff1;
  set final15;
  if replicate ne 0;
run;

proc univariate data = ff1 noprint;
  var exp expp;
  output out=ff2 pctlpts=2.5 97.5  pctlpre = exp expp pctlname = _025 _975;
run;

data ff3;
  merge ff0 ff2;
run;

proc sort data = temp40;
  by replicate id;
run;

data out21;
  merge temp46 (in = a keep = replicate id) temp40;
  by replicate id;
  if a;
  drop p1 p2 nknot method;
run;

data out22;
  set temp43;
  drop p1 p2 nknot method;
run;

data out23;
  set final15;
run;

data out24;
  set ff3;
run;

/*proc datasets lib = work;*/
/*  save _dat05 temp60 temp64 temp43 temp3 temp63 aa1 out21 out22 out23 out24 temp66;*/
/*run;*/

data aa3;
  set aa1;
  if replicate = 0;
run;

data mm00;
  set temp3;
  by replicate;
  if replicate = 0;
  if last.replicate;
  ind = _n_;
run;

data _null_;
  set mm00;
  call symput('MMAX', IND);
run;

%put &mmax;

%if &mmax gt 300 %then %do;

  data mm01;
    set temp3;
    by replicate;
    if replicate = 0;
    if first.replicate or last.replicate;
  run;

  proc transpose data = mm01 out = mm02 (drop = _name_ rename = (col1 = min_surv col2 = max_surv));
    by replicate;
    var SURVtime;
  run;

  data mm03;
    merge temp3 mm02;
    by replicate;
    if replicate = 0;
    do i = 1 to 300;
      if survtime ge (min_surv + (i - 1) * (max_surv - min_surv)/300) and survtime lt (min_surv + i * (max_surv - min_surv)/300) 
      then survtime =  min_surv + (i - 1) * (max_surv - min_surv)/300;
    end;
    drop i max_surv min_surv;
  run;

  proc sort data = mm03 out = mm05 nodupkey;
    by survtime;
  run;

  data temp3;
    set mm05;
  run;

  data mm04;
    merge aa3 mm02;
    by replicate;
    if replicate = 0;
    do i = 1 to 300;
      if &SURVLtime ge (min_surv + (i - 1) * (max_surv - min_surv)/300) and &SURVLtime lt (min_surv + i * (max_surv - min_surv)/300) 
      then &SURVLtime =  min_surv + (i - 1) * (max_surv - min_surv)/300;
      if &SURVRtime ge (min_surv + (i - 1) * (max_surv - min_surv)/300) and &SURVRtime lt (min_surv + i * (max_surv - min_surv)/300) 
      then &SURVRtime =  min_surv + (i - 1) * (max_surv - min_surv)/300;
    end;
    drop i max_surv min_surv;
  run;

  data aa3;
    set mm04;
  run;

  proc datasets lib = work;
    drop mm01 mm02 mm03 mm04 mm05;
  run;

%end;

/* Initial value for Gamma */

data temp31;
  set temp3;
  by replicate;
  retain ind;
  if first.replicate then ind = 0;
  ind = ind + 1;
  &SURVLtime = survtime;
  &SURVRtime = survtime;
  if replicate = 0;
  keep replicate ind &SURVLtime &SURVRtime;
run;

proc sort data = aa3 out = aa32;
  by replicate &SURVLtime;
run;

data aa33;
  merge aa32(in = a) temp31 (drop = &SURVRtime);
  by replicate &SURVLtime;
  if a;
  rename ind = indL;
run;

proc sort data = aa33 out = aa34;
  by replicate &SURVRtime;
run;

data aa35;
  merge aa34(in = a) temp31 (drop = &SURVLtime);
  by replicate &SURVRtime;
  if a;
  rename ind = indR;
run;

proc sort data = aa35;
  by replicate subjid;
run;

/* Add indL and indR as survival time stepwise indicators */

proc phreg data = temp3;
  model SURVtime * event (0) = ;
  baseline out = temp91 survival = _all_ CUMHAZ=_ALL_;
  by replicate;
run;

data temp92;
  set temp91;
  by replicate descending survival;
  retain ind;
  if first.replicate then ind = 0;
  ind = ind + 1;
  rename survival = survivalL;
  keep ind survival replicate;
run;

data temp93;
  set temp91;
  by replicate descending survival;
  retain ind;
  if first.replicate then ind = -1;
  ind = ind + 1;
  rename survival = survivalR;
  keep ind survival replicate;
run;

data temp94;
  merge temp92 temp93;
  by replicate ind;
  if survivalL = . or survivalR = . then delete;
run;

data temp95;
  set temp94;
  ga = log(log(survivalL) - log(survivalR));
  if ga gt 0 then ga = 0;
run;

data out31 out32 out33 out34 out35;
  replicate = -1;
run;

%do i = 1 %to 1;

  data temp95r;
    set temp95;
    if replicate = (&i - 1);
  run;

  data temp64r;
    set temp64;
    if replicate = (&i - 1);
  run;

  data temp60r;
    set temp60;
    if replicate = (&i - 1);
  run;

  data temp43r;
    set temp43;
    if replicate = (&i - 1);
  run;

  data aa35r;
    set aa35;
    if replicate = (&i - 1);
  run;

  data tt71;
    set temp95r;
    np3 = _n_;
  run;

  data tt72;
    set tt71;
    by replicate;
    if last.replicate;
    keep np3;
  run; 

  data _null_;
    set tt72;
    call symput('NP3', np3);
  run;

  *Survival baseline hazard parameters;

  data base1; 
    length model $4 parameter $5; 
    set temp95r;
    model = 'Base';
    by replicate;
    format row 3.0;
    row = _n_;
    parameter = "g" || left(row);
    keep model parameter replicate ga;
    rename ga = estimate;
  run;

  data base11;
    set base1;
    by replicate estimate;
    if last.replicate;
  run;

  data base12;
    set base1;
    by replicate estimate;
    if first.replicate;
  run;

  data _null_;
    set base11;
    call symput('gammaR', parameter);
  run;

  data _null_;
    set base12;
    call symput('gammaL', parameter);
  run;

  data base2;
    retain replicate model parameter estimate;
    set base1;
  run;

  proc transpose data = base2 out = base3 (drop = _name_);
    id  parameter;
    var estimate;
    by replicate;
  run;

  data init;
    set temp43r base2;
    drop p1 p2 id nknot StandardError;
    format model $4.;
    if substr(parameter, 1, 1) = 'a' then model = 'Cure';
    if substr(parameter, 1, 1) = 'b' then model = 'Surv';
    if substr(parameter, 1, 1) not in ('a', 'b', 'g') then delete;
    if substr(parameter, 1, 2)in ('ga') then delete;
  run;

  proc nlmixed data=aa35r cov maxfunc = 10000 maxiter = 1000;
    parms /bydata data=init;
	bounds &gammaL - &gammaR <= 0;

	cureprop = exp(min(&linpredCURE, 600))/(1+exp(min(&linpredCURE, 600)));

    array x[&np3] &gammaL - &gammaR;
    sumL = 0;
    do i = 1 to indL;
      sumL = sumL + exp(min((x[i] + &linpredSURV), 600));
    end;
	if &SURVevent = 1 then do;
      sumR = 0;
      do j = 1 to indR;
        sumR = sumR + exp(min((x[j] + &linpredSURV), 600));
      end;
      loglike = log(cureprop) + log(exp((-1) * sumL) - exp((-1) * sumR)); 
	end;
	if &SURVevent = 0 then loglike = log((1 - cureprop) + cureprop * exp((-1) * sumL)); 

	model &SURVevent ~ general(loglike);
 
    ods output ParameterEstimates=temp25;
    ods output COVMatParmEst = temp26; 
    ods output FitStatistics = temp29; 
    ods output ConvergenceStatus = temp27 (keep = replicate status rename = (status = value));
	by replicate;
  run;

  data temp30;
    set temp29 temp27 (in = a);
    if a then Descr = 'Con';
  run;

  proc sort data = temp30;
    by replicate;
  run;

  proc transpose data = temp30 out = temp32 (rename = (AIC__smaller_is_better_ = AIC AICC__smaller_is_better_ = AICC BIC__smaller_is_better_ = BIC));
    id descr;
    var value;
    by replicate;
  run;

  data temp28;
    ite = 1;
  run;

  data aa42 aa43 aa44;
    set temp25;
    if substr(parameter, 1, 1) = 'a' then output aa42;
    if substr(parameter, 1, 1) = 'b' then output aa43;
    if substr(parameter, 1, 1) = 'g' then output aa44;
  run;

  data temp261;
    set temp26;
    drop replicate row parameter;
  run;

  data _dat05r;
    set _dat05;
    if replicate = (&i - 1);
  run;

  data temp3r;
    set temp3;
    if replicate = (&i - 1);
  run;

  data temp63r;
    set temp63;
    if replicate = (&i - 1);
  run;

  data temp66r;
    set temp66;
    if replicate = (&i - 1);
	t_access = &t_access;
	if t_access = 999999 then teind = 0;
	else if t_access le max_time then teind = 1;
	else if t_access gt max_time then teind = 2;
  	if t_access ne 999999 then max_time = t_access;
	keep t_access max_time teind;
  run;

  data _null_;
    set temp66r;
    call symput('rmst', max_time);
    call symput('teind', teind);
  run;

  %put &rmst &teind;

  /* If assessment is less than last observed time */

  %if &teind = 1 %then %do;

    %put 1;

	data temp31r;
	  SURVtime = &rmst;
	  replicate = 0;
	  event = 1;
	run;

	data temp32r;
	  set temp3r temp31r;
	run;

	proc sort data = temp32r;
	  by replicate survtime;
	run;

	data temp33r;
	  merge temp32r temp31r (keep = replicate survtime rename = (survtime = survtime1));
	  by replicate;
	  if survtime gt survtime1 then delete;
	  drop survtime1;
	run;

	proc sort data = temp33r out = temp34r nodupkey;
	  by replicate survtime;
	run;

	data temp35r;
	  set temp34r;
	  by replicate;
	  if last.replicate;
	  ind = _n_;
	  keep replicate ind;
	run;

	data temp36r;
	  merge aa44 temp35r;
	  by replicate;
	  if _n_ gt ind then delete;
	  drop ind;
	run;

	proc iml;
	  use temp261;
	  read all into cov;
      use aa43;
      read all var{estimate} into survest;
      use aa42;
      read all var{estimate} into cureest;
      use temp36r;
      read all var{replicate} into baseest;
      ncure = nrow(cureest);
      nsurv = nrow(survest);
      nbase = nrow(baseest);

	  cov1 = cov[1:(ncure + nsurv + nbase), 1:(ncure + nsurv + nbase)];
      create temp37r from cov1;
      append from cov1;
	  quit;
	run;

	data temp261;
	  set temp37r;
	run;

	data aa44;
	  set temp36r;
	run;

	data temp3r;
	  set temp34r;
	run;

   	proc datasets lib = work;
	  delete temp31r temp32r temp33r temp34r temp35r temp36r temp37r;
	run;

  %end;

  /* If assessment is equal or less than last observed time */

  %if &teind = 1 |  &teind = 0 %then %do;

	%put 0;
    proc iml;
      use temp261;
      read all into cov;
      use aa43;
      read all var{estimate} into survest;
      use aa42;
      read all var{estimate} into cureest;
      use aa44;
      read all var{estimate} into baseest;
      use temp3r;
      read all var{SURVtime} into int;
      ntime = nrow(int);
      max_time = &rmst;

	  %if &stat = 1 %then %do;
	    %put 1;
	    use _dat05r;
        read all var{id} into id;
        read all var{percent} into per;
        read all var{replicate} into replicate;
	    nobs = replicate[,1];
        obs = nrow(nobs);
        read all var{&SURVCov} into SURVcov;
        read all var{&CURECov} into CUREcov;

        cure1 = j(obs,1,1)|| j(obs,1,1)|| CUREcov;
        cure0 = j(obs,1,1)|| j(obs,1,0)|| CUREcov;
        surv1 = j(obs,1,1)|| SURVcov;
        surv0 = j(obs,1,0)|| SURVcov;
      %end;

      %else %if &stat = 2 %then %do;
	    %put 2;
	    use _dat05r;
        read all var{id} into id;
        read all var{percent} into per;
        read all var{replicate} into replicate;
	    nobs = replicate[,1];
        obs = nrow(nobs);
        read all var{&SURVCov} into SURVcov;

        cure1 = j(obs,1,1)|| j(obs,1,1);
        cure0 = j(obs,1,1)|| j(obs,1,0);
        surv1 = j(obs,1,1)|| SURVcov;
        surv0 = j(obs,1,0)|| SURVcov;
      %end;

      %else %if &stat = 3 %then %do;
	    %put 3;
	    use _dat05r;
        read all var{id} into id;
        read all var{percent} into per;
        read all var{replicate} into replicate;
	    nobs = replicate[,1];
        obs = nrow(nobs);
        read all var{&CURECov} into CUREcov;

        cure1 = j(obs,1,1)|| j(obs,1,1)|| CUREcov;
        cure0 = j(obs,1,1)|| j(obs,1,0)|| CUREcov;
        surv1 = j(obs,1,1);
        surv0 = j(obs,1,0);
      %end;

      %else %if &stat = 4 %then %do;
	    %put 4;
	    use _dat05r;
        read all var{id} into id;
        read all var{percent} into per;
        read all var{replicate} into replicate;
	    nobs = replicate[,1];
        obs = nrow(nobs);

        cure1 = j(obs,1,1)|| j(obs,1,1);
        cure0 = j(obs,1,1)|| j(obs,1,0);
        surv1 = j(obs,1,1);
        surv0 = j(obs,1,0);
      %end;

	  /* remove missing covariance value in covariance matrix */
	  cov1 = cov[loc(cov[1,] ^=.),];
      cov2 = cov1[,loc(cov[1,] ^=.)];

      int1 = 0//int[1:(nrow(int)-1)];
      interv = int - int1;

      min_b = j(obs, 1, 600);
      temp1 = j(obs, ntime, .);
      tem1 = j(obs, ntime, .);
      te1 = j(obs, ntime, .);
      do i = 1 to ntime;
        temp1[,i] = exp(min_b >< (j(obs, 1, baseest[i]) + surv1 * survest));
    	te1[,i] = exp(min_b >< (j(obs, 1, baseest[i]) + surv1 * survest));
	  end;
      do p = 2 to ntime;
	    temp1[,p] = temp1[,p] + temp1[,p-1];
	  end;
	  /* Shift one column for Restricted Mean Survival Time Calculation */
	  tem1 = j(obs, 1, 0)||temp1[, 1:(ntime-1)];
	  s1 = j(obs, ntime, .);
	  s1 = exp((-1) * tem1);
      /* No need to shift for survival probability at Tmax Calculation at step down */
	  ss1 = j(obs, 1, .);
      /* No need to shift for survival probability at Tmax Calculation at step up */
      ss1 = (exp((-1) * tem1))[, ntime];
      tee1 = te1[, 1:(ntime-1)]||j(obs, 1, 0);

      temp0 = j(obs, ntime, .);
      tem0 = j(obs, ntime, .);
      te0 = j(obs, ntime, .);
      do m = 1 to ntime;
      	temp0[,m] = exp(min_b >< (j(obs, 1, baseest[m]) + surv0 * survest));
	    te0[,m] = exp(min_b >< (j(obs, 1, baseest[m]) + surv0 * survest));
	  end;
      do n = 2 to ntime;
	    temp0[,n] = temp0[,n] + temp0[,n-1];
	  end;
      /* Shift one column for Restricted Mean Survival Time Calculation */
	  tem0 = j(obs, 1, 0)||temp0[, 1:(ntime-1)];
	  s0 = j(obs, ntime, .);
	  s0 = exp((-1) * tem0);
      /* No need to shift for survival probability at Tmax Calculation */
	  ss0 = j(obs, 1, .);
      /* No need to shift for survival probability at Tmax Calculation at step up */
      ss0 = (exp((-1) * tem0))[, ntime];
      tee0 = te0[, 1:(ntime-1)]||j(obs, 1, 0);

      rmst1 = j(obs, 1, .);
      rmst0 = j(obs, 1, .);
      rmst1 = s1 * interv;
      rmst0 = s0 * interv;

      cureprop1 = exp(min_b >< (cure1 * cureest)) /(1 + exp(min_b >< (cure1 * cureest)));
      cureprop0 = exp(min_b >< (cure0 * cureest)) /(1 + exp(min_b >< (cure0 * cureest)));

	  /* Overall Exposure Effect on RMST */
      exp = j(obs, 1, .);
      exp = (((1 - cureprop1) * max_time + cureprop1 # rmst1) - ((1 - cureprop0) * max_time + cureprop0 # rmst0))#per/100;
      apv = sum(exp); 
	  /* Overall Exposure Effect on survival probability at T* */
      expp = j(obs, 1, .);
      expp = (((1 - cureprop1) + cureprop1 # ss1) - ((1 - cureprop0) + cureprop0 # ss0))#per/100;
      apvv = sum(expp); 
	  print apv apvv;

      ncure = nrow(cureest);
      nsurv = nrow(survest);
      nbase = nrow(baseest);

      k = nrow(survest) + nrow(cureest) + nrow(baseest);
      /* Delta method for variance calculation for exposure effect on RMST */
      g = j(1,k,0);
      /* Delta method for variance calculation for exposure effect on probability at Tmax */
      gg = j(1,k,0);
      lincure1 = cure1 * cureest;
      lincure0 = cure0 * cureest;
      linsurv1 = surv1 * survest;
      linsurv0 = surv0 * survest;

      do i = 1 to obs; 
        t11 = exp(min(600, lincure1[i]))/(1+exp(min(600, lincure1[i])))**2 * cure1[i,];
        t10 = exp(min(600, lincure0[i]))/(1+exp(min(600, lincure0[i])))**2 * cure0[i,];
        g[1:ncure] = g[1:ncure] + t((((-1) * max_time * t11 + rmst1[i] * t11) - ((-1) * max_time * t10 + rmst0[i] * t10))*per[i]/100);
        gg[1:ncure] = gg[1:ncure] + t((((-1) * t11 + ss1[i] * t11) - ((-1) * t10 + ss0[i] * t10))*per[i]/100);

        t21 = (s1 # tem1 * (-1) * interv);
        t20 = (s0 # tem0 * (-1) * interv);
        g[(ncure + 1):(ncure + nsurv)] = g[(ncure + 1):(ncure + nsurv)] + t((cureprop1[i] * t21[i] * surv1[i,] - cureprop0[i] * t20[i] * surv0[i,])*per[i]/100);
        tt21 = (ss1 # tem1[, ntime] * (-1));
        tt20 = (ss0 # tem0[, ntime] * (-1));
        gg[(ncure + 1):(ncure + nsurv)] = gg[(ncure + 1):(ncure + nsurv)] + t((cureprop1[i] * tt21[i] * surv1[i,] - cureprop0[i] * tt20[i] * surv0[i,])*per[i]/100);

        /* Generate Triangle Matrix */
        temp11 = j(ntime, ntime, 1);
        r = row(temp11);      /* create helper matrices */
        c = col(temp11);
        lowerIdx = loc(r >= c);
        temp11[lowerIdx] = 0;

        temp12 = t(t(s1[i,]) * j(1, ntime, 1));
        temp13 = t(te1[i,]) * j(1, ntime, 1) * (-1);
        t31 = temp11 # temp12 # temp13 * interv;

        temp02 = t(t(s0[i,]) * j(1, ntime, 1));
        temp03 = t(te0[i,]) * j(1, ntime, 1) * (-1);
        t30 = temp11 # temp02 # temp03 * interv;
        g[(ncure + nsurv + 1):k] = g[(ncure + nsurv + 1):k] + ((cureprop1[i] * t31 - cureprop0[i] * t30)*per[i]/100);

        tt31 = t(ss1[i] * tee1[i, ] * (-1));
        tt30 = t(ss0[i] * tee0[i, ] * (-1));
        gg[(ncure + nsurv + 1):k] = gg[(ncure + nsurv + 1):k] + ((cureprop1[i] * tt31 - cureprop0[i] * tt30)*per[i]/100);
      end;
      g1 = t(g[loc(cov[1,] ^=.)]);
      gg1 = t(gg[loc(cov[1,] ^=.)]);
      finalvar = g1 * cov2 * t(g1);
      finalsd = sqrt(finalvar);
      finalvarr = gg1 * cov2 * t(gg1);
      finalsdd = sqrt(finalvarr);

      postprobs=replicate||id||rmst1||rmst0||ss1||ss0||cureprop1||cureprop0;
      cname = {"replicate" "id" "rmst1" "rmst0" "s1" "s0" "pi1" "pi0"};

      create final21 from postprobs  [ colname=cname ];
      append from postprobs;

      postprobs=replicate[1]||apv||apvv||finalvar||finalvarr||finalsd||finalsdd;
      cname = {"replicate" "apv" "apvv" "var" "varr" "sd" "sdd"};

      create final3 from postprobs  [ colname=cname ];
      append from postprobs;
      quit;
    run;
  %end;

  /* If assessment is greater than last observed time */

  %if &teind = 2 %then %do;

  	%put 2;
	proc iml;
	  use temp261;
      read all into cov;
      use aa43;
      read all var{estimate} into survest;
      use aa42;
      read all var{estimate} into cureest;
      use aa44;
      read all var{estimate} into baseest;
      use temp3r;
      read all var{SURVtime} into int;
      ntime = nrow(int);
      max_time = &rmst;

	  %if &stat = 1 %then %do;
	    %put 1;
	    use _dat05r;
        read all var{id} into id;
        read all var{percent} into per;
        read all var{replicate} into replicate;
	    nobs = replicate[,1];
        obs = nrow(nobs);
        read all var{&SURVCov} into SURVcov;
        read all var{&CURECov} into CUREcov;

        cure1 = j(obs,1,1)|| j(obs,1,1)|| CUREcov;
        cure0 = j(obs,1,1)|| j(obs,1,0)|| CUREcov;
        surv1 = j(obs,1,1)|| SURVcov;
        surv0 = j(obs,1,0)|| SURVcov;
      %end;

      %else %if &stat = 2 %then %do;
	    %put 2;
	    use _dat05r;
        read all var{id} into id;
        read all var{percent} into per;
        read all var{replicate} into replicate;
	    nobs = replicate[,1];
        obs = nrow(nobs);
        read all var{&SURVCov} into SURVcov;

        cure1 = j(obs,1,1)|| j(obs,1,1);
        cure0 = j(obs,1,1)|| j(obs,1,0);
        surv1 = j(obs,1,1)|| SURVcov;
        surv0 = j(obs,1,0)|| SURVcov;
      %end;

      %else %if &stat = 3 %then %do;
	    %put 3;
	    use _dat05r;
        read all var{id} into id;
        read all var{percent} into per;
        read all var{replicate} into replicate;
	    nobs = replicate[,1];
        obs = nrow(nobs);
        read all var{&CURECov} into CUREcov;

        cure1 = j(obs,1,1)|| j(obs,1,1)|| CUREcov;
        cure0 = j(obs,1,1)|| j(obs,1,0)|| CUREcov;
        surv1 = j(obs,1,1);
        surv0 = j(obs,1,0);
      %end;

      %else %if &stat = 4 %then %do;
	    %put 4;
	    use _dat05r;
        read all var{id} into id;
        read all var{percent} into per;
        read all var{replicate} into replicate;
	    nobs = replicate[,1];
        obs = nrow(nobs);

        cure1 = j(obs,1,1)|| j(obs,1,1);
        cure0 = j(obs,1,1)|| j(obs,1,0);
        surv1 = j(obs,1,1);
        surv0 = j(obs,1,0);
      %end;

      /* remove missing covariance value in covariance matrix */
      cov1 = cov[loc(cov[1,] ^=.),];
      cov2 = cov1[,loc(cov[1,] ^=.)];

      int1 = 0//int[1:(nrow(int))];
      int2 = int[1:(nrow(int))]//&rmst;
      interv = int2 - int1;

      min_b = j(obs, 1, 600);
      temp1 = j(obs, ntime, .);
      tem1 = j(obs, (ntime+1), .);
      te1 = j(obs, ntime, .);
      do i = 1 to ntime;
        temp1[,i] = exp(min_b >< (j(obs, 1, baseest[i]) + surv1 * survest));
	    te1[,i] = exp(min_b >< (j(obs, 1, baseest[i]) + surv1 * survest));
	  end;
      do p = 2 to ntime;
	    temp1[,p] = temp1[,p] + temp1[,p-1];
	  end;
	  /* Add one column 1 for Restricted Mean Survival Time Calculation */
	  tem1 = j(obs, 1, 0)||temp1[, 1:ntime];
	  s1 = j(obs, (ntime+1), .);
	  s1 = exp((-1) * tem1);
      /* No need to shift for survival probability at Tmax Calculation at step down */
	  ss1 = j(obs, 1, .);
      /* No need to shift for survival probability at Tmax Calculation at step up */
      ss1 = (exp((-1) * tem1))[, (ntime+1)];

      temp0 = j(obs, ntime, .);
      tem0 = j(obs, (ntime+1), .);
      te0 = j(obs, ntime, .);
      do m = 1 to ntime;
      	temp0[,m] = exp(min_b >< (j(obs, 1, baseest[m]) + surv0 * survest));
    	te0[,m] = exp(min_b >< (j(obs, 1, baseest[m]) + surv0 * survest));
  	  end;
      do n = 2 to ntime;
	    temp0[,n] = temp0[,n] + temp0[,n-1];
	  end;
      /* Shift one column for Restricted Mean Survival Time Calculation */
	  tem0 = j(obs, 1, 0)||temp0[, 1:ntime];
	  s0 = j(obs, (ntime+1), .);
	  s0 = exp((-1) * tem0);
      /* No need to shift for survival probability at Tmax Calculation */
	  ss0 = j(obs, 1, .);
      /* No need to shift for survival probability at Tmax Calculation at step up */
      ss0 = (exp((-1) * tem0))[, (ntime+1)];

      rmst1 = j(obs, 1, .);
      rmst0 = j(obs, 1, .);
      rmst1 = s1 * interv;
      rmst0 = s0 * interv;

      cureprop1 = exp(min_b >< (cure1 * cureest)) /(1 + exp(min_b >< (cure1 * cureest)));
      cureprop0 = exp(min_b >< (cure0 * cureest)) /(1 + exp(min_b >< (cure0 * cureest)));
      /* Overall Exposure Effect on RMST */
      exp = j(obs, 1, .);
      exp = (((1 - cureprop1) * max_time + cureprop1 # rmst1) - ((1 - cureprop0) * max_time + cureprop0 # rmst0))#per/100;
      apv = sum(exp); 
      /* Overall Exposure Effect on survival probability at Tmax */
      expp = j(obs, 1, .);
      expp = (((1 - cureprop1) + cureprop1 # ss1) - ((1 - cureprop0) + cureprop0 # ss0))#per/100;
      apvv = sum(expp); 
      print apv apvv;

      ncure = nrow(cureest);
      nsurv = nrow(survest);
      nbase = nrow(baseest);

      k = nrow(survest) + nrow(cureest) + nrow(baseest);
      /* Delta method for variance calculation for exposure effect on RMST */
      g = j(1,k,0);
      /* Delta method for variance calculation for exposure effect on probability at Tmax */
      gg = j(1,k,0);
      lincure1 = cure1 * cureest;
      lincure0 = cure0 * cureest;
      linsurv1 = surv1 * survest;
      linsurv0 = surv0 * survest;

      do i = 1 to obs; 
        t11 = exp(min(600, lincure1[i]))/(1+exp(min(600, lincure1[i])))**2 * cure1[i,];
        t10 = exp(min(600, lincure0[i]))/(1+exp(min(600, lincure0[i])))**2 * cure0[i,];
        g[1:ncure] = g[1:ncure] + t((((-1) * max_time * t11 + rmst1[i] * t11) - ((-1) * max_time * t10 + rmst0[i] * t10))*per[i]/100);
        gg[1:ncure] = gg[1:ncure] + t((((-1) * t11 + ss1[i] * t11) - ((-1) * t10 + ss0[i] * t10))*per[i]/100);
        t21 = (s1 # tem1 * (-1) * interv);
        t20 = (s0 # tem0 * (-1) * interv);
        g[(ncure + 1):(ncure + nsurv)] = g[(ncure + 1):(ncure + nsurv)] + t((cureprop1[i] * t21[i] * surv1[i,] - cureprop0[i] * t20[i] * surv0[i,])*per[i]/100);
        tt21 = (ss1 # temp1[, ntime] * (-1));
        tt20 = (ss0 # temp0[, ntime] * (-1));
        gg[(ncure + 1):(ncure + nsurv)] = gg[(ncure + 1):(ncure + nsurv)] + t((cureprop1[i] * tt21[i] * surv1[i,] - cureprop0[i] * tt20[i] * surv0[i,])*per[i]/100);

        /* Generate Triangle Matrix */
        temp11 = j(ntime, ntime, 1);
        r = row(temp11);      /* create helper matrices */
        c = col(temp11);
        lowerIdx = loc(r > c);
        temp11[lowerIdx] = 0;

        temp12 = t(t(s1[i,2:(ntime + 1)]) * j(1, ntime, 1));
        temp13 = t(te1[i,]) * j(1, ntime, 1) * (-1);
        temp14 = temp11 # temp12 # temp13;
        temp15 = j(ntime, 1, 0)||temp14;
        t31 = temp15 * interv;

        temp02 = t(t(s0[i,2:(ntime + 1)]) * j(1, ntime, 1));
        temp03 = t(te0[i,]) * j(1, ntime, 1) * (-1);
        temp04 = temp11 # temp02 # temp03;
        temp05 = j(ntime, 1, 0)||temp04;
        t30 =  temp05 * interv;
        g[(ncure + nsurv + 1):k] = g[(ncure + nsurv + 1):k] + ((cureprop1[i] * t31 - cureprop0[i] * t30)*per[i]/100);

        tt31 = t(ss1[i] * te1[i, ] * (-1));
        tt30 = t(ss0[i] * te0[i, ] * (-1));
        gg[(ncure + nsurv + 1):k] = gg[(ncure + nsurv + 1):k] + ((cureprop1[i] * tt31 - cureprop0[i] * tt30)*per[i]/100);
      end;
      g1 = t(g[loc(cov[1,] ^=.)]);
      gg1 = t(gg[loc(cov[1,] ^=.)]);
      finalvar = g1 * cov2 * t(g1);
      finalsd = sqrt(finalvar);
      finalvarr = gg1 * cov2 * t(gg1);
      finalsdd = sqrt(finalvarr);

	  postprobs=replicate||id||rmst1||rmst0||ss1||ss0||cureprop1||cureprop0;
      cname = {"replicate" "id" "rmst1" "rmst0" "s1" "s0" "pi1" "pi0"};

      create final21 from postprobs  [ colname=cname ];
      append from postprobs;

      postprobs=replicate[1]||apv||apvv||finalvar||finalvarr||finalsd||finalsdd;
      cname = {"replicate" "apv" "apvv" "var" "varr" "sd" "sdd"};

      create final3 from postprobs  [ colname=cname ];
      append from postprobs;

      quit;
    run;
  %end;

  data out31;
    set out31 temp32 (in = a);
    if a then replicate = &i - 1;
	drop _name_;
  run;

  data out32;
    set out32 temp25 (in = a);
    if a then replicate = &i - 1;
  run;

  data out33;
    set out33 final21 (in = a);
    if a then replicate = &i - 1;
  run;

  data out34;
    set out34 final3 (in = a);
    if a then replicate = &i - 1;
  run;

  data out35;
    set out35 temp28 (in = a);
    if a then replicate = &i - 1;
  run;

  proc datasets lib = work;
    delete temp95r temp64r temp60r temp43r tt71 tt72 base1 base11 base12 base2 base3
    temp25 temp28 _dat05r temp3r aa42 aa43 aa44 temp261 temp63r init  
    temp26 temp27 temp29 temp30 temp32 final3 final21 mm00;
  run;
%end;

data final2;
  set out21(in = b) out31 (in = c) out35(in = d);
  format method $10.;
  if b then method = 'Select';
  if c then method = 'Znonpara';
  if d then method = 'Znonpara';
  if replicate ne -1;
run;

data final3;
  set out22(in = b) out32 (in = c);
  format method $10.;
  if b then method = 'Select';
  if c then method = 'Znonpara';
  if replicate ne -1;
run;

data final4;
  set out23(in = b) out34 (in = c);
  format method $10.;
  if b then method = 'Select';
  if c then method = 'Znonpara';
  if replicate ne -1;
run;

data final41;
  set out23;
  if replicate = 0;
run;

proc univariate data = out23 noprint;
  var exp expp;
  output out=final42 pctlpts=2.5 97.5  pctlpre = exp expp pctlname = _025 _975;
  where replicate ne 0;
run;

data final43;
  merge final41 final42;
run;

data final44;
  set out34;
  if replicate ne -1;
  exp_025 = apv - 1.959964 * sd;
  exp_975 = apv + 1.959964 * sd;
  rename apv = exp;
  rename sd = std_exp;
  expp_025 = apvv - 1.959964 * sdd;
  expp_975 = apvv + 1.959964 * sdd;
  rename apvv = expp;
  rename sdd = std_expp;
run;

data final45;
  set final43 (in = b) final44 (in = c);
  if b then method = 'FP or RCS Parametric Approximation';
  if c then method = 'Piecewise Linear Approximation';
  keep method exp exp_025 exp_975 expp expp_025 expp_975;
run;

data temp67;
  set temp66;
  if replicate = 0;
  t_access = &t_access;
  if t_access ne 999999 and t_access gt max_time then sta = 1;
  else sta = 0;
run;

data _null_;
  set temp67;
  call symput('stata', sta);
run;

%if &stata = 1 %then %do;
  proc iml;
    print,"WARNING: THE USER REQUIRED EXPOSURE EFFECT ASSESSMENT AFTER THE LAST OBSERVED TIME (NOT SUGGESTED, AND THE REASON IS THAT BASELINE DISTRIUBTION APPROXIMATION WHICH IS REQUIRED FOR THE EXPOSURE EFFECT ESTIMATION MAY NOT BE ACCURATE, SEE REFERENCE FOR DETAILS)";
  quit;
%end;

proc datasets lib = work;
  save final2 final3 final4 final5 final6 final45;
run; 

data &out;
  retain method exp exp_025 exp_975;
  set final45;
  drop expp expp_025 expp_975;
run;

data final51;
set final5;
if replicate = 0;
keep model parameter effect;
rename effect = variable;
run;

data final61;
set final6;
if replicate = 0;
model = 'Cure';
keep model parameter variable;
run;

data final7;
set final61 final51 ;
run;

proc sort data = final7;
by parameter;
run;

data final31;
set final3;
keep id method parameter estimate standarderror probt;
if replicate = 0 and substr(parameter, 1, 1) in ('a', 'b');
run;

proc sort data = final31;
by parameter;
run;

data final32;
merge final31 final7;
by parameter;
if method = 'Select' then method1 = 'FP or RCS Parametric Approximation';
if method = 'Znonpara' then method1 = 'Piecewise Linear Approximation';
if id ge 1 and id le 44 then select_model = 'FP Parametric';
if id = 45 then select_model = 'RCS 1 knot';
if id = 46 then select_model = 'RCS 2 knos';
if id = 47 then select_model = 'RCS 3 knots';
run;

proc sort data = final32;
by method parameter;
run;

data final33;
retain method1 select_model model variable; 
set final32;
drop id method;
rename method1 = method;
run;

ods listing;

%if &printall=T %then %do;

/*  dm  'log;clear;out;clear;';*/

  Title 'Overall Exposure Effects Estimates';

  proc print data = &out noobs label;
    label exp = 'Exposure Effect on RMST';
    label exp_025 = '95% Lower Limit for RMST Exposure';
    label exp_975 = '95% Upper Limit for RMST Exposure';
/*    label expp = 'Exposure Effect on Survival Probability';*/
/*    label expp_025 = '95% Lower Limit for Survival Probability Exposure';*/
/*    label expp_975 = '95% Upper Limit for Survival Probability Exposure';*/
	var method exp exp_025 exp_975;
  run;

Title 'Cox Proportional Hazards Cure Models Parameter Estimates';

  proc print data = final33 noobs label;
  run;

%end;  *-end printall=T option;

%end;

%mend;

data ref;
input sexF Duration F10Cigs nsub;
datalines;
1 10 10 10
0 10 10 10
;
run;

proc datasets lib = work;
  save ref;
run; 

/*libname te 'M:\Various\Wei\SMMR Sub Sep 2019\Program';*/
/**/
/*data aa11;*/
/*set te.cessation;*/
/*drop zip;*/
/*if _N_ le 200;*/
/*run;*/
/**/
/*proc export */
/*  data=work.aa11 */
/*  dbms=xlsx */
/*  outfile="M:\Various\Wei\SMMR Sub Sep 2019\Program\example.xlsx" */
/*  replace;*/
/*run;*/

proc import 
  datafile="M:\Various\Wei\SMMR Sub Sep 2019\Program\example.xlsx" 
  dbms=xlsx 
  out=work.aa11 
  replace;
run;

libname temp 'M:\Various\Wei\SMMR Sub Sep 2019\Program\Output';

%CUREEXP   (data             = aa11, 
            subjid           = ObsLHS,
            CUREcov          = sexF Duration F10Cigs,
            SURVcov          = sexF Duration F10Cigs,
            EXP              = SI_UC, 
            nboot            = 1, 
            SURVLtime        = timept1,
            SURVRtime        = timept2,
            SURVevent        = relapse,
			T_ACCESS         = 999999,
			REFERENCE        = 1,
            seedboot         = 1617893,
            Printall         = T,
            Out             = out02
            );

ods listing;
