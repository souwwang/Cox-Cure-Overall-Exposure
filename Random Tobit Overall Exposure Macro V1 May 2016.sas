************************************************************************************************** ;
*	ESTIMATING OVERALL EXPOSURE EFFECTS FOR RANDOM EFFECT TOBIT REGRESSION MODELS                  ;
*   Version 1, Uploaded 05/06/2016           										               ;
*																					               ;
*	Code written by Wei Wang (2016)													               ;
* 																					               ;
*	Reference:																		               ;
* 	Wang W, and Griswold ME (2016). Estimating overall exposure effects for the clustered and      ; 
*   censored outcome using random effect Tobit regression models.                                  ;
* 	Statistics in Medicine. Submitted.	                                    		               ;
* 																					               ;
***************************************************************************************************;
*																					               ;
*	PROGRAM DESCRIPTION:															               ;
*   -----------------------------------------------------------------------------------------------;
*	THIS PROGRAM PRODUCES THE ESTIMATION OF OVERALL EXPSOURE EFFECTS FOR RANDOM EFFECT TOBIT       ; 
*	REGRESSION MODELS. IT CAN ACCOMODATE SINGLE BOUNDARY AND DOUBLE BOUNDARY                       ; 
*   RESPONSE VARIABLES AND ALLOW MULTIPLE COVARIATES IN THE DATA SET. THIS MACRO                   ;
*   WILL USE TWO MAJOR APPROACHES TO ESTIMATE THE OVERALL EXPOSURE EFFECTS FOR RANDOM EFFECT TOBIT ;
*   MODELS, DIRECT-MARGINALIZATION APPROACH AND AVERAGE PREDICTED VALUE APPROACH. THE              ;
*   DIRECT-MARGINALIZATION APPROACH AND AVERAGE PREDICTED VALUE APPROACH ASSUME DIFFERNT MODEL     ;
*   SETTINGS, IN WHICH THE FORMER ASSUMES HOMOGENUOUS EXPOSURE AND COVARIATE                       ;
* 	EFFECTS FOR THE MARGINAL MEANS																   ;
*	AND THE LATTER ASSUMES HOMOGENUOUS EXPOSURE AND COVARIATE EFFECTS FOR THE LATENT DEPENDENT     ;
*   VARIABLE GIVEN THE SPECIFIED VALUE OF THE RANDOM EFFECT. GOODNESS-OF-FIT STATISTICS AIC WILL BE; 
*   PROVIDED FOR EACH METHOD, AND THE OVERALL EXPOSURE EFFECTS WITH LEAST AIC                      ;
*	IS SUGGESTED TO BE CHOSEN. FOR AVERAGE PREDICTED VALUE APPROACH, EMPIRICAL DISTRIBUTION        ;
*	OF THE BASELINE COVARIATES FROM THE TARGET OR REFERENCE POPULATION ARE USED FOR                ;  
*   THE OVERALL EXPOSURE EFFECTS ESTIMATION. THE EXPOSED GROUP                                     ;
*	WAS USEED AS TARGET POPULATION IN THIS MACRO.                                                  ;
*	THIS MACRO IS MAINLY TO ASSESS THE OVEARALL EXPOSURE EFFECTS FOR EXPOSED VS.                   ;
*	NON-EXPOSED GROUP WITH  CENSORED AND CLUSTER-CORRELATED DATA. SIMPLE RANDOM INTERCEPT IS       ;
*	ASSUMED IN TOBIT MODEL TO ACCOUNT FOR CLUSTERING ISSUE IN THE DATA     					       ;
*																					               ;
*	DATASET REQUIREMENTS:															               ;
*   -----------------------------------------------------------------------------------------------;
*	THE DATASET, MUST HAVE ONE LINE PER SUBJECT WHERE EACH SUBJECT MUST CONTAIN 	               ;
*	ONE BINARY EXPOSURE VARIABLE, ANY NUBMER OF COVARIATES (ZERO TO INFINITY)                      ;
*   AND ONE TRUNCATED NORMALLY DISTRIBUTED OUTCOME WITH INFLATION AT BOUNDARY VALUE.               ;
*   OF NOTE, BASELINE COVARIATES CAN ONLY BE BINARY OR CONTINUOUS, SO CATEGORICAL                  ;
*   BASELINE COVARIATES WITH MORE THAN TWO CATEGORIES SHOULD BE TRANSFORMED                        ;
*   TO DUMMY VARIABLES FIRST BEFORE APPLYINGCAN THIS MACRO. FOR EXAMPLE, CATEGORICAL               ;
*   COVARIATES WITH THREE DISTINCT VALUES SHOULD BE TRANSFORMED TO TWO DUMMY VARIABLES             ;
*   BEFORE THE CURRENT MACRO CAN BE USED.                                                          ;
*	THIS MACRO CAN ONLY WORK FOR DATASET WITHOUT MISSING VALUES. THE MISSING VALUES SHOULD BE      ;  
*   IMPUTED OR ANY SUBJECT WITH MISSING DATA SHOULD BE DELETED. THE SAMPLE                         ;
*	DATA SET WITH TWO COVARIATES AND OUTCOME WITH LOWER MEASUREMENT LIMIT 3.0 IS                   ;
*	LISTED AS FOLLOWING.																	       ;
*																					               ;
*			CLUSTER_ID SUBJECT  EXPOSURE COVARIATE1 COVARIATE2 OUTCOME	   			               ;
*				 1        1	    	1		21		    1       13.2        		               ;
*				 1        2		    0		28		    0        3.0         		               ;
*				 1        3		    1		25		    1       12.5	       		               ;
*				 1        4		    0		24		    1        3.0	         	               ;
*				 2        5		    0		26		    1       10.1	      	                   ;
*				 2        6		    0		32		    1        7.8	    	                   ;
*				 2        7		    1		16		    1        6.5	    	                   ;
*																					               ;
*																					               ;
*	MODEL SPECIFICATION																               ;
*   -----------------------------------------------------------------------------------------------;
*	T: BINARY EXPOSURE Y*: NORMALLY DISTRIBUTED LATENT DEPENDENT VARIALBE WITHOUT                  ;
*   BOUNDARY Y: TRUNCATED NORMALLY DISTRIBUTED OUTCOME WITH INFLATION AT BOUNDAY                   ;
*   VALUE W: BASELINE COVARIATE. Y* CAN ONLY BE OBSERVED IF WITHIN DETECTABLE RANGE                ;
*																					               ;
*   LATENT RANDOM EFFECT TOBIT REGRESSION MODEL									                   ;
*	Y* = BETA0 + BETA1 * W + BETA2 * T + B + EPSILON         B IS RANDOM EFFECT                    ;
*	WHERE B~N(0, SIGMAUSQUARE), EPSILON~N(0, SIGMASQUARE)							               ;
*																					               ;
*   MARGINAL RANDOM EFFECT TOBIT REGRESSION MODEL								                   ;
*	E(Y) = GAMMA0 + GAMMA1 * W + GAMMA2 * T                                                        ;
*																					               ;
*													                 				               ;
*	MACRO VARIABLES:																               ;
*   -----------------------------------------------------------------------------------------------;
*																					               ;
*	DATASET:        MAIN INPUT DATASET												               ;
*																					               ;
*	Y:		SHOULD BE CONTINUOUS VARIABLE WITH OR WITHOUT BOUNDAIES              	               ;
*																					               ;
*	X:		COVARIATE VECTOR INCLUDING BOTH BASELINE COVARIATES AS WELL AS BINARY                  ; 
*           EXPOSURE VARIABLE. OF NOTE, THE MACRO ALLOWS ANY NUMBER OF BASELINE                    ; 
*           COVARIATES, BUT CATEGORICAL COVARIATES WITH MORE THAN TWO DISTINCT                     ; 
*           VALUES SHOULD BE TRANSFORMED TO DUMMY VARIABLES FIRST. AND ALSO BASELINE               ;
*			COVARIATES SHOULD BE SPECIFIED FIRST AND LAST VARIABLE IN COVARIATE                    ;
*           VECTOR X SHOULD BE THE BINARY EXPOSURE VARIABLE.						               ;
*																					               ;
*	CUTOFF1:LOWER DETECTION LIMIT IF APPLICABLE 									               ;
*																					               ;
*	CUTOFF2:UPPER DETECTION LIMIT IF APPLICABLE 									               ;
*																					               ;
*	OF NOTE, IF CUTOFF1 OR CUTOFF2 DOES NOT EXIST (SINGLE BOUNDARY RESPONSE VARIABLE               ;
*   ), ONLY THE AVAILABLE BOUNDARY NEEDS TO BE SPECIFIED, AND THE ONE THAT DOES NOT                ;
*   EXIST CAN BE OMITTED (SEE EXAMPLE BELOW).                                                      ; 
*																					               ;
*	CLUSTER: CLUSTER ID																               ;
*																					               ;
*	OUT: NAME OF THE OUTPUT DATASET CONTAINING THE OVERALL EXPOSURE EFFECTS ESTIMATE               ;
*																					               ;
*	OF NOTE, WHEN THE MAGNITUDE OF SOME VARIABLES ARE LARGE, THE ESTIMATION METHOD IN THIS MACRO   ;
*	MAY BE EASY TO FAIL BY REACHING ABSGCONV CONVERGENCE CRITERION WITH STARTING VALUES. WHEN THIS ;
*	PROBLEM HAPPENS, ONE METHOD IS TO RESCALE THE INPUT VARIABLES BEFORE APPLYING THE MACRO.       ;
*	FOR EXAMPLE, ONE BASELINE COVARIATE WITH MEAN, 100 STD 10 CAN BE DIVIDED BY 100 BEFORE APPLYING;
*	THE MACRO. WITH THE TRANSFORMED VARIABLE, THE FINAL COEFFICIENT SHOULD BE TRANSFORMED BACK TO  ;
*	ORIGINAL SCALE BY DIVIDING 100.      											               ;
*																					               ;
*	EXAMPLE CODE:																	               ;
*																					               ;
*   %include 'TOBITCORR.sas' 														           ;
*																					               ;
*	DOUBLE BOUNDARY:    															               ;
*   %TOBITCORR(dataset=DSN, y=RESPONSE, x=COVARIATE1 COVARIATE2 EXPOSURE,                          ; 
*   cutoff1 = 0, cutoff2 =100, cluster = cluster_id, out =OUT1)                                    ;
*																					               ;
*	SINGLE BOUNDARY:    															               ;
*   %TOBITCORR(dataset=DSN, y=RESPONSE, x=COVARIATE1 COVARIATE2 EXPOSURE,                          ; 
*   cutoff1 = 0, cluster = cluster_id, out =OUT1)                                                  ;
*																					               ;
***************************************************************************************************;

***************************************************************************************************;
*************************************Final Macro***************************************************;
***************************************************************************************************;

%macro tobitcorr(dataset=' ', y=' ', x=' ', cutoff1 =., cutoff2 =., cluster = '', out =' ');

/* No Boundary using 40-point Gaussian-Hermite Quadrature*/

%macro tobitcorr0(dsn=, y=,x=, cutoff1=, cutoff2 =, theta= 0);
proc iml;
reset noname;
use &dsn;
read all var {&x} into x;
read all var {&y} into y;
yx=y || x;
* a = lower truncation point;
use tt3;
read all into est;
use abswei;
read all var {abs weight} into w;
limit = y;

start LL(theta) global (yx, limit, w);
* split data into censored and uncensored ohs;
yx0 = yx[loc((limit ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

k = ncol(theta) - 1;
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma * sigma;
sigmau = theta [k + 1];
sigmau2 = sigmau * sigmau;
sum1 = 0;
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do e = 1 to n0;
   mu = x0[e,]*beta;
   yi = y0[e];
   chk3 = t(1/sqrt(constant('pi')) /sqrt(2*constant('pi')*sigma**2) * exp((-1)*(yi - mu - sqrt(2) * sigmau * w[,1])
             #(yi - mu - sqrt(2) * sigmau * w[,1])/(2*sigma**2)))* w[,2];
   if chk3 <= 0 then chk3 = 0.1e8;
   sum1 = sum1 + log(chk3);
end;

f = sum1;
return(f);
finish LL ;

optn={1 2};
*tc={3000 15000 . 1E-12};
tc={3000 15000 . 1E-6};
*tc={3000 15000};
* define starting va1ues for theta: either user-provided or O1S;
theta0 = &theta;
if theta0 = 0 then
do;
theta0 = t(est);
print 'Mixed Effects Model ignoring Boundary starting values for theta = '
theta0;
end;
else print 'User-supplied starting values for
theta = ' theta0;
t1 = j(1,ncol(theta0),.);
t2 = t1;
t2[ncol(theta0)] = 0;
t2[ncol(theta0)-1] = 0;
t3 = t2//t1;
con = t3;
call nlpqn(rc, theta, 'LL',theta0,optn, con, tc);
*call nlpqn(rc, theta, 'LL',theta0,optn, con, tc) grd='GRAD';
call nlpfdd(f, g, h,'LL', theta);
var = inv(-h);
sd = t(j(1,ncol(theta0),.));
temp = vecdiag(var);
do i = 1 to nrow(temp);
if temp[i] >= 0 then sd[i] = sqrt(temp[i]);
else sd[i] = .;
end; 
*sd = sqrt(vecdiag(var));
print 'Hessian = ' h ' Covariance matrix = ' var
' Standard errors = ' sd;
*print 'Hessian = ' h ' Covariance matrix = ' var
;
*print rc theta;

yx0 = yx[loc((limit ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

k = ncol(theta) - 1;
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma * sigma;
sigmau = theta [k + 1];
sigmau2 = sigmau * sigmau;
sum1 = 0;
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do e = 1 to n0;
   mu = x0[e,]*beta;
   yi = y0[e];
   chk3 = t(1/sqrt(constant('pi')) /sqrt(2*constant('pi')*sigma**2) * exp((-1)*(yi - mu - sqrt(2) * sigmau * w[,1])
             #(yi - mu - sqrt(2) * sigmau * w[,1])/(2*sigma**2)))* w[,2];
   if chk3 <= 0 then chk3 = 0.1e8;
   sum1 = sum1 + log(chk3);
end;

like = sum1*(-2);
aic = sum1*(-2) + 2 * ncol(theta);
bic = ncol(theta) * log(nrow(yx0)) + sum1*(-2);

sum1 = theta[k-1];
finalsd = sd[k-1];

postprobs=like||aic||bic||sum1||finalsd||rc||theta||t(sd);
 cname = {"RLike" "RAIC" "RBIC" "REST" "RSTD" "RSTA" };

create ff11 from postprobs  [ colname=cname ];
append from postprobs;

quit;
run;
%mend;

/* Two Boundary APV Method using 40-point Gaussian-Hermite Quadrature, T = 1 as reference*/

%macro tobitcorr1(dsn=, y=,x=, cutoff1=, cutoff2 =, theta= 0);
proc iml;
reset noname;
use &dsn;
read all var {&x} into x;
read all var {&y} into y;
yx=y || x;
* a = lower truncation point;
use tt3;
read all into est;
use abswei;
read all var {abs weight} into w;
a = &cutoff1;
b = &cutoff2;
* define var indicating censored ohs;
limit = y;
limit1 = y;
limit2 = y;
limit1 [loc(y=a)] = .;
limit2 [loc(y=b)] = .;
limit3 = limit1;
limit3 [loc(limit1=b)] = .;
n0 = nrow(y);
n1 = ncol(loc( limit1=.));
n2 = ncol(loc( limit2=.));
print n1 "observations are censored at &cutoff1 and" n2 "observations are censored at &cutoff2 out of " n0 " obs";
* define log likelihood function;
*print limit1, limit2, limit3;

start LL(theta) global (yx, limit, limit1, limit2, limit3, a, b, w);
* split data into censored and uncensored ohs;
yx0 = yx[loc((limit3 ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx1 = yx[loc(( limit1 = .)),];
n1 = nrow(yx1);
y1 = yx1[,1];
x1 = j(n1,1,1) || yx1[,2:k];

yx2 = yx[loc(( limit2 = .)),];
n2 = nrow(yx2);
y2 = yx2[,1];
x2 = j(n2,1,1) || yx2[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta) - 1;
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma * sigma;
sigmau = theta [k + 1];
sigmau2 = sigmau * sigmau;
sum1 = 0;
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do i = 1 to n1;
   mu = x1[i,]*beta;
   chk1 = t(1/sqrt(constant('pi')) * probnorm((a - mu - sqrt(2) * sigmau * w[,1])/sigma))* w[,2];
   if chk1 <= 0 then chk1 = 0.1e8;
   sum1 = sum1 + log(chk1);
end;

do j = 1 to n2;
   mu = x2[j,]*beta;
   chk2 = t(1/sqrt(constant('pi')) * probnorm((mu + sqrt(2) * sigmau * w[,1] - b)/sigma))* w[,2];
   if chk2 <= 0 then chk2 = 0.1e8;
   sum1 = sum1 + log(chk2);
end;

do e = 1 to n0;
   mu = x0[e,]*beta;
   yi = y0[e];
   chk3 = t(1/sqrt(constant('pi')) /sqrt(2*constant('pi')*sigma**2) * exp((-1)*(yi - mu - sqrt(2) * sigmau * w[,1])
             #(yi - mu - sqrt(2) * sigmau * w[,1])/(2*sigma**2)))* w[,2];
   if chk3 <= 0 then chk3 = 0.1e8;
   sum1 = sum1 + log(chk3);
end;

f = sum1;
return(f);
finish LL ;

optn={1 2};
*tc={3000 15000 . 1E-12};
tc={3000 15000 . 1E-6};
*tc={3000 15000};
* define starting va1ues for theta: either user-provided or O1S;
theta0 = &theta;
if theta0 = 0 then
do;
theta0 = t(est);
print 'Mixed Effects Model ignoring Boundary starting values for theta = '
theta0;
end;
else print 'User-supplied starting values for
theta = ' theta0;
t1 = j(1,ncol(theta0),.);
t2 = t1;
t2[ncol(theta0)] = 0;
t2[ncol(theta0)-1] = 0;
t3 = t2//t1;
con = t3;
call nlpqn(rc, theta, 'LL',theta0,optn, con, tc);
*call nlpqn(rc, theta, 'LL',theta0,optn, con, tc) grd='GRAD';
call nlpfdd(f, g, h,'LL', theta);
var = inv(-h);
sd = t(j(1,ncol(theta0),.));
temp = vecdiag(var);
do i = 1 to nrow(temp);
if temp[i] >= 0 then sd[i] = sqrt(temp[i]);
else sd[i] = .;
end; 
*sd = sqrt(vecdiag(var));
print 'Hessian = ' h ' Covariance matrix = ' var
' Standard errors = ' sd;
*print 'Hessian = ' h ' Covariance matrix = ' var
;
*print rc theta;

yx0 = yx[loc((limit3 ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx1 = yx[loc(( limit1 = .)),];
n1 = nrow(yx1);
y1 = yx1[,1];
x1 = j(n1,1,1) || yx1[,2:k];

yx2 = yx[loc(( limit2 = .)),];
n2 = nrow(yx2);
y2 = yx2[,1];
x2 = j(n2,1,1) || yx2[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta) - 1;
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma * sigma;
sigmau = theta [k + 1];
sigmau2 = sigmau * sigmau;
sum1 = 0;
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do i = 1 to n1;
   mu = x1[i,]*beta;
   chk1 = t(1/sqrt(constant('pi')) * probnorm((a - mu - sqrt(2) * sigmau * w[,1])/sigma))* w[,2];
   if chk1 <= 0 then chk1 = 0.1e8;
   sum1 = sum1 + log(chk1);
end;

do j = 1 to n2;
   mu = x2[j,]*beta;
   chk2 = t(1/sqrt(constant('pi')) * probnorm((mu + sqrt(2) * sigmau * w[,1] - b)/sigma))* w[,2];
   if chk2 <= 0 then chk2 = 0.1e8;
   sum1 = sum1 + log(chk2);
end;

do e = 1 to n0;
   mu = x0[e,]*beta;
   yi = y0[e];
   chk3 = t(1/sqrt(constant('pi')) /sqrt(2*constant('pi')*sigma**2) * exp((-1)*(yi - mu - sqrt(2) * sigmau * w[,1])
             #(yi - mu - sqrt(2) * sigmau * w[,1])/(2*sigma**2)))* w[,2];
   if chk3 <= 0 then chk3 = 0.1e8;
   sum1 = sum1 + log(chk3);
end;

like = sum1*(-2);
aic = sum1*(-2) + 2 * ncol(theta);
bic = ncol(theta) * log(nrow(yx0)) + sum1*(-2);

k = ncol(theta) - 1;
k1 = k-1;
beta = theta [1:k1];
x1 = x[loc((x[,k-2] = 1)),];
n1 = nrow(x1);
if k >= 4 then x0 = j(n1,1,1) || x1[,1:(k-3)]||j(n1,1,0);
else if k = 3 then x0 = j(n1,1,1) ||j(n1,1,0);
x1 = j(n1,1,1) || x1;
sigma = theta [k];
sigma2 = sigma * sigma;
sigmau = theta [k + 1];
sigmau2 = sigmau * sigmau;
sum1 = 0;
g = j(1,k+1,0);

do i = 1 to n1;
   mu = x1[i,]*beta;
   s1 = probnorm((b - (mu + sqrt(2) * sigmau * w[,1]))/sigma);
   s2 = probnorm((a - (mu + sqrt(2) * sigmau * w[,1]))/sigma);
   s3 = pdf ('normal', (b - (mu + sqrt(2) * sigmau * w[,1]))/sigma, 0, 1);
   s4 = pdf ('normal', (a - (mu + sqrt(2) * sigmau * w[,1]))/sigma, 0, 1);
   ss1 = (-1) * s3/sigma;
   ss2 = (-1) * s4/sigma;
   ss3 = s3 # (b - (mu + sqrt(2) * sigmau * w[,1])) / sigma**2;
   ss4 = s4 # (a - (mu + sqrt(2) * sigmau * w[,1])) / sigma**2;
   sss1 = (-1) * (b - (mu + sqrt(2) * sigmau * w[,1])) # s3/sigma**2;
   sss2 = (-1) * (a - (mu + sqrt(2) * sigmau * w[,1])) # s4/sigma**2;
   sss3 = s3 # (b - (mu + sqrt(2) * sigmau * w[,1]))# (b - (mu + sqrt(2) * sigmau * w[,1]))/ sigma**3;
   sss4 = s4 # (a - (mu + sqrt(2) * sigmau * w[,1]))# (a - (mu + sqrt(2) * sigmau * w[,1]))/ sigma**3;
   f11 = t(1/sqrt(constant('pi')) * ((s1 - s2) # (mu + sqrt(2) * sigmau * w[,1]) + (s4 - s3) * sigma + s2 * a + (1 - s1) * b))* w[,2]; 
   f12 = t(1/sqrt(constant('pi')) * ((s1 - s2) + (ss1 - ss2) # (mu + sqrt(2) * sigmau * w[,1]) + (ss4 - ss3) * sigma + ss2 * a - ss1 * b))* w[,2]; 
   f13 = t(1/sqrt(constant('pi')) * ((sss1 - sss2) # (mu + sqrt(2) * sigmau * w[,1]) + (s4 - s3) + (sss4 - sss3) * sigma + sss2 * a - sss1 * b))* w[,2]; 
   f14 = t(1/sqrt(constant('pi')) * ((s1 - s2) # (mu + sqrt(2) * sigmau * w[,1]) + (s4 - s3) * sigma + s2 * a + (1 - s1) * b)#((sqrt(2) * sigmau * w[,1])#(sqrt(2) * sigmau * w[,1])/sigmau**3 - 1/sigmau))* w[,2];

   mu = x0[i,]*beta;
   s1 = probnorm((b - (mu + sqrt(2) * sigmau * w[,1]))/sigma);
   s2 = probnorm((a - (mu + sqrt(2) * sigmau * w[,1]))/sigma);
   s3 = pdf ('normal', (b - (mu + sqrt(2) * sigmau * w[,1]))/sigma, 0, 1);
   s4 = pdf ('normal', (a - (mu + sqrt(2) * sigmau * w[,1]))/sigma, 0, 1);
   ss1 = (-1) * s3/sigma;
   ss2 = (-1) * s4/sigma;
   ss3 = s3 # (b - (mu + sqrt(2) * sigmau * w[,1])) / sigma**2;
   ss4 = s4 # (a - (mu + sqrt(2) * sigmau * w[,1])) / sigma**2;
   sss1 = (-1) * (b - (mu + sqrt(2) * sigmau * w[,1])) # s3/sigma**2;
   sss2 = (-1) * (a - (mu + sqrt(2) * sigmau * w[,1])) # s4/sigma**2;
   sss3 = s3 # (b - (mu + sqrt(2) * sigmau * w[,1]))# (b - (mu + sqrt(2) * sigmau * w[,1]))/ sigma**3;
   sss4 = s4 # (a - (mu + sqrt(2) * sigmau * w[,1]))# (a - (mu + sqrt(2) * sigmau * w[,1]))/ sigma**3;
   f01 = t(1/sqrt(constant('pi')) * ((s1 - s2) # (mu + sqrt(2) * sigmau * w[,1]) + (s4 - s3) * sigma + s2 * a + (1 - s1) * b))* w[,2]; 
   f02 = t(1/sqrt(constant('pi')) * ((s1 - s2) + (ss1 - ss2) # (mu + sqrt(2) * sigmau * w[,1]) + (ss4 - ss3) * sigma + ss2 * a - ss1 * b))* w[,2]; 
   f03 = t(1/sqrt(constant('pi')) * ((sss1 - sss2) # (mu + sqrt(2) * sigmau * w[,1]) + (s4 - s3) + (sss4 - sss3) * sigma + sss2 * a - sss1 * b))* w[,2]; 
   f04 = t(1/sqrt(constant('pi')) * ((s1 - s2) # (mu + sqrt(2) * sigmau * w[,1]) + (s4 - s3) * sigma + s2 * a + (1 - s1) * b)#((sqrt(2) * sigmau * w[,1])#(sqrt(2) * sigmau * w[,1])/sigmau**3 - 1/sigmau))* w[,2];

   sum1 = sum1 + (f11 - f01);
   g[1:k1] = g[1:k1] + (f12  * t(x1[i,]) - f02 * t(x0[i,])); 
   g[k] = g[k] + (f13 - f03);
   g[k+1] = g[k+1] + (f14 - f04);
end;

sum1 = sum1/n1;
g = g/n1;
finalvar = g * var * t(g);
if finalvar >= 0 then finalsd = sqrt(finalvar);
print sum1 finalvar finalsd;

postprobs=like||aic||bic||sum1||finalsd||rc||theta||t(sd);
 cname = {"Like" "AIC" "BIC" "EST" "STD" "STA" };

create ff11 from postprobs  [ colname=cname ];
append from postprobs;

quit;
run;
%mend;

/* Two Boundary Direct-Marginalization Approach with Newton-Raphson Method using 40-point Gaussian-Hermite Quadrature with Gradient Specification*/

%macro tobitcorr2(dsn=, y=,x=, cutoff1=, cutoff2 =, theta= 0);
proc iml;
reset noname;
use &dsn;
read all var {&x} into x;
read all var {&y} into y;
yx=y || x;
* a = lower truncation point;
use tt3;
read all into est;
use abswei;
read all var {abs weight} into w;
a = &cutoff1;
b = &cutoff2;
* define var indicating censored ohs;
limit = y;
limit1 = y;
limit2 = y;
limit1 [loc(y=a)] = .;
limit2 [loc(y=b)] = .;
limit3 = limit1;
limit3 [loc(limit1=b)] = .;
n0 = nrow(y);
n1 = ncol(loc( limit1=.));
n2 = ncol(loc( limit2=.));
print n1 "observations are censored at &cutoff1 and" n2 "observations are censored at &cutoff2 out of " n0 " obs";
* define log likelihood function;
*print limit1, limit2, limit3;

start LL(theta) global (yx, limit, limit1, limit2, limit3, a, b, w);
* split data into censored and uncensored ohs;
yx0 = yx[loc((limit3 ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx1 = yx[loc(( limit1 = .)),];
n1 = nrow(yx1);
y1 = yx1[,1];
x1 = j(n1,1,1) || yx1[,2:k];

yx2 = yx[loc(( limit2 = .)),];
n2 = nrow(yx2);
y2 = yx2[,1];
x2 = j(n2,1,1) || yx2[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta) - 1;
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma * sigma;
sigmau = theta [k + 1];
sigmau2 = sigmau * sigmau;
sum1 = 0;
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do i = 1 to n1;
   mu1 = x1[i,]*beta;
   mu = mu1;
   diff = 2;
   do new = 1 to 200 while (abs(diff) > 1E-8);
   temp = mu;

   z1 = probnorm((b - (mu + sqrt(2) * sigmau * w[,1]))/sigma);
   z2 = probnorm((a - (mu + sqrt(2) * sigmau * w[,1]))/sigma);
   z3 = pdf ('normal', (b - (mu + sqrt(2) * sigmau * w[,1]))/sigma, 0, 1);
   z4 = pdf ('normal', (a - (mu + sqrt(2) * sigmau * w[,1]))/sigma, 0, 1);
   zz1 = (-1) * z3/sigma;
   zz2 = (-1) * z4/sigma;
   zz3 = z3 # (b - (mu + sqrt(2) * sigmau * w[,1])) / sigma**2;
   zz4 = z4 # (a - (mu + sqrt(2) * sigmau * w[,1])) / sigma**2;
   f1 = t(1/sqrt(constant('pi')) * ((z1 - z2) # (mu + sqrt(2) * sigmau * w[,1]) + (z4 - z3) * sigma + z2 * a + (1 - z1) * b))* w[,2]; 
   f2 = t(1/sqrt(constant('pi')) * ((z1 - z2) + (zz1 - zz2) # (mu + sqrt(2) * sigmau * w[,1]) + (zz4 - zz3) * sigma + zz2 * a - zz1 * b))* w[,2]; 

   mu = mu - sign(f2) * (f1 - mu1)/( max(abs(f2), exp(-150)));

   diff = mu - temp;
   end;
   chk1 = t(1/sqrt(constant('pi')) * probnorm((a - mu - sqrt(2) * sigmau * w[,1])/sigma))* w[,2];
   if chk1 <= 0 then chk1 = 0.1e8;
   sum1 = sum1 + log(chk1);
end;

do j = 1 to n2;
   mu1 = x2[j,]*beta;
   mu = mu1;
   diff = 2;
   do new = 1 to 200 while (abs(diff) > 1E-8);
   temp = mu;

   z1 = probnorm((b - (mu + sqrt(2) * sigmau * w[,1]))/sigma);
   z2 = probnorm((a - (mu + sqrt(2) * sigmau * w[,1]))/sigma);
   z3 = pdf ('normal', (b - (mu + sqrt(2) * sigmau * w[,1]))/sigma, 0, 1);
   z4 = pdf ('normal', (a - (mu + sqrt(2) * sigmau * w[,1]))/sigma, 0, 1);
   zz1 = (-1) * z3/sigma;
   zz2 = (-1) * z4/sigma;
   zz3 = z3 # (b - (mu + sqrt(2) * sigmau * w[,1])) / sigma**2;
   zz4 = z4 # (a - (mu + sqrt(2) * sigmau * w[,1])) / sigma**2;
   f1 = t(1/sqrt(constant('pi')) * ((z1 - z2) # (mu + sqrt(2) * sigmau * w[,1]) + (z4 - z3) * sigma + z2 * a + (1 - z1) * b))* w[,2]; 
   f2 = t(1/sqrt(constant('pi')) * ((z1 - z2) + (zz1 - zz2) # (mu + sqrt(2) * sigmau * w[,1]) + (zz4 - zz3) * sigma + zz2 * a - zz1 * b))* w[,2]; 

   mu = mu - sign(f2) * (f1 - mu1)/( max(abs(f2), exp(-150)));

   diff = mu - temp;
   end;
   chk2 = t(1/sqrt(constant('pi')) * probnorm((mu + sqrt(2) * sigmau * w[,1] - b)/sigma))* w[,2];
   if chk2 <= 0 then chk2 = 0.1e8;
   sum1 = sum1 + log(chk2);
end;

do e = 1 to n0;
   mu1 = x0[e,]*beta;
   mu = mu1;
   diff = 2;
   do new = 1 to 200 while (abs(diff) > 1E-8);
   temp = mu;

   z1 = probnorm((b - (mu + sqrt(2) * sigmau * w[,1]))/sigma);
   z2 = probnorm((a - (mu + sqrt(2) * sigmau * w[,1]))/sigma);
   z3 = pdf ('normal', (b - (mu + sqrt(2) * sigmau * w[,1]))/sigma, 0, 1);
   z4 = pdf ('normal', (a - (mu + sqrt(2) * sigmau * w[,1]))/sigma, 0, 1);
   zz1 = (-1) * z3/sigma;
   zz2 = (-1) * z4/sigma;
   zz3 = z3 # (b - (mu + sqrt(2) * sigmau * w[,1])) / sigma**2;
   zz4 = z4 # (a - (mu + sqrt(2) * sigmau * w[,1])) / sigma**2;
   f1 = t(1/sqrt(constant('pi')) * ((z1 - z2) # (mu + sqrt(2) * sigmau * w[,1]) + (z4 - z3) * sigma + z2 * a + (1 - z1) * b))* w[,2]; 
   f2 = t(1/sqrt(constant('pi')) * ((z1 - z2) + (zz1 - zz2) # (mu + sqrt(2) * sigmau * w[,1]) + (zz4 - zz3) * sigma + zz2 * a - zz1 * b))* w[,2]; 

   mu = mu - sign(f2) * (f1 - mu1)/( max(abs(f2), exp(-150)));

   diff = mu - temp;
   end;
   yi = y0[e];
   chk3 = t(1/sqrt(constant('pi')) /sqrt(2*constant('pi')*sigma**2) * exp((-1)*(yi - mu - sqrt(2) * sigmau * w[,1])
             #(yi - mu - sqrt(2) * sigmau * w[,1])/(2*sigma**2)))* w[,2];
   if chk3 <= 0 then chk3 = 0.1e8;
   sum1 = sum1 + log(chk3);
end;

f = sum1;
return(f);
finish LL ;

optn={1 2};
*tc={3000 15000 . 1E-12};
tc={3000 15000 . 1E-6};
*tc={3000 15000};
* define starting va1ues for theta: either user-provided or O1S;
theta0 = &theta;
if theta0 = 0 then
do;
theta0 = t(est);
print 'Mixed Effects Model ignoring Boundary starting values for theta = '
theta0;
end;
else print 'User-supplied starting values for
theta = ' theta0;
t1 = j(1,ncol(theta0),.);
t2 = t1;
t2[ncol(theta0)] = 0;
t2[ncol(theta0)-1] = 0;
t3 = t2//t1;
con = t3;
call nlpqn(rc, theta, 'LL',theta0,optn, con, tc);
*call nlpqn(rc, theta, 'LL',theta0,optn, con, tc) grd='GRAD';
call nlpfdd(f, g, h,'LL', theta);
var = inv(-h);
sd = t(j(1,ncol(theta0),.));
temp = vecdiag(var);
do i = 1 to nrow(temp);
if temp[i] >= 0 then sd[i] = sqrt(temp[i]);
else sd[i] = .;
end; 
*sd = sqrt(vecdiag(var));
print 'Hessian = ' h ' Covariance matrix = ' var
' Standard errors = ' sd;
*print 'Hessian = ' h ' Covariance matrix = ' var
;
*print rc theta;

yx0 = yx[loc((limit3 ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx1 = yx[loc(( limit1 = .)),];
n1 = nrow(yx1);
y1 = yx1[,1];
x1 = j(n1,1,1) || yx1[,2:k];

yx2 = yx[loc(( limit2 = .)),];
n2 = nrow(yx2);
y2 = yx2[,1];
x2 = j(n2,1,1) || yx2[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta) - 1;
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma * sigma;
sigmau = theta [k + 1];
sigmau2 = sigmau * sigmau;
sum1 = 0;
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do i = 1 to n1;
   mu1 = x1[i,]*beta;
   mu = mu1;
   diff = 2;
   do new = 1 to 200 while (abs(diff) > 1E-8);
   temp = mu;

   z1 = probnorm((b - (mu + sqrt(2) * sigmau * w[,1]))/sigma);
   z2 = probnorm((a - (mu + sqrt(2) * sigmau * w[,1]))/sigma);
   z3 = pdf ('normal', (b - (mu + sqrt(2) * sigmau * w[,1]))/sigma, 0, 1);
   z4 = pdf ('normal', (a - (mu + sqrt(2) * sigmau * w[,1]))/sigma, 0, 1);
   zz1 = (-1) * z3/sigma;
   zz2 = (-1) * z4/sigma;
   zz3 = z3 # (b - (mu + sqrt(2) * sigmau * w[,1])) / sigma**2;
   zz4 = z4 # (a - (mu + sqrt(2) * sigmau * w[,1])) / sigma**2;
   f1 = t(1/sqrt(constant('pi')) * ((z1 - z2) # (mu + sqrt(2) * sigmau * w[,1]) + (z4 - z3) * sigma + z2 * a + (1 - z1) * b))* w[,2]; 
   f2 = t(1/sqrt(constant('pi')) * ((z1 - z2) + (zz1 - zz2) # (mu + sqrt(2) * sigmau * w[,1]) + (zz4 - zz3) * sigma + zz2 * a - zz1 * b))* w[,2]; 

   mu = mu - sign(f2) * (f1 - mu1)/( max(abs(f2), exp(-150)));

   diff = mu - temp;
   end;
   chk1 = t(1/sqrt(constant('pi')) * probnorm((a - mu - sqrt(2) * sigmau * w[,1])/sigma))* w[,2];
   if chk1 <= 0 then chk1 = 0.1e8;
   sum1 = sum1 + log(chk1);
end;

do j = 1 to n2;
   mu1 = x2[j,]*beta;
   mu = mu1;
   diff = 2;
   do new = 1 to 200 while (abs(diff) > 1E-8);
   temp = mu;

   z1 = probnorm((b - (mu + sqrt(2) * sigmau * w[,1]))/sigma);
   z2 = probnorm((a - (mu + sqrt(2) * sigmau * w[,1]))/sigma);
   z3 = pdf ('normal', (b - (mu + sqrt(2) * sigmau * w[,1]))/sigma, 0, 1);
   z4 = pdf ('normal', (a - (mu + sqrt(2) * sigmau * w[,1]))/sigma, 0, 1);
   zz1 = (-1) * z3/sigma;
   zz2 = (-1) * z4/sigma;
   zz3 = z3 # (b - (mu + sqrt(2) * sigmau * w[,1])) / sigma**2;
   zz4 = z4 # (a - (mu + sqrt(2) * sigmau * w[,1])) / sigma**2;
   f1 = t(1/sqrt(constant('pi')) * ((z1 - z2) # (mu + sqrt(2) * sigmau * w[,1]) + (z4 - z3) * sigma + z2 * a + (1 - z1) * b))* w[,2]; 
   f2 = t(1/sqrt(constant('pi')) * ((z1 - z2) + (zz1 - zz2) # (mu + sqrt(2) * sigmau * w[,1]) + (zz4 - zz3) * sigma + zz2 * a - zz1 * b))* w[,2]; 

   mu = mu - sign(f2) * (f1 - mu1)/( max(abs(f2), exp(-150)));

   diff = mu - temp;
   end;
   chk2 = t(1/sqrt(constant('pi')) * probnorm((mu + sqrt(2) * sigmau * w[,1] - b)/sigma))* w[,2];
   if chk2 <= 0 then chk2 = 0.1e8;
   sum1 = sum1 + log(chk2);
end;

do e = 1 to n0;
   mu1 = x0[e,]*beta;
   mu = mu1;
   diff = 2;
   do new = 1 to 200 while (abs(diff) > 1E-8);
   temp = mu;

   z1 = probnorm((b - (mu + sqrt(2) * sigmau * w[,1]))/sigma);
   z2 = probnorm((a - (mu + sqrt(2) * sigmau * w[,1]))/sigma);
   z3 = pdf ('normal', (b - (mu + sqrt(2) * sigmau * w[,1]))/sigma, 0, 1);
   z4 = pdf ('normal', (a - (mu + sqrt(2) * sigmau * w[,1]))/sigma, 0, 1);
   zz1 = (-1) * z3/sigma;
   zz2 = (-1) * z4/sigma;
   zz3 = z3 # (b - (mu + sqrt(2) * sigmau * w[,1])) / sigma**2;
   zz4 = z4 # (a - (mu + sqrt(2) * sigmau * w[,1])) / sigma**2;
   f1 = t(1/sqrt(constant('pi')) * ((z1 - z2) # (mu + sqrt(2) * sigmau * w[,1]) + (z4 - z3) * sigma + z2 * a + (1 - z1) * b))* w[,2]; 
   f2 = t(1/sqrt(constant('pi')) * ((z1 - z2) + (zz1 - zz2) # (mu + sqrt(2) * sigmau * w[,1]) + (zz4 - zz3) * sigma + zz2 * a - zz1 * b))* w[,2]; 

   mu = mu - sign(f2) * (f1 - mu1)/( max(abs(f2), exp(-150)));

   diff = mu - temp;
   end;
   yi = y0[e];
   chk3 = t(1/sqrt(constant('pi')) /sqrt(2*constant('pi')*sigma**2) * exp((-1)*(yi - mu - sqrt(2) * sigmau * w[,1])
             #(yi - mu - sqrt(2) * sigmau * w[,1])/(2*sigma**2)))* w[,2];
   if chk3 <= 0 then chk3 = 0.1e8;
   sum1 = sum1 + log(chk3);
end;

like = sum1*(-2);
aic = sum1*(-2) + 2 * ncol(theta);
bic = ncol(theta) * log(nrow(yx0)) + sum1*(-2);

sum1 = theta[k-1];
finalsd = sd[k-1];

postprobs=like||aic||bic||sum1||finalsd||rc||theta||t(sd);
 cname = {"Like" "AIC" "BIC" "EST" "STD" "STA" };

create ff12 from postprobs  [ colname=cname ];
append from postprobs;

quit;
run;
%mend;

/* One Upper Boundary APV Method using 40-point Gaussian-Hermite Quadrature, T = 1 as reference*/

%macro tobitcorr3(dsn=, y=,x=, cutoff2 = , theta= 0);
proc iml;
reset noname;
use &dsn;
read all var {&x} into x;
read all var {&y} into y;
yx=y || x;
use tt3;
read all into est;
use abswei;
read all var {abs weight} into w;
* a = lower truncation point;
b = &cutoff2;
* define var indicating censored ohs;
limit = y;
limit [loc(y=b)] = .;
n0 = nrow(y);
n2 = ncol(loc( limit=.));
print n2 "observations are censored at &cutoff2 out of " n0 " obs";
* define log likelihood function;
*print limit1, limit2, limit3;
start LL(theta) global (yx, limit, b, w);
* split data into censored and uncensored ohs;
yx0 = yx[loc((limit ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx2 = yx[loc(( limit = .)),];
n2 = nrow(yx2);
y2 = yx2[,1];
x2 = j(n2,1,1) || yx2[,2:k];

k = ncol(theta) - 1;
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma * sigma;
sigmau = theta [k + 1];
sigmau2 = sigmau * sigmau;
sum1 = 0;
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do j = 1 to n2;
   mu = x2[j,]*beta;
   chk2 = t(1/sqrt(constant('pi')) * probnorm((mu + sqrt(2) * sigmau * w[,1] - b)/sigma))* w[,2];
   if chk2 <= 0 then chk2 = 0.1e8;
   sum1 = sum1 + log(chk2);
end;

do e = 1 to n0;
   mu = x0[e,]*beta;
   yi = y0[e];
   chk3 = t(1/sqrt(constant('pi')) /sqrt(2*constant('pi')*sigma**2) * exp((-1)*(yi - mu - sqrt(2) * sigmau * w[,1])
             #(yi - mu - sqrt(2) * sigmau * w[,1])/(2*sigma**2)))* w[,2];
   if chk3 <= 0 then chk3 = 0.1e8;
   sum1 = sum1 + log(chk3);
end;

f = sum1;
return(f);
finish LL ;

optn={1 2};
*tc={3000 15000 . 1E-12};
tc={3000 15000 . 1E-6};
*tc={3000 15000};
* define starting va1ues for theta: either user-provided or O1S;
theta0 = &theta;
if theta0 = 0 then
do;
theta0 = t(est);
print 'Mixed Effects Model ignoring Boundary starting values for theta = '
theta0;
end;
else print 'User-supplied starting values for
theta = ' theta0;
t1 = j(1,ncol(theta0),.);
t2 = t1;
t2[ncol(theta0)] = 0;
t2[ncol(theta0)-1] = 0;
t3 = t2//t1;
con = t3;
call nlpqn(rc, theta, 'LL',theta0,optn, con, tc);
*call nlpqn(rc, theta, 'LL',theta0,optn, con, tc) grd='GRAD';
call nlpfdd(f, g, h,'LL', theta);
var = inv(-h);
sd = t(j(1,ncol(theta0),.));
temp = vecdiag(var);
do i = 1 to nrow(temp);
if temp[i] >= 0 then sd[i] = sqrt(temp[i]);
else sd[i] = .;
end; 
*sd = sqrt(vecdiag(var));
print 'Hessian = ' h ' Covariance matrix = ' var
' Standard errors = ' sd;
*print 'Hessian = ' h ' Covariance matrix = ' var
;
*print rc theta;

yx0 = yx[loc((limit ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx2 = yx[loc(( limit = .)),];
n2 = nrow(yx2);
y2 = yx2[,1];
x2 = j(n2,1,1) || yx2[,2:k];

k = ncol(theta) - 1;
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma * sigma;
sigmau = theta [k + 1];
sigmau2 = sigmau * sigmau;
sum1 = 0;
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do j = 1 to n2;
   mu = x2[j,]*beta;
   chk2 = t(1/sqrt(constant('pi')) * probnorm((mu + sqrt(2) * sigmau * w[,1] - b)/sigma))* w[,2];
   if chk2 <= 0 then chk2 = 0.1e8;
   sum1 = sum1 + log(chk2);
end;

do e = 1 to n0;
   mu = x0[e,]*beta;
   yi = y0[e];
   chk3 = t(1/sqrt(constant('pi')) /sqrt(2*constant('pi')*sigma**2) * exp((-1)*(yi - mu - sqrt(2) * sigmau * w[,1])
             #(yi - mu - sqrt(2) * sigmau * w[,1])/(2*sigma**2)))* w[,2];
   if chk3 <= 0 then chk3 = 0.1e8;
   sum1 = sum1 + log(chk3);
end;

like = sum1*(-2);
aic = sum1*(-2) + 2 * ncol(theta);
bic = ncol(theta) * log(nrow(yx0)) + sum1*(-2);

k = ncol(theta) - 1;
k1 = k-1;
beta = theta [1:k1];
x1 = x[loc((x[,k-2] = 1)),];
n1 = nrow(x1);
if k >= 4 then x0 = j(n1,1,1) || x1[,1:(k-3)]||j(n1,1,0);
else if k = 3 then x0 = j(n1,1,1) ||j(n1,1,0);
x1 = j(n1,1,1) || x1;
sigma = theta [k];
sigma2 = sigma * sigma;
sigmau = theta [k + 1];
sigmau2 = sigmau * sigmau;
sum1 = 0;
g = j(1,k+1,0);

do i = 1 to n1;
   mu = x1[i,]*beta;
   s1 = probnorm((b - (mu + sqrt(2) * sigmau * w[,1]))/sigma);
   s2 = 0;
   s3 = pdf ('normal', (b - (mu + sqrt(2) * sigmau * w[,1]))/sigma, 0, 1);
   s4 = 0;
   ss1 = (-1) * s3/sigma;
   ss2 = 0;
   ss3 = s3 # (b - (mu + sqrt(2) * sigmau * w[,1])) / sigma**2;
   ss4 = 0;
   sss1 = (-1) * (b - (mu + sqrt(2) * sigmau * w[,1])) # s3/sigma**2;
   sss2 = 0;
   sss3 = s3 # (b - (mu + sqrt(2) * sigmau * w[,1]))# (b - (mu + sqrt(2) * sigmau * w[,1]))/ sigma**3;
   sss4 = 0;
   f11 = t(1/sqrt(constant('pi')) * ((s1 - s2) # (mu + sqrt(2) * sigmau * w[,1]) + (s4 - s3) * sigma + (1 - s1) * b))* w[,2]; 
   f12 = t(1/sqrt(constant('pi')) * ((s1 - s2) + (ss1 - ss2) # (mu + sqrt(2) * sigmau * w[,1]) + (ss4 - ss3) * sigma - ss1 * b))* w[,2]; 
   f13 = t(1/sqrt(constant('pi')) * ((sss1 - sss2) # (mu + sqrt(2) * sigmau * w[,1]) + (s4 - s3) + (sss4 - sss3) * sigma - sss1 * b))* w[,2]; 
   f14 = t(1/sqrt(constant('pi')) * ((s1 - s2) # (mu + sqrt(2) * sigmau * w[,1]) + (s4 - s3) * sigma + (1 - s1) * b)#((sqrt(2) * sigmau * w[,1])#(sqrt(2) * sigmau * w[,1])/sigmau**3 - 1/sigmau))* w[,2];

   mu = x0[i,]*beta;
   s1 = probnorm((b - (mu + sqrt(2) * sigmau * w[,1]))/sigma);
   s2 = 0;
   s3 = pdf ('normal', (b - (mu + sqrt(2) * sigmau * w[,1]))/sigma, 0, 1);
   s4 = 0;
   ss1 = (-1) * s3/sigma;
   ss2 = 0;
   ss3 = s3 # (b - (mu + sqrt(2) * sigmau * w[,1])) / sigma**2;
   ss4 = 0;
   sss1 = (-1) * (b - (mu + sqrt(2) * sigmau * w[,1])) # s3/sigma**2;
   sss2 = 0;
   sss3 = s3 # (b - (mu + sqrt(2) * sigmau * w[,1]))# (b - (mu + sqrt(2) * sigmau * w[,1]))/ sigma**3;
   sss4 = 0;
   f01 = t(1/sqrt(constant('pi')) * ((s1 - s2) # (mu + sqrt(2) * sigmau * w[,1]) + (s4 - s3) * sigma + (1 - s1) * b))* w[,2]; 
   f02 = t(1/sqrt(constant('pi')) * ((s1 - s2) + (ss1 - ss2) # (mu + sqrt(2) * sigmau * w[,1]) + (ss4 - ss3) * sigma - ss1 * b))* w[,2]; 
   f03 = t(1/sqrt(constant('pi')) * ((sss1 - sss2) # (mu + sqrt(2) * sigmau * w[,1]) + (s4 - s3) + (sss4 - sss3) * sigma - sss1 * b))* w[,2]; 
   f04 = t(1/sqrt(constant('pi')) * ((s1 - s2) # (mu + sqrt(2) * sigmau * w[,1]) + (s4 - s3) * sigma + (1 - s1) * b)#((sqrt(2) * sigmau * w[,1])#(sqrt(2) * sigmau * w[,1])/sigmau**3 - 1/sigmau))* w[,2];

   sum1 = sum1 + (f11 - f01);
   g[1:k1] = g[1:k1] + (f12  * t(x1[i,]) - f02 * t(x0[i,])); 
   g[k] = g[k] + (f13 - f03);
   g[k+1] = g[k+1] + (f14 - f04);
end;

sum1 = sum1/n1;
g = g/n1;
finalvar = g * var * t(g);
if finalvar >= 0 then finalsd = sqrt(finalvar);
print sum1 finalvar finalsd;

postprobs=like||aic||bic||sum1||finalsd||rc||theta||t(sd);
 cname = {"Like" "AIC" "BIC" "EST" "STD" "STA" };

create ff11 from postprobs  [ colname=cname ];
append from postprobs;

quit;
run;
%mend;

/* One Upper Boundary Direct-Marginalization Approach with Newton-Raphson Method using 40-point Gaussian-Hermite Quadrature*/

%macro tobitcorr4(dsn=, y=,x=, cutoff2 =, theta= 0);
proc iml;
reset noname;
use &dsn;
read all var {&x} into x;
read all var {&y} into y;
yx=y || x;
* a = lower truncation point;
use tt3;
read all into est;
use abswei;
read all var {abs weight} into w;
* a = lower truncation point;
b = &cutoff2;
* define var indicating censored ohs;
limit = y;
limit [loc(y=b)] = .;
n0 = nrow(y);
n2 = ncol(loc( limit=.));
print n2 "observations are censored at &cutoff2 out of " n0 " obs";
* define log likelihood function;
*print limit1, limit2, limit3;
start LL(theta) global (yx, limit, b, w);
* split data into censored and uncensored ohs;
yx0 = yx[loc((limit ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx2 = yx[loc(( limit = .)),];
n2 = nrow(yx2);
y2 = yx2[,1];
x2 = j(n2,1,1) || yx2[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta) - 1;
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma * sigma;
sigmau = theta [k + 1];
sigmau2 = sigmau * sigmau;
sum1 = 0;
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do j = 1 to n2;
   mu1 = x2[j,]*beta;
   mu = mu1;
   diff = 2;
   do new = 1 to 200 while (abs(diff) > 1E-8);
   temp = mu;

   z1 = probnorm((b - (mu + sqrt(2) * sigmau * w[,1]))/sigma);
   z2 = 0;
   z3 = pdf ('normal', (b - (mu + sqrt(2) * sigmau * w[,1]))/sigma, 0, 1);
   z4 = 0;
   zz1 = (-1) * z3/sigma;
   zz2 = 0;
   zz3 = z3 # (b - (mu + sqrt(2) * sigmau * w[,1])) / sigma**2;
   zz4 = 0;
   f1 = t(1/sqrt(constant('pi')) * ((z1 - z2) # (mu + sqrt(2) * sigmau * w[,1]) + (z4 - z3) * sigma + (1 - z1) * b))* w[,2]; 
   f2 = t(1/sqrt(constant('pi')) * ((z1 - z2) + (zz1 - zz2) # (mu + sqrt(2) * sigmau * w[,1]) + (zz4 - zz3) * sigma - zz1 * b))* w[,2]; 

   mu = mu - sign(f2) * (f1 - mu1)/( max(abs(f2), exp(-150)));

   diff = mu - temp;
   end;
   chk2 = t(1/sqrt(constant('pi')) * probnorm((mu + sqrt(2) * sigmau * w[,1] - b)/sigma))* w[,2];
   if chk2 <= 0 then chk2 = 0.1e8;
   sum1 = sum1 + log(chk2);
end;

do e = 1 to n0;
   mu1 = x0[e,]*beta;
   mu = mu1;
   diff = 2;
   do new = 1 to 200 while (abs(diff) > 1E-8);
   temp = mu;

   z1 = probnorm((b - (mu + sqrt(2) * sigmau * w[,1]))/sigma);
   z2 = 0;
   z3 = pdf ('normal', (b - (mu + sqrt(2) * sigmau * w[,1]))/sigma, 0, 1);
   z4 = 0;
   zz1 = (-1) * z3/sigma;
   zz2 = 0;
   zz3 = z3 # (b - (mu + sqrt(2) * sigmau * w[,1])) / sigma**2;
   zz4 = 0;
   f1 = t(1/sqrt(constant('pi')) * ((z1 - z2) # (mu + sqrt(2) * sigmau * w[,1]) + (z4 - z3) * sigma + (1 - z1) * b))* w[,2]; 
   f2 = t(1/sqrt(constant('pi')) * ((z1 - z2) + (zz1 - zz2) # (mu + sqrt(2) * sigmau * w[,1]) + (zz4 - zz3) * sigma - zz1 * b))* w[,2]; 

   mu = mu - sign(f2) * (f1 - mu1)/( max(abs(f2), exp(-150)));

   diff = mu - temp;
   end;
   yi = y0[e];
   chk3 = t(1/sqrt(constant('pi')) /sqrt(2*constant('pi')*sigma**2) * exp((-1)*(yi - mu - sqrt(2) * sigmau * w[,1])
             #(yi - mu - sqrt(2) * sigmau * w[,1])/(2*sigma**2)))* w[,2];
   if chk3 <= 0 then chk3 = 0.1e8;
   sum1 = sum1 + log(chk3);
end;

f = sum1;
return(f);
finish LL ;

optn={1 2};
*tc={3000 15000 . 1E-12};
tc={3000 15000 . 1E-6};
*tc={3000 15000};
* define starting va1ues for theta: either user-provided or O1S;
theta0 = &theta;
if theta0 = 0 then
do;
theta0 = t(est);
print 'Mixed Effects Model ignoring Boundary starting values for theta = '
theta0;
end;
else print 'User-supplied starting values for
theta = ' theta0;
t1 = j(1,ncol(theta0),.);
t2 = t1;
t2[ncol(theta0)] = 0;
t2[ncol(theta0)-1] = 0;
t3 = t2//t1;
con = t3;
call nlpqn(rc, theta, 'LL',theta0,optn, con, tc);
*call nlpqn(rc, theta, 'LL',theta0,optn, con, tc) grd='GRAD';
call nlpfdd(f, g, h,'LL', theta);
var = inv(-h);
sd = t(j(1,ncol(theta0),.));
temp = vecdiag(var);
do i = 1 to nrow(temp);
if temp[i] >= 0 then sd[i] = sqrt(temp[i]);
else sd[i] = .;
end; 
*sd = sqrt(vecdiag(var));
print 'Hessian = ' h ' Covariance matrix = ' var
' Standard errors = ' sd;
*print 'Hessian = ' h ' Covariance matrix = ' var
;
*print rc theta;

yx0 = yx[loc((limit ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx2 = yx[loc(( limit = .)),];
n2 = nrow(yx2);
y2 = yx2[,1];
x2 = j(n2,1,1) || yx2[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta) - 1;
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma * sigma;
sigmau = theta [k + 1];
sigmau2 = sigmau * sigmau;
sum1 = 0;
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do j = 1 to n2;
   mu1 = x2[j,]*beta;
   mu = mu1;
   diff = 2;
   do new = 1 to 200 while (abs(diff) > 1E-8);
   temp = mu;

   z1 = probnorm((b - (mu + sqrt(2) * sigmau * w[,1]))/sigma);
   z2 = 0;
   z3 = pdf ('normal', (b - (mu + sqrt(2) * sigmau * w[,1]))/sigma, 0, 1);
   z4 = 0;
   zz1 = (-1) * z3/sigma;
   zz2 = 0;
   zz3 = z3 # (b - (mu + sqrt(2) * sigmau * w[,1])) / sigma**2;
   zz4 = 0;
   f1 = t(1/sqrt(constant('pi')) * ((z1 - z2) # (mu + sqrt(2) * sigmau * w[,1]) + (z4 - z3) * sigma + (1 - z1) * b))* w[,2]; 
   f2 = t(1/sqrt(constant('pi')) * ((z1 - z2) + (zz1 - zz2) # (mu + sqrt(2) * sigmau * w[,1]) + (zz4 - zz3) * sigma - zz1 * b))* w[,2]; 

   mu = mu - sign(f2) * (f1 - mu1)/( max(abs(f2), exp(-150)));

   diff = mu - temp;
   end;
   chk2 = t(1/sqrt(constant('pi')) * probnorm((mu + sqrt(2) * sigmau * w[,1] - b)/sigma))* w[,2];
   if chk2 <= 0 then chk2 = 0.1e8;
   sum1 = sum1 + log(chk2);
end;

do e = 1 to n0;
   mu1 = x0[e,]*beta;
   mu = mu1;
   diff = 2;
   do new = 1 to 200 while (abs(diff) > 1E-8);
   temp = mu;

   z1 = probnorm((b - (mu + sqrt(2) * sigmau * w[,1]))/sigma);
   z2 = 0;
   z3 = pdf ('normal', (b - (mu + sqrt(2) * sigmau * w[,1]))/sigma, 0, 1);
   z4 = 0;
   zz1 = (-1) * z3/sigma;
   zz2 = 0;
   zz3 = z3 # (b - (mu + sqrt(2) * sigmau * w[,1])) / sigma**2;
   zz4 = 0;
   f1 = t(1/sqrt(constant('pi')) * ((z1 - z2) # (mu + sqrt(2) * sigmau * w[,1]) + (z4 - z3) * sigma + (1 - z1) * b))* w[,2]; 
   f2 = t(1/sqrt(constant('pi')) * ((z1 - z2) + (zz1 - zz2) # (mu + sqrt(2) * sigmau * w[,1]) + (zz4 - zz3) * sigma - zz1 * b))* w[,2]; 

   mu = mu - sign(f2) * (f1 - mu1)/( max(abs(f2), exp(-150)));

   diff = mu - temp;
   end;
   yi = y0[e];
   chk3 = t(1/sqrt(constant('pi')) /sqrt(2*constant('pi')*sigma**2) * exp((-1)*(yi - mu - sqrt(2) * sigmau * w[,1])
             #(yi - mu - sqrt(2) * sigmau * w[,1])/(2*sigma**2)))* w[,2];
   if chk3 <= 0 then chk3 = 0.1e8;
   sum1 = sum1 + log(chk3);
end;

like = sum1*(-2);
aic = sum1*(-2) + 2 * ncol(theta);
bic = ncol(theta) * log(nrow(yx0)) + sum1*(-2);

sum1 = theta[k-1];
finalsd = sd[k-1];
postprobs=like||aic||bic||sum1||finalsd||rc||theta||t(sd);
 cname = {"Like" "AIC" "BIC" "EST" "STD" "STA" };

create ff12 from postprobs  [ colname=cname ];
append from postprobs;

quit;
run;
%mend;

/* One Lower Boundary APV Method using 40-point Gaussian-Hermite Quadrature, T = 1 as reference*/

%macro tobitcorr5(dsn=, y=,x=, cutoff1=, theta= 0);
proc iml;
reset noname;
use &dsn;
read all var {&x} into x;
read all var {&y} into y;
yx=y || x;
* a = lower truncation point;
use tt3;
read all into est;
use abswei;
read all var {abs weight} into w;
a = &cutoff1;
* define var indicating censored ohs;
limit = y;
limit [loc(y=a)] = .;
n0 = nrow(y);
n1 = ncol(loc( limit=.));
print n1 "observations are censored at &cutoff1 out of" n0 " obs";
* define log likelihood function;
*print limit1, limit2, limit3;
start LL(theta) global (yx, limit, a, w);
* split data into censored and uncensored ohs;
yx0 = yx[loc((limit ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx1 = yx[loc(( limit = .)),];
n1 = nrow(yx1);
y1 = yx1[,1];
x1 = j(n1,1,1) || yx1[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta) - 1;
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma * sigma;
sigmau = theta [k + 1];
sigmau2 = sigmau * sigmau;
sum1 = 0;
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do i = 1 to n1;
   mu = x1[i,]*beta;
   chk1 = t(1/sqrt(constant('pi')) * probnorm((a - mu - sqrt(2) * sigmau * w[,1])/sigma))* w[,2];
   if chk1 <= 0 then chk1 = 0.1e8;
   sum1 = sum1 + log(chk1);
end;

do e = 1 to n0;
   mu = x0[e,]*beta;
   yi = y0[e];
   chk3 = t(1/sqrt(constant('pi')) /sqrt(2*constant('pi')*sigma**2) * exp((-1)*(yi - mu - sqrt(2) * sigmau * w[,1])
             #(yi - mu - sqrt(2) * sigmau * w[,1])/(2*sigma**2)))* w[,2];
   if chk3 <= 0 then chk3 = 0.1e8;
   sum1 = sum1 + log(chk3);
end;

f = sum1;
return(f);
finish LL ;

optn={1 2};
*tc={3000 15000 . 1E-12};
tc={3000 15000 . 1E-6};
*tc={3000 15000};
* define starting va1ues for theta: either user-provided or O1S;
theta0 = &theta;
if theta0 = 0 then
do;
theta0 = t(est);
print 'Mixed Effects Model ignoring Boundary starting values for theta = '
theta0;
end;
else print 'User-supplied starting values for
theta = ' theta0;
t1 = j(1,ncol(theta0),.);
t2 = t1;
t2[ncol(theta0)] = 0;
t2[ncol(theta0)-1] = 0;
t3 = t2//t1;
con = t3;
call nlpqn(rc, theta, 'LL',theta0,optn, con, tc);
*call nlpqn(rc, theta, 'LL',theta0,optn, con, tc) grd='GRAD';
call nlpfdd(f, g, h,'LL', theta);
var = inv(-h);
sd = t(j(1,ncol(theta0),.));
temp = vecdiag(var);
do i = 1 to nrow(temp);
if temp[i] >= 0 then sd[i] = sqrt(temp[i]);
else sd[i] = .;
end; 
*sd = sqrt(vecdiag(var));
print 'Hessian = ' h ' Covariance matrix = ' var
' Standard errors = ' sd;
*print 'Hessian = ' h ' Covariance matrix = ' var
;
*print rc theta;

yx0 = yx[loc((limit ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx1 = yx[loc(( limit = .)),];
n1 = nrow(yx1);
y1 = yx1[,1];
x1 = j(n1,1,1) || yx1[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta) - 1;
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma * sigma;
sigmau = theta [k + 1];
sigmau2 = sigmau * sigmau;
sum1 = 0;
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do i = 1 to n1;
   mu = x1[i,]*beta;
   chk1 = t(1/sqrt(constant('pi')) * probnorm((a - mu - sqrt(2) * sigmau * w[,1])/sigma))* w[,2];
   if chk1 <= 0 then chk1 = 0.1e8;
   sum1 = sum1 + log(chk1);
end;

do e = 1 to n0;
   mu = x0[e,]*beta;
   yi = y0[e];
   chk3 = t(1/sqrt(constant('pi')) /sqrt(2*constant('pi')*sigma**2) * exp((-1)*(yi - mu - sqrt(2) * sigmau * w[,1])
             #(yi - mu - sqrt(2) * sigmau * w[,1])/(2*sigma**2)))* w[,2];
   if chk3 <= 0 then chk3 = 0.1e8;
   sum1 = sum1 + log(chk3);
end;

like = sum1*(-2);
aic = sum1*(-2) + 2 * ncol(theta);
bic = ncol(theta) * log(nrow(yx0)) + sum1*(-2);

k = ncol(theta) - 1;
k1 = k-1;
beta = theta [1:k1];
x1 = x[loc((x[,k-2] = 1)),];
n1 = nrow(x1);
if k >= 4 then x0 = j(n1,1,1) || x1[,1:(k-3)]||j(n1,1,0);
else if k = 3 then x0 = j(n1,1,1) ||j(n1,1,0);
x1 = j(n1,1,1) || x1;
sigma = theta [k];
sigma2 = sigma * sigma;
sigmau = theta [k + 1];
sigmau2 = sigmau * sigmau;
sum1 = 0;
g = j(1,k+1,0);

do i = 1 to n1;
   mu = x1[i,]*beta;
   s1 = 1;
   s2 = probnorm((a - (mu + sqrt(2) * sigmau * w[,1]))/sigma);
   s3 = 0;
   s4 = pdf ('normal', (a - (mu + sqrt(2) * sigmau * w[,1]))/sigma, 0, 1);
   ss1 = 0;
   ss2 = (-1) * s4/sigma;
   ss3 = 0;
   ss4 = s4 # (a - (mu + sqrt(2) * sigmau * w[,1])) / sigma**2;
   sss1 = 0;
   sss2 = (-1) * (a - (mu + sqrt(2) * sigmau * w[,1])) # s4/sigma**2;
   sss3 = 0;
   sss4 = s4 # (a - (mu + sqrt(2) * sigmau * w[,1]))# (a - (mu + sqrt(2) * sigmau * w[,1]))/ sigma**3;
   f11 = t(1/sqrt(constant('pi')) * ((s1 - s2) # (mu + sqrt(2) * sigmau * w[,1]) + (s4 - s3) * sigma + s2 * a))* w[,2]; 
   f12 = t(1/sqrt(constant('pi')) * ((s1 - s2) + (ss1 - ss2) # (mu + sqrt(2) * sigmau * w[,1]) + (ss4 - ss3) * sigma + ss2 * a))* w[,2]; 
   f13 = t(1/sqrt(constant('pi')) * ((sss1 - sss2) # (mu + sqrt(2) * sigmau * w[,1]) + (s4 - s3) + (sss4 - sss3) * sigma + sss2 * a))* w[,2]; 
   f14 = t(1/sqrt(constant('pi')) * ((s1 - s2) # (mu + sqrt(2) * sigmau * w[,1]) + (s4 - s3) * sigma + s2 * a)#((sqrt(2) * sigmau * w[,1])#(sqrt(2) * sigmau * w[,1])/sigmau**3 - 1/sigmau))* w[,2];

   mu = x0[i,]*beta;
   s1 = 1;
   s2 = probnorm((a - (mu + sqrt(2) * sigmau * w[,1]))/sigma);
   s3 = 0;
   s4 = pdf ('normal', (a - (mu + sqrt(2) * sigmau * w[,1]))/sigma, 0, 1);
   ss1 = 0;
   ss2 = (-1) * s4/sigma;
   ss3 = 0;
   ss4 = s4 # (a - (mu + sqrt(2) * sigmau * w[,1])) / sigma**2;
   sss1 = 0;
   sss2 = (-1) * (a - (mu + sqrt(2) * sigmau * w[,1])) # s4/sigma**2;
   sss3 = 0;
   sss4 = s4 # (a - (mu + sqrt(2) * sigmau * w[,1]))# (a - (mu + sqrt(2) * sigmau * w[,1]))/ sigma**3;
   f01 = t(1/sqrt(constant('pi')) * ((s1 - s2) # (mu + sqrt(2) * sigmau * w[,1]) + (s4 - s3) * sigma + s2 * a))* w[,2]; 
   f02 = t(1/sqrt(constant('pi')) * ((s1 - s2) + (ss1 - ss2) # (mu + sqrt(2) * sigmau * w[,1]) + (ss4 - ss3) * sigma + ss2 * a))* w[,2]; 
   f03 = t(1/sqrt(constant('pi')) * ((sss1 - sss2) # (mu + sqrt(2) * sigmau * w[,1]) + (s4 - s3) + (sss4 - sss3) * sigma + sss2 * a))* w[,2]; 
   f04 = t(1/sqrt(constant('pi')) * ((s1 - s2) # (mu + sqrt(2) * sigmau * w[,1]) + (s4 - s3) * sigma + s2 * a)#((sqrt(2) * sigmau * w[,1])#(sqrt(2) * sigmau * w[,1])/sigmau**3 - 1/sigmau))* w[,2];

   sum1 = sum1 + (f11 - f01);
   g[1:k1] = g[1:k1] + (f12  * t(x1[i,]) - f02 * t(x0[i,])); 
   g[k] = g[k] + (f13 - f03);
   g[k+1] = g[k+1] + (f14 - f04);
end;

sum1 = sum1/n1;
g = g/n1;
finalvar = g * var * t(g);
if finalvar >= 0 then finalsd = sqrt(finalvar);
print sum1 finalvar finalsd;

postprobs=like||aic||bic||sum1||finalsd||rc||theta||t(sd);
 cname = {"Like" "AIC" "BIC" "EST" "STD" "STA" };

create ff11 from postprobs  [ colname=cname ];
append from postprobs;

quit;
run;
%mend;

/* One Lower Boundary Direct-Marginalization Approach with Newton-Raphson Method using 40-point Gaussian-Hermite Quadrature*/

%macro tobitcorr6(dsn=, y=,x=, cutoff1=, theta= 0);
proc iml;
reset noname;
use &dsn;
read all var {&x} into x;
read all var {&y} into y;
yx=y || x;
* a = lower truncation point;
use tt3;
read all into est;
use abswei;
read all var {abs weight} into w;
a = &cutoff1;
* define var indicating censored ohs;
limit = y;
limit [loc(y=a)] = .;
n0 = nrow(y);
n1 = ncol(loc( limit=.));
print n1 "observations are censored at &cutoff1 out of" n0 " obs";
* define log likelihood function;
*print limit1, limit2, limit3;
start LL(theta) global (yx, limit, a, w);
* split data into censored and uncensored ohs;
yx0 = yx[loc((limit ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx1 = yx[loc(( limit = .)),];
n1 = nrow(yx1);
y1 = yx1[,1];
x1 = j(n1,1,1) || yx1[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta) - 1;
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma * sigma;
sigmau = theta [k + 1];
sigmau2 = sigmau * sigmau;
sum1 = 0;
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do i = 1 to n1;
   mu1 = x1[i,]*beta;
   mu = mu1;
   diff = 2;
   do new = 1 to 200 while (abs(diff) > 1E-8);
   temp = mu;

   z1 = 1;
   z2 = probnorm((a - (mu + sqrt(2) * sigmau * w[,1]))/sigma);
   z3 = 0;
   z4 = pdf ('normal', (a - (mu + sqrt(2) * sigmau * w[,1]))/sigma, 0, 1);
   zz1 = 0;
   zz2 = (-1) * z4/sigma;
   zz3 = 0;
   zz4 = z4 # (a - (mu + sqrt(2) * sigmau * w[,1])) / sigma**2;
   f1 = t(1/sqrt(constant('pi')) * ((z1 - z2) # (mu + sqrt(2) * sigmau * w[,1]) + (z4 - z3) * sigma + z2 * a))* w[,2]; 
   f2 = t(1/sqrt(constant('pi')) * ((z1 - z2) + (zz1 - zz2) # (mu + sqrt(2) * sigmau * w[,1]) + (zz4 - zz3) * sigma + zz2 * a))* w[,2]; 

   mu = mu - sign(f2) * (f1 - mu1)/( max(abs(f2), exp(-150)));

   diff = mu - temp;
   end;
   chk1 = t(1/sqrt(constant('pi')) * probnorm((a - mu - sqrt(2) * sigmau * w[,1])/sigma))* w[,2];
   if chk1 <= 0 then chk1 = 0.1e8;
   sum1 = sum1 + log(chk1);
end;

do e = 1 to n0;
   mu1 = x0[e,]*beta;
   mu = mu1;
   diff = 2;
   do new = 1 to 200 while (abs(diff) > 1E-8);
   temp = mu;

   z1 = 1;
   z2 = probnorm((a - (mu + sqrt(2) * sigmau * w[,1]))/sigma);
   z3 = 0;
   z4 = pdf ('normal', (a - (mu + sqrt(2) * sigmau * w[,1]))/sigma, 0, 1);
   zz1 = 0;
   zz2 = (-1) * z4/sigma;
   zz3 = 0;
   zz4 = z4 # (a - (mu + sqrt(2) * sigmau * w[,1])) / sigma**2;
   f1 = t(1/sqrt(constant('pi')) * ((z1 - z2) # (mu + sqrt(2) * sigmau * w[,1]) + (z4 - z3) * sigma + z2 * a))* w[,2]; 
   f2 = t(1/sqrt(constant('pi')) * ((z1 - z2) + (zz1 - zz2) # (mu + sqrt(2) * sigmau * w[,1]) + (zz4 - zz3) * sigma + zz2 * a))* w[,2]; 

   mu = mu - sign(f2) * (f1 - mu1)/( max(abs(f2), exp(-150)));

   diff = mu - temp;
   end;
   yi = y0[e];
   chk3 = t(1/sqrt(constant('pi')) /sqrt(2*constant('pi')*sigma**2) * exp((-1)*(yi - mu - sqrt(2) * sigmau * w[,1])
             #(yi - mu - sqrt(2) * sigmau * w[,1])/(2*sigma**2)))* w[,2];
   if chk3 <= 0 then chk3 = 0.1e8;
   sum1 = sum1 + log(chk3);
end;

f = sum1;
return(f);
finish LL ;

optn={1 2};
*tc={3000 15000 . 1E-12};
tc={3000 15000 . 1E-6};
*tc={3000 15000};
* define starting va1ues for theta: either user-provided or O1S;
theta0 = &theta;
if theta0 = 0 then
do;
theta0 = t(est);
print 'Mixed Effects Model ignoring Boundary starting values for theta = '
theta0;
end;
else print 'User-supplied starting values for
theta = ' theta0;
t1 = j(1,ncol(theta0),.);
t2 = t1;
t2[ncol(theta0)] = 0;
t2[ncol(theta0)-1] = 0;
t3 = t2//t1;
con = t3;
call nlpqn(rc, theta, 'LL',theta0,optn, con, tc);
*call nlpqn(rc, theta, 'LL',theta0,optn, con, tc) grd='GRAD';
call nlpfdd(f, g, h,'LL', theta);
var = inv(-h);
sd = t(j(1,ncol(theta0),.));
temp = vecdiag(var);
do i = 1 to nrow(temp);
if temp[i] >= 0 then sd[i] = sqrt(temp[i]);
else sd[i] = .;
end; 
*sd = sqrt(vecdiag(var));
print 'Hessian = ' h ' Covariance matrix = ' var
' Standard errors = ' sd;
*print 'Hessian = ' h ' Covariance matrix = ' var
;
*print rc theta;

yx0 = yx[loc((limit ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx1 = yx[loc(( limit = .)),];
n1 = nrow(yx1);
y1 = yx1[,1];
x1 = j(n1,1,1) || yx1[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta) - 1;
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma * sigma;
sigmau = theta [k + 1];
sigmau2 = sigmau * sigmau;
sum1 = 0;
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do i = 1 to n1;
   mu1 = x1[i,]*beta;
   mu = mu1;
   diff = 2;
   do new = 1 to 200 while (abs(diff) > 1E-8);
   temp = mu;

   z1 = 1;
   z2 = probnorm((a - (mu + sqrt(2) * sigmau * w[,1]))/sigma);
   z3 = 0;
   z4 = pdf ('normal', (a - (mu + sqrt(2) * sigmau * w[,1]))/sigma, 0, 1);
   zz1 = 0;
   zz2 = (-1) * z4/sigma;
   zz3 = 0;
   zz4 = z4 # (a - (mu + sqrt(2) * sigmau * w[,1])) / sigma**2;
   f1 = t(1/sqrt(constant('pi')) * ((z1 - z2) # (mu + sqrt(2) * sigmau * w[,1]) + (z4 - z3) * sigma + z2 * a))* w[,2]; 
   f2 = t(1/sqrt(constant('pi')) * ((z1 - z2) + (zz1 - zz2) # (mu + sqrt(2) * sigmau * w[,1]) + (zz4 - zz3) * sigma + zz2 * a))* w[,2]; 

   mu = mu - sign(f2) * (f1 - mu1)/( max(abs(f2), exp(-150)));

   diff = mu - temp;
   end;
   chk1 = t(1/sqrt(constant('pi')) * probnorm((a - mu - sqrt(2) * sigmau * w[,1])/sigma))* w[,2];
   if chk1 <= 0 then chk1 = 0.1e8;
   sum1 = sum1 + log(chk1);
end;

do e = 1 to n0;
   mu1 = x0[e,]*beta;
   mu = mu1;
   diff = 2;
   do new = 1 to 200 while (abs(diff) > 1E-8);
   temp = mu;

   z1 = 1;
   z2 = probnorm((a - (mu + sqrt(2) * sigmau * w[,1]))/sigma);
   z3 = 0;
   z4 = pdf ('normal', (a - (mu + sqrt(2) * sigmau * w[,1]))/sigma, 0, 1);
   zz1 = 0;
   zz2 = (-1) * z4/sigma;
   zz3 = 0;
   zz4 = z4 # (a - (mu + sqrt(2) * sigmau * w[,1])) / sigma**2;
   f1 = t(1/sqrt(constant('pi')) * ((z1 - z2) # (mu + sqrt(2) * sigmau * w[,1]) + (z4 - z3) * sigma + z2 * a))* w[,2]; 
   f2 = t(1/sqrt(constant('pi')) * ((z1 - z2) + (zz1 - zz2) # (mu + sqrt(2) * sigmau * w[,1]) + (zz4 - zz3) * sigma + zz2 * a))* w[,2]; 

   mu = mu - sign(f2) * (f1 - mu1)/( max(abs(f2), exp(-150)));

   diff = mu - temp;
   end;
   yi = y0[e];
   chk3 = t(1/sqrt(constant('pi')) /sqrt(2*constant('pi')*sigma**2) * exp((-1)*(yi - mu - sqrt(2) * sigmau * w[,1])
             #(yi - mu - sqrt(2) * sigmau * w[,1])/(2*sigma**2)))* w[,2];
   if chk3 <= 0 then chk3 = 0.1e8;
   sum1 = sum1 + log(chk3);
end;

like = sum1*(-2);
aic = sum1*(-2) + 2 * ncol(theta);
bic = ncol(theta) * log(nrow(yx0)) + sum1*(-2);

sum1 = theta[k-1];
finalsd = sd[k-1];
postprobs=like||aic||bic||sum1||finalsd||rc||theta||t(sd);
 cname = {"Like" "AIC" "BIC" "EST" "STD" "STA" };

create ff12 from postprobs  [ colname=cname ];
append from postprobs;

quit;
run;
%mend;

%if &dataset= | &y= | &x= | &cluster =  %then %do;
    proc iml;
      print,"WARNING: NEED TO SPECIFY DATASET, COVARIATE, OUTCOME, AND CLUSTER_ID","PROGRAM WILL TERMINATE",;
	quit;
%end;

%else %do;

data tt1;
set &dataset;
if &y = &cutoff1 then vis1 = 1; else vis1 = 0;
if &y = &cutoff2 then vis2 = 1; else vis2 = 0;
run;

proc means data = tt1 noprint;
var vis1;
output out = sim1 sum(vis1) = na sum(vis2) = nb;
run;

data _null_;
set sim1;
call symput('NA', na);
call symput('NB', nb);
run;

/***********PROC MIXED Ignoring Boundary to Determinet Starting Values*********************/

ods output solutionf = tt1
CovParms = tt2;
proc mixed data = &dataset method = ml;
class &cluster;
model &y = &x/s;
random intercept/subject = &cluster;
run;

data tt21 tt22;
set tt2;
if CovParm = 'Intercept' and estimate = 0 then estimate = 1;
if _n_ = 1 then output tt22;
if _n_ = 2 then output tt21;
run;

data tt23;
set tt21 tt22;
run;

data tt3;
set tt1 tt23 (in = a);
keep estimate;
if a then estimate = sqrt(estimate);
run;

proc transpose data = tt3 out = tt4 (drop = _name_);
run;

%if  &NA gt 0  and &NB gt 0 %then %do;

%tobitcorr1(dsn=&dataset, y=&y, x=&x, cutoff1 = &cutoff1, cutoff2 = &cutoff2);
%tobitcorr2(dsn=&dataset, y=&y, x=&x, cutoff1 = &cutoff1, cutoff2 = &cutoff2);

%end;

%else %if  &NA = 0  and &NB gt 0 %then %do;

%tobitcorr3(dsn=&dataset, y=&y, x=&x, cutoff2 = &cutoff2);
%tobitcorr4(dsn=&dataset, y=&y, x=&x, cutoff2 = &cutoff2);

%end;

%else %if  &NA gt 0  and &NB = 0 %then %do;

%tobitcorr5(dsn=&dataset, y=&y, x=&x, cutoff1 = &cutoff1);
%tobitcorr6(dsn=&dataset, y=&y, x=&x, cutoff1 = &cutoff1);

%end;

%else %if  &NA = 0  and &NB = 0 %then %do;

%tobitcorr0(dsn=&dataset, y=&y, x=&x);

%end;

 %if  &NA = 0  and &NB = 0 %then %do;

data temp;
merge ff11;
if Rsta lt 0 then do;
raic = .;
rest = .;
rstd = .;
end;
keep Method raic rest rstd est_lower est_upper rp;
rest_lower = trim(left(put((rest - 1.959964 * rstd), 10.4)));
rest_upper = trim(left(put((rest + 1.959964 * rstd), 10.4)));
*rci =trim(left(put((rest - 1.959964 * rstd), 10.4)))||'-'||trim(left(put((rest + 1.959964 * rstd), 10.4))); 
rp = put((2 - 2 * probnorm(abs(rest/rstd))), 5.3);
format method $50.;
Method = 'Mixed effects model with Noboundary Data';
run;

data &out;
retain Method raic rest rstd rest_lower rest_upper rp;
SET TEMP;
RUN;

proc print data = &out noobs label;
label 
raic = "AIC"
rest = "Exposure Effects Estimate"
rstd = "Standard Error of Exposure Effects Estimate"
rest_lower = 'Estimate 95% CI Lower Limit'
rest_upper = 'Estimate 95% CI Upper Limit'
/*rci = "95% Confidence Interval of Exposure Effects Estimate"*/
rp = "p-Value";
run;

%end;

%else %do;

data temp;
set ff11 (in = a) ff12 (in = b);
format method $50.;
if a then Method = 'APV Approach';
if b then Method = 'Direct-Marginalization Approach';
if sta lt 0 then do;
aic = .;
est = .;
std = .;
end;
est_lower = trim(left(put((est - 1.959964 * std), 10.4)));
est_upper = trim(left(put((est + 1.959964 * std), 10.4)));
*CI =trim(left(put((EST - 1.959964 * STD), 10.4)))||'-'||trim(left(put((EST + 1.959964 * STD), 10.4))); 
P = put((2 - 2* probnorm(abs(EST/STD))), 5.3);
keep method AIC EST STD est_lower est_upper P;
run;

data &out;
retain Method AIC est std est_lower est_upper P;
set temp;
run;

proc print data = &out noobs label;
label AIC = "AIC"
Est = "Overall Exposure Effects Estimate"
StD = "Standard Error of Overall Exposure Effects Estimate"
est_lower = 'Estimate 95% CI Lower Limit'
est_upper = 'Estimate 95% CI Upper Limit'
/*CI = "95% Confidence Interval of Overall Exposure Effects Estimate"*/
P = "p-Value";
run;

%end;

%end;

%mend;

/* Gauss-Hermite quadrature Weights and Abscissae Data Set */

data abswei;
input k abs weight;
datalines;
1	-8.098761139250850052013	2.5910437138470814735E-29
2	-7.411582531485468809439	8.544056963775510774E-25
3	-6.84023730524935541785	2.5675933654116696605E-21
4	-6.328255351220081955657	1.9891810121165024856E-18
5	-5.85409505603040010804	6.0083587894908166903E-16
6	-5.4066542479701276084	8.80570764521613225662E-14
7	-4.979260978545255871627	7.1565280526903187084E-12
8	-4.567502072844394855169	3.52562079136541190292E-10
9	-4.16825706683250020154	1.12123608322758101745E-8
10	-3.779206753435223493119	2.41114416367052344179E-7
11	-3.398558265859628346294	3.6315761506930235118E-6
12	-3.024879883901284437677	3.93693398109249277043E-5
13	-2.656995998442895794981	3.13853594541331475647E-4
14	-2.293917141875083421885	0.001871496829597952779484
15	-1.934791472282295793298	0.00846088800825813243994
16	-1.578869894931613886258	0.0293125655361723698457
17	-1.225480109046289030949	0.0784746058654043913089
18	-0.8740066123570880774379	0.1633787327132714570815
19	-0.5238747138322771926149	0.265728251877377076143
20	-0.1745372145975823834895	0.338643277425589218202
21	0.174537214597582383489	0.3386432774255892182024
22	0.523874713832277192615	0.265728251877377076143
23	0.874006612357088077438	0.163378732713271457082
24	1.225480109046289030949	0.0784746058654043913089
25	1.578869894931613886258	0.0293125655361723698457
26	1.934791472282295793298	0.00846088800825813243994
27	2.293917141875083421885	0.00187149682959795277948
28	2.656995998442895794981	3.13853594541331475647E-4
29	3.024879883901284437677	3.93693398109249277043E-5
30	3.398558265859628346294	3.6315761506930235118E-6
31	3.77920675343522349312	2.41114416367052344179E-7
32	4.168257066832500201536	1.12123608322758101745E-8
33	4.56750207284439485517	3.5256207913654119029E-10
34	4.97926097854525587163	7.1565280526903187084E-12
35	5.4066542479701276084	8.8057076452161322566E-14
36	5.85409505603040010804	6.0083587894908166903E-16
37	6.328255351220081955657	1.9891810121165024856E-18
38	6.84023730524935541785	2.5675933654116696605E-21
39	7.41158253148546880944	8.54405696377551077388E-25
40	8.098761139250850052013	2.5910437138470814735E-29

;
run;

proc import datafile = 'C:\Users\Wei Wang\Box Sync\Rand Tob Overall\programs\SAS Macro Upload\example.xls' out = ttemp1 replace;
run; 

proc mixed data = ttemp1 method = ml;
class cluster_id;
model outcome = covariate1 exposure/s;
random intercept/subject = cluster_id;
run;

proc means data = ttemp1;
var outcome;
run;

%tobitcorr(dataset=ttemp1, y=outcome, x=covariate1 exposure, cutoff1 =., cutoff2 =100, cluster = cluster_id, out =out1);
%tobitcorr(dataset=ttemp1, y=outcome, x=covariate1 exposure, cutoff1 =50, cutoff2 =., cluster = cluster_id, out =out2);
%tobitcorr(dataset=ttemp1, y=outcome, x=covariate1 exposure, cutoff1 =50, cutoff2 =100, cluster = cluster_id, out =out3);
