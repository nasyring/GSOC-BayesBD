// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-


// we only include RcppArmadillo.h which pulls Rcpp.h in for us

#include "RcppArmadillo.h"


using namespace Rcpp;
using namespace arma;
using namespace std;
// via the depends attribute we tell Rcpp to create hooks for

// RcppArmadillo so that the build process will know what to do
//

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]

Rcpp::List BayesBDbinary(SEXP & obs, SEXP & inimean, SEXP & nrun, SEXP & nburn, SEXP & J, SEXP & ordering, SEXP & slice, SEXP & outputAll) { 

// Including additional user-defined Rcpp functions
RNGScope scp;
Rcpp::Function eigenfun("eigenfun");
Rcpp::Function besselIs("besselIs");
Rcpp::Function unisliceL("unisliceL");

// Extracting the inputs
Rcpp::List 		      	obsL(obs);
arma::colvec thetaobs 	 	= Rcpp::as<arma::colvec>(obsL["theta.obs"]);
arma::colvec demean 	 	= Rcpp::as<arma::colvec>(obsL["r.obs"]);
arma::colvec intensityobs       = Rcpp::as<arma::colvec>(obsL["intensity"]);
Rcpp::CharacterVector order     = Rcpp::as<CharacterVector>(ordering);
Rcpp::CharacterVector O     = Rcpp::CharacterVector(1);O(0)="O";
Rcpp::CharacterVector I     = Rcpp::CharacterVector(1);I(0)="I";
int i_nrun 				= Rcpp::as<int>(nrun);
int i_nburn				= Rcpp::as<int>(nburn);
int i_J				= Rcpp::as<int>(J);
double d_inimean 			= Rcpp::as<double>(inimean); 
Rcpp::LogicalVector output_All = Rcpp::as<LogicalVector>(outputAll);
Rcpp::LogicalVector slice_ = Rcpp::as<LogicalVector>(slice);

// Declaring variables 
Rcpp::List temp2del;
Rcpp::List lambdalist;
Rcpp::List result;

int s					= thetaobs.size();
int L                  		= 2*i_J+1;
int betatau           		= 1;
double alphatau                 = 100.0;
int alpha_a      		= 2;
int beta_a       		= 1;
double alpha_1                  = 0.0;
double beta_1                   = 0.0;
int nin = 0;
int ninone = 0;
int nout = 0;
int noutone = 0;

arma::colvec tmpmatk		= arma::colvec(s); tmpmatk.fill(0.0);
arma::colvec mu			= arma::colvec(s); mu.fill(d_inimean);
arma::colvec diffini 		= demean - mu; 
arma::colvec eigencnk   	= arma::colvec(1); eigencnk.fill(0.0);
arma::colvec ank   	      = arma::colvec(1); ank.fill(0.0); 
arma::colvec astar   	      = arma::colvec(1); astar(0)=0.0;
arma::colvec bstar   	      = arma::colvec(1); bstar(0)=0.0;
arma::colvec tauinirg   	= arma::colvec(1); tauinirg(0)=100.0;
arma::colvec x1 			= arma::colvec(1); x1.fill(0.0); 
arma::colvec x2 			= arma::colvec(1); x2.fill(0.0);
arma::colvec piin1 		= arma::colvec(1); piin1.fill(0.0);
arma::colvec piin2 		= arma::colvec(1); piin2.fill(0.0);
arma::colvec piout1 		= arma::colvec(1); piout1.fill(0.0);
arma::colvec piout2 		= arma::colvec(1); piout2.fill(0.0);
arma::colvec piinini 		= arma::colvec(1); piinini.fill(0.0);
arma::colvec pioutini 		= arma::colvec(1); pioutini.fill(0.0);
arma::colvec calc1 		= arma::colvec(1); calc1.fill(0.0);
arma::colvec alphaini 		= arma::colvec(1); alphaini.fill(0.0);
arma::colvec betaini 		= arma::colvec(1); betaini.fill(0.0);
arma::colvec lambdaini 		= arma::colvec(1); lambdaini(0)=1.0;
arma::colvec anini 		= arma::colvec(L); anini.fill(0.0);
arma::colvec cnini            = arma::colvec(i_J+1); cnini.fill(0.0);
arma::colvec eigencnini       = arma::colvec(L); eigencnini.fill(0.0);
arma::colvec tmp2del2         = arma::colvec(1); tmp2del2.fill(0.0);
arma::colvec interim          = arma::colvec(1); interim.fill(0.0);
arma::colvec rmat             = arma::colvec(s); rmat.fill(0.0);
arma::colvec tmpp             = arma::colvec(s); tmpp.fill(0.0);
arma::colvec tmp              = arma::colvec(s); tmp.fill(0.0);
arma::colvec gx0              = arma::colvec(1); gx0.fill(0.0);
arma::colvec ansum		= arma::colvec(L); ansum.fill(0.0);
arma::mat ansmp			= arma::mat(L,i_nrun); ansmp.fill(0.0);
arma::mat pismp			= arma::mat(2,i_nrun); pismp.fill(0.0);
arma::colvec pivec		= arma::colvec(2);
arma::mat tmpmat              = arma::mat(s,L); tmpmat.fill(0.0);

// Variables used for slice and MH sampling
arma::colvec m = arma::colvec(1);m.fill(10000.0);
NumericVector rexpval;
arma::colvec logy = arma::colvec(1); 
arma::colvec xx1 = arma::colvec(1); xx1.fill(0.0);
arma::colvec LL = arma::colvec(1); 
arma::colvec R = arma::colvec(1); 
arma::colvec g = arma::colvec(1); g.fill(0.0);
arma::colvec diffnew = arma::colvec(s);diffnew.fill(0.0);
arma::colvec Lv = arma::colvec(s);Lv.fill(LL(0));
arma::colvec Rv = arma::colvec(s);Rv.fill(R(0));
arma::colvec xv = arma::colvec(s);xv.fill(0.0);
arma::colvec ankv = arma::colvec(s);ankv.fill(ank(0));
arma::colvec count = arma::colvec(1); count.fill(0.0);
bool res1 = TRUE;
bool res2 = FALSE;
bool res3 = FALSE;
int g1 = 0;
int g3 = 0;
NumericVector anok; NumericVector ankk;
arma::colvec lLnew = arma::colvec(1);
arma::colvec r = arma::colvec(1);
arma::colvec uu = arma::colvec(1);


// Filling tmpmat with values for each of L eigenfunctions for each of s observations
for (int i=0; i<s; i++) {
	for(int j = 0; j<L;j++){
            tmpmat(i,j) = Rcpp::as<double>(eigenfun(j+1,thetaobs[i]));
	}
}


for(int i = 0; i<s; i++){
	if(demean(i) < mu(i)){
		if(intensityobs(i) == 1.0){
			piin1(0) = piin1(0) + 1.0;
		}
	}
}

for(int i = 0; i<s; i++){
	if(demean(i) < mu(i)){
		piin2(0) = piin2(0) + 1.0;
	}
}

for(int i = 0; i<s; i++){
	if(demean(i) >= mu(i)){
		if(intensityobs(i) == 1.0){
			piout1(0) = piout1(0) + 1.0;
		}
	}
}

for(int i = 0; i<s; i++){
	if(demean(i) >= mu(i)){
		piout2(0) = piout2(0) + 1.0;
	}
}

piinini = piin1 / piin2;
pioutini = piout1 / piout2;


calc1 = 2*M_PI*std::exp(-2.0*lambdaini(0));
for(int i=0; i<(i_J+1); i++){
	cnini(i) = calc1(0)*Rcpp::as<double>(besselIs(max(2*lambdaini,0.0),i));
}

eigencnini(0) = cnini(0);
for(int i=1; i<(i_J+1); i++){
	eigencnini(2*i-1) = cnini(i);
	eigencnini(2*i) = cnini(i);
} 
alphaini = log(piinini*(1-pioutini)/(pioutini*(1-piinini)));
betaini = log((1 - piinini)/(1 - pioutini));
for(int i = 0; i<s; i++){
	if(diffini(i)<0){
		tmp(i) = 1;
	}
}

nin = sum(tmp);
for(int i=0; i<s; i++){
	if(intensityobs(i)==1){
	ninone = ninone + tmp(i);
	}
}

//main loop
for(int i=0; i<(i_nrun+i_nburn); i++){
for(int k=0; k<L; k++){
	eigencnk(0) = eigencnini(k);
	ank(0) = anini(k);
	tmpmatk = tmpmat.col(k);
	gx0 = ninone*alphaini+nin*betaini-(((pow(ank,2))/eigencnk)*tauinirg/2);
	
	//Slice sampling of kth element of z
	if(slice_[0] == TRUE){
 		rexpval = Rcpp::rexp(1);
		logy(0) = gx0(0) - rexpval(0);
		xx1.fill(0.0);
		rexpval = Rcpp::runif(1, 0, 1.0);
		LL(0) = ank(0) - rexpval[0];
		R = LL + 1.0;
		g.fill(0.0);
		diffnew.fill(0.0);
		Lv.fill(LL(0));
		Rv.fill(R(0));
		xv.fill(0.0);
		ankv.fill(ank(0));
		count.fill(0.0);
		res1 = TRUE;
		res2 = FALSE;
		res3 = FALSE;
		g1 = 0;
		g3 = 0;
		while( res1 ){
            	count = count+1;
			res1 = count(0) < 100.0;
			res2 = LL(0)<=-10000;
			if(res2){
				LL(0) = LL(0)+1;
				res1 = FALSE;
				}
			else {
				Lv.fill(LL(0));
				diffnew = diffini - ((Lv - ankv) % tmpmatk);
				g1 = 0; g3 = 0;
				for(int i=0; i<s;i++){
					if(diffnew(i)<=0){
						g3 = g3+1;
						if(intensityobs(i)==1){
							g1 = g1+1;
						}
					}
				}
			g = g1*alphaini+g3*betaini-((pow(LL,2))/eigencnk)*tauinirg/2;
       		res3 = g(0) <= logy(0);
			}
          		if (res3){ 
			LL(0) = LL(0) + 1.0;
			res1=FALSE;}
           		LL(0) = LL(0) - 1.0;
		}

		res2 = FALSE;
		res3 = FALSE;
		count = 0.0;
		res1 = TRUE;
      	while( res1 ){
            	count=count+1.0;
			res1 = count(0) < 100.0;
			res2 = R(0)>=10000.0;
			if(res2){
				R(0) = R(0) - 1.0;
				res1 = FALSE;
			}
			else {
				Rv.fill(R(0));
				diffnew = diffini - ((Rv - ankv) % tmpmatk);
				g1 = 0; g3 = 0;
				for(int i=0; i<s;i++){
					if(diffnew(i)<=0){
						g3 = g3+1;
						if(intensityobs(i)==1){
							g1 = g1+1;
						}
					}
				}
				g = g1*alphaini+g3*betaini-((pow(R,2))/eigencnk)*tauinirg/2;
				res3 = g(0) <= logy(0);
      		}      
			if (res3){ 
            		R(0)=R(0)-1.0;
				res1 = FALSE;}
            		R(0) = R(0) + 1.0;
          		}
			res1 = (LL(0)<-10000.0);
    			if (res1) {
       			LL(0) = -10000.0;
   			}
			res1 = (R(0)>10000.0);
    			if (res1) {
        			R(0) = 10000.0;
    		}
		count = 0.0;
		res1 = TRUE;
		res2 = FALSE;
		res3 = FALSE;
    		while(res1) {
			count=count+1;
			res1 = count(0) < 200.0 ;
			if(res2){
				res1 = FALSE;
			}
			else{
        			xx1 = Rcpp::runif(1, LL(0), R(0));
				xv.fill(xx1(0));
				diffnew = diffini - ((xv - ankv) % tmpmatk);
	     			g1 =0; g3=0;
              		for(int i=0; i<s;i++){
					if(diffnew(i)<=0){
						g3 = g3+1;
					if(intensityobs(i)==1){
						g1 = g1+1;
					}
				}
			}
        		g = g1*alphaini+g3*betaini-((pow(xx1,2))/eigencnk)*tauinirg/2;
	   		res3 = (xx1(0) > ank(0));
			}
 			res2 = (g(0) >= logy(0));
        		if (res3) {
            		R(0) = xx1(0);
        		}
        		else {
            		LL(0) = xx1(0);
        		}
    		}  

		anini(k) = xx1(0);
		ninone = g1;
		nin = g3;
		diffini = diffnew;
	}// begin MH sampling of kth element of z
	else{
		g1 = 0; g3 = 0;											
		anok[0] = anini(k);
		ankk = Rcpp::rnorm(1, anini(k), 0.05); 
		Lv.fill(ankk[0]);
		ankv.fill(anok[0]);
		diffnew = diffini - ((Lv - ankv) % tmpmatk);
		g1 = 0; g3 = 0;
		for(int i=0; i<s;i++){
			if(diffnew(i)<=0){
				g3 = g3+1;
				if(intensityobs(i)==1){
					g1 = g1+1;
				}
			}
		}
		lLnew = g1*alphaini+g3*betaini-((pow(ankk[0],2))/eigencnk)*tauinirg/2;
		r = Rcpp::dnorm(ankk, anok[0],.05) - Rcpp::dnorm(anok,ankk[0],.05);
		r = r + lLnew - gx0;
      	r(0) = fmin(std::exp(r(0)), 1.0);
		uu = Rcpp::runif(1);
      	if(uu(0) <= r(0)) {
			ninone = g1;
			nin = g3;
			diffini = diffnew;
			anini(k) = ankk[0];									// update z, nin, and loglikelihood if sample is accepted
			gx0 = lLnew;
      	}
	}
}// end of sampling of z

astar = alphatau + L/2;
interim = anini.t()*(anini / eigencnini);
bstar = betatau + interim/2;
tauinirg = Rcpp::rgamma(1,astar(0),bstar(0));
rmat.fill(0);
rmat = tmpmat*anini + mu;
for(int j = 0; j< s; j++){
	if(demean(j)<=rmat(j)){tmpp(j)=1.0;}
	else{tmpp(j)=0.0;}
}
nin = sum(tmpp);
ninone = sum(intensityobs%tmpp);
nout = s - nin;
noutone = sum(intensityobs) - ninone;

x1 = Rcpp::rbeta(1, alpha_1 + ninone, beta_1 + nin - ninone);
x2 = Rcpp::rbeta(1, alpha_1 + noutone, beta_1 + nout - noutone);

if(order(0) == I(0)){
piinini = arma::max(x1, x2);
pioutini = arma::min(x1, x2);
}else if(order(0) == O(0)){
piinini = arma::min(x1, x2);
pioutini = arma::max(x1, x2);
} else {
piinini = x1;
pioutini = x2;
}

alphaini = log(piinini * (1 - pioutini)/(pioutini *(1 - piinini)));
betaini = log((1 - piinini)/(1 - pioutini));

gx0 = -1/2 * (sum(log(eigencnini)) + tauinirg * interim) + (alpha_a - 1) * log(lambdaini) - beta_a;

lambdalist = unisliceL(lambdaini, gx0, i_J, tauinirg, anini,  alpha_a,  beta_a, lambdaini);
lambdaini = Rcpp::as<NumericVector>(lambdalist["x1"]);

if(i > i_nburn){
ansum = ansum + anini;
ansmp.col(i-i_nburn) = anini;
pivec(0) = piinini(0); pivec(1) = pioutini(0);
pismp.col(i-i_nburn) = pivec; 
}

}

ansum = ansum/i_nrun;

arma::colvec thetaplot = arma::colvec(200);
for(int k = 0; k<200; k++){
	thetaplot(k) = 2.0*M_PI*(k/200.0); 
}
arma::mat estfunc = arma::mat(200,L);
arma::colvec esttheta = arma::colvec(200);
double est = 0.0;

for(int k = 0; k<200; k++){
		est=0.0;
  		for(int j = 0; j<L; j++){
			estfunc(k,j) = Rcpp::as<double>(eigenfun(j+1,thetaplot(k)));
			est = est + Rcpp::as<double>(eigenfun(j+1,thetaplot(k)))*ansum(j)+(d_inimean/L);
		}
		esttheta(k) =est;
}

arma::mat estthetapts = arma::mat(200,i_nrun);
arma::mat estpts = arma::mat(200,i_nrun);
arma::colvec mutheta = arma::colvec(200);mutheta.fill(d_inimean);

for(int i = 0; i<i_nrun; i++){
		estpts.col(i) = estfunc*ansmp.col(i)+mutheta;
		estthetapts.col(i) = esttheta;
}
arma::colvec variance = arma::colvec(200);

variance = sum(pow(estpts-estthetapts,2),1)/(i_nrun-1.0);
arma::colvec maxed = arma::colvec(1);
arma::colvec diffed = arma::colvec(1);
arma::colvec maxes = arma::colvec(i_nrun);
arma::colvec sortval = arma::colvec(200);
arma::colvec sorted = arma::colvec(200);
for(int j = 0; j<i_nrun; j++){
maxed(0)=0.0;
for(int i = 0; i<200; i++){
diffed = abs(estpts(i,j)-estthetapts(i,j))/pow(variance(i),0.5);
maxed(0) = fmax(maxed(0),diffed(0));
}
maxes(j) = maxed(0);
}
sorted = sort(maxes);
sortval.fill(sorted(floor(0.95*i_nrun)));
arma::colvec lower = arma::colvec(200);
arma::colvec upper = arma::colvec(200);
lower = esttheta - sortval%pow(variance,0.5);
upper = esttheta + sortval%pow(variance,0.5);

if(output_All[0] == TRUE){result = Rcpp::List::create(Rcpp::Named("estimate") = esttheta,Rcpp::Named("theta") = thetaplot,Rcpp::Named("lower") = lower,Rcpp::Named("upper") = upper, Rcpp::Named("pi.smp") = pismp, Rcpp::Named("coef.smp") = ansmp);}
else{result = Rcpp::List::create(Rcpp::Named("estimate") = esttheta,Rcpp::Named("theta") = thetaplot,Rcpp::Named("lower") = lower,Rcpp::Named("upper") = upper);}

return result;
}

// [[Rcpp::export]]

Rcpp::List unisliceL(SEXP & ix0, SEXP & igx0, SEXP & ii_J, SEXP & itauini, SEXP & ianini, SEXP & ialpha_a, SEXP & ibeta_a, SEXP & ilambdaini){
RNGScope scope;
Rcpp::Function besselIs("besselIs");
arma::colvec x0 = Rcpp::as<arma::colvec>(ix0);
arma::colvec gx0 = Rcpp::as<arma::colvec>(igx0);
int i_J = Rcpp::as<int>(ii_J);
arma::colvec tauini = Rcpp::as<arma::colvec>(itauini);
arma::colvec anini = Rcpp::as<arma::colvec>(ianini);
int alpha_a = Rcpp::as<int>(ialpha_a);
int beta_a = Rcpp::as<int>(ibeta_a);
arma::colvec lambdaini = Rcpp::as<arma::colvec>(ilambdaini);
Rcpp::NumericVector rndsmp = Rcpp::NumericVector( Rcpp::Dimension(1)); rndsmp = Rcpp::rexp(1);
arma::colvec logy = arma::colvec(1); logy(0) = gx0(0) - rndsmp[0];
rndsmp = Rcpp::runif(1, 0, 1.0);
arma::colvec L = arma::colvec(1); L(0) = lambdaini(0) - rndsmp[0];
arma::colvec R = arma::colvec(1); R(0) = L(0) + 1.0;
arma::colvec x1 = arma::colvec(1);x1(0)=0.0;
arma::colvec g = arma::colvec(1); g(0)=0.0;
arma::colvec count = arma::colvec(1); count(0)=0.0;
arma::colvec interim = arma::colvec(1);interim(0) = 0.0;
arma::colvec cnini = arma::colvec(i_J+1);
arma::colvec eigencnini = arma::colvec(2*i_J+1);
arma::colvec calc1 = arma::colvec(1);calc1(0) = 0.0;
bool res1 = TRUE;
bool res2 = FALSE;
bool res3 = FALSE;


while( res1 ){
		interim(0)=0.0;
            count = count+1.0;
		res1 =  count(0) < 20.0 ;
	      res2 = L(0)<=0.0;
            if (res2){ 
			L(0) = L(0)+1.0;
			res1= FALSE;}
		else {
			calc1 = 2*M_PI*std::exp(-2.0*L(0));
			for(int i=0; i<(i_J+1); i++){
				cnini(i) = calc1(0)*Rcpp::as<double>(besselIs(2*L,i));
			}
			eigencnini(0) = cnini(0);
			for(int i=1; i<(i_J+1); i++){
				eigencnini(2*i-1) = cnini(i);
				eigencnini(2*i) = cnini(i);
			} 
			interim = anini.t()*(anini / eigencnini);
			g = -0.5 * (sum(log(eigencnini)) + tauini * interim) + (alpha_a - 1) * log(L) - beta_a;
			res3 = g(0) <= logy(0);
		}
            if (res3){ 
			L(0) = L(0)+1.0;
			res1=FALSE;}
            L(0) = L(0) - 1.0;
           }
		 count = 0.0;
		res1 = TRUE;res2 = FALSE;res3 = FALSE;
       while( res1 ){
		interim(0)=0.0;
            count = count+1.0;
		res1 =  count(0) < 20.0 ;
		res2 = R(0)>=500.0;
            if (res2){ 
			R(0) = R(0)-1.0;
			res1= FALSE;}
		else {
			calc1 = 2*M_PI*std::exp(-2.0*R(0));
			for(int i=0; i<(i_J+1); i++){
				cnini(i) = calc1(0)*Rcpp::as<double>(besselIs(2*R,i));
			}
			eigencnini(0) = cnini(0);
			for(int i=1; i<(i_J+1); i++){
				eigencnini(2*i-1) = cnini(i);
				eigencnini(2*i) = cnini(i);
			}
			interim = anini.t()*(anini / eigencnini);
			g = -0.5 * (sum(log(eigencnini)) + tauini * interim) + (alpha_a - 1) * log(R) - beta_a;
			res3 = (g(0) <= logy(0));
		}
            if (res3){ 
			R(0) = R(0)-1.0;
			res1=FALSE;}
            R(0) = R(0) + 1.0;
            }
		
	res1 = L(0)<=0.0;
    if (res1) {
        L(0) = 0.0;
    }
	res1 = R(0)>500.0;
    if (res1) {
        R(0) = 500.0;
    }
	count = 0.0;
	res1 = TRUE;res2 = FALSE;res3 = FALSE;
    while(res1) {
		interim(0)=0.0;
		count=count+1;
		res1 =  count(0) < 100.0  ;
       	if (res2){ 
			res1 = FALSE;
			}
		else {
			x1 = Rcpp::runif(1, L(0), R(0));
			calc1 = 2*M_PI*std::exp(-2.0*x1(0));
			for(int i=0; i<(i_J+1); i++){
				cnini(i) = calc1(0)*Rcpp::as<double>(besselIs(2*x1,i));
			}
			eigencnini(0) = cnini(0);
			for(int i=1; i<(i_J+1); i++){
				eigencnini(2*i-1) = cnini(i);
				eigencnini(2*i) = cnini(i);
			}
			interim = anini.t()*(anini / eigencnini);
			g = -0.5 * (sum(log(eigencnini)) + tauini * interim) + (alpha_a - 1) * log(x1) - beta_a;
			res3 = (x1(0) > x0(0));
		}
	  	  	res2 = (g(0) >= logy(0));
        if (res3) {
            R(0) = x1(0);
        }
        else {
            L(0) = x1(0);
        }
    }  
 	return Rcpp::List::create(Rcpp::Named("x1") = x1, Rcpp::Named("count") = count);
}

// [[Rcpp::export]]

double eigenfun(SEXP & in, SEXP & ix) { 
	int n = Rcpp::as<int>(in); double x = Rcpp::as<double>(ix);
	int k1 = n%2;
    	double k2 = (n - k1)/2;
	double ret = 0;
   	if (n == 1) {
       		 ret = 1/sqrt(2 * M_PI) + 0 * x;
    	}
    	if (n > 1) {
        	if (k1 == 0) {
            		 ret = 1/sqrt(M_PI) * cos(k2 * x);
        	}
       		if (k1 == 1) {
            		 ret = 1/sqrt(M_PI) * sin(k2 * x);
       		}
    }
    return ret;
}

// [[Rcpp::export]]

Rcpp::List BayesBDnormal(SEXP & obs, SEXP & inimean, SEXP & nrun, SEXP & nburn, SEXP & J, SEXP & ordering_mu, SEXP & ordering_sigma, SEXP & slice, SEXP & outputAll) { 

// Including additional user-defined Rcpp functions
RNGScope scp;
Rcpp::Function eigenfun("eigenfun");
Rcpp::Function besselIs("besselIs");


// Extracting the inputs
Rcpp::List 		      	obsL(obs);								// the observed (or simulated) data
arma::colvec thetaobs 	 	= Rcpp::as<arma::colvec>(obsL["theta.obs"]);			// angle in polar coords of datapoints
arma::colvec demean 	 	= Rcpp::as<arma::colvec>(obsL["r.obs"]);			// radius in polar coords of datapoints
arma::colvec intensityobs       = Rcpp::as<arma::colvec>(obsL["intensity"]);		// image intensity of datapoints
Rcpp::CharacterVector order_mu  = Rcpp::as<CharacterVector>(ordering_mu);
Rcpp::CharacterVector order_sd = Rcpp::as<CharacterVector>(ordering_sigma);
Rcpp::CharacterVector O     = Rcpp::CharacterVector(1);O(0)="O";
Rcpp::CharacterVector I     = Rcpp::CharacterVector(1);I(0)="I";
int i_nrun 				= Rcpp::as<int>(nrun);						// number of posterior samples to be retained
int i_nburn				= Rcpp::as<int>(nburn);						// number of posterior samples to be burned
int i_J				= Rcpp::as<int>(J);						// 2*J+1 is number of eigenfunctions to use in basis expansion of image boundary
double d_inimean 			= Rcpp::as<double>(inimean);					// radius of a circle used for initial value of image boundary
Rcpp::LogicalVector output_All = Rcpp::as<LogicalVector>(outputAll);
Rcpp::LogicalVector slice_ = Rcpp::as<LogicalVector>(slice);

// Declaring variables 
Rcpp::List result;
int s					= thetaobs.size();						// number of datapoints, observations						
int L                  		= 2*i_J+1;								// 2*J+1 is number of eigenfunctions to use in basis expansion of image boundary
int betatau           		= 1;
double alphatau                 = 500.0;									// hyperparameter for updating tau
int alpha_a      		= 2;									// hyperparameters for updataing lambda
int beta_a       		= 1;
int nin = 0;												// count of the number of datapoints inside the image boundary
int nout = 0;	
double alpha_2 = 0.01; double beta_2 = 0.01;											// count of the number of datapoints outside the image boundary


arma::colvec mu			= arma::colvec(s); mu.fill(d_inimean);			// vector with initial guess of image boundary (just a circle of radius d_inimean)
arma::colvec diffini          = demean - mu;
arma::colvec ank   	      = arma::colvec(1); ank.fill(0.0); 				
arma::colvec astar   	      = arma::colvec(1); astar(0)=0.0;				// hyperparameters for sampling tauinirg
arma::colvec bstar   	      = arma::colvec(1); bstar(0)=0.0;
arma::colvec tauinirg   	= arma::colvec(1); tauinirg(0)=500.0;			// tau parameter from Li and Ghosal
arma::colvec mu0 		      = arma::colvec(1); mu0 = sum(intensityobs)/s;		// Next several parameters all relate to calculation of mu1, mu2, sd1, sd2, the parameters of the Gaussian likelihoods for intensities within and without the image boundary
arma::colvec sd0 		      = arma::colvec(1); sd0.fill(1000.0);
arma::colvec mu1 		      = arma::colvec(1); mu1(0) = mu0(0);
arma::colvec mu2 		      = arma::colvec(1); mu2(0) = mu0(0);
arma::colvec sd1 		      = arma::colvec(1); sd1(0) = sd0(0);
arma::colvec sd2 		      = arma::colvec(1); sd2(0) = sd0(0);
arma::colvec tempmu1 		      = arma::colvec(1); tempmu1(0) = 0.0;
arma::colvec tempmu2 		      = arma::colvec(1); tempmu2(0) = 0.0;
arma::colvec tempsd1 		      = arma::colvec(1); tempsd1(0) = 0.0;
arma::colvec tempsd2 		      = arma::colvec(1); tempsd2(0) = 0.0;
arma::colvec xbarin 		= arma::colvec(1); xbarin.fill(0.0);
arma::colvec xbarout 		= arma::colvec(1); xbarout.fill(0.0);
arma::colvec xinsum 		= arma::colvec(1); xinsum.fill(0.0);
arma::colvec xoutsum 		= arma::colvec(1); xoutsum.fill(0.0);
arma::colvec xinsumsq 		= arma::colvec(1); xinsumsq.fill(0.0);
arma::colvec xoutsumsq 		= arma::colvec(1); xoutsumsq.fill(0.0);
arma::colvec ssin 		= arma::colvec(1); ssin.fill(0.0);
arma::colvec ssout 		= arma::colvec(1); ssout.fill(0.0);
arma::colvec calc1 		= arma::colvec(1); calc1.fill(0.0);				// intermediary calculation of Bessel functions for eigencnini
arma::colvec std_sum          = arma::colvec(1); std_sum.fill(0.0);                 
arma::colvec lambdaini        = arma::colvec(1); lambdaini(0) = 1.0;                // the a parameter, see Section 5: Sampling Algorithms in Li and Ghosal 2015 for details
arma::colvec anini 		= arma::colvec(L); anini.fill(0.0);				// the z parameter (eigenfunction coefficients); see Section 5: Sampling Algorithms in Li and Ghosal 2015 for details
arma::colvec cnini            = arma::colvec(i_J+1); cnini.fill(0.0);			// computation of eigencnini
arma::colvec eigencnini       = arma::colvec(L); eigencnini.fill(0.0);			// computation of Bessel functions
arma::colvec interim          = arma::colvec(1); interim.fill(0.0);			// intermediary calculation for gx0
arma::colvec rmat             = arma::colvec(s); rmat.fill(0.0);				// matrix of radii at each observed datapoint for an estimate of image boundary
arma::colvec tmpp             = arma::colvec(s); tmpp.fill(0.0);				// used to count nin and nout
arma::colvec sdcomp           = arma::colvec(s); sdcomp.fill(0.0);			// computation of gxZ
arma::colvec gx0              = arma::colvec(1); gx0.fill(0.0);				// loglikelihood of lambdaini parameter
arma::colvec gxZ              = arma::colvec(1); gxZ.fill(0.0);				// loglikelihood of eigenfunction coefficients, z parameter, called anini
arma::colvec ansum		= arma::colvec(L); ansum.fill(0.0);				// averaging the eigenfunction coefficients (anini) over posterior samples
arma::mat ansmp			= arma::mat(L,i_nrun); ansmp.fill(0.0);
arma::mat musigsmp		= arma::mat(4,i_nrun); musigsmp.fill(0.0);
arma::colvec musigvec		= arma::colvec(4);
arma::mat tmpmat              = arma::mat(s,L); tmpmat.fill(0.0);				// Fixed matrix of eigenfunction values at each datapoint location

Rcpp::NumericVector rndsmp    = Rcpp::NumericVector( Rcpp::Dimension(1)); 		// variables for slice sampling of z (anini) and a (lambdaini)
arma::colvec logy             = arma::colvec(1); 
arma::colvec LL               = arma::colvec(1); 
arma::colvec R                = arma::colvec(1); 
arma::colvec x1               = arma::colvec(1);
arma::colvec g                = arma::colvec(1); 
arma::colvec gxF              = arma::colvec(1); 
arma::colvec count            = arma::colvec(1); 
int ninnew                    = 0; 
arma::colvec annew; 
bool res1                     = TRUE;
bool res2                     = FALSE;
bool res3                     = FALSE;

NumericVector ankk;											// variables for MH sampling of z (anini)
NumericVector anok; 
arma::colvec lLnew = arma::colvec(1); 
arma::colvec r = arma::colvec(1); 
arma::colvec uu = arma::colvec(1);  



// Filling tmpmat with values for each of L eigenfunctions for each of s obs
for (int i=0; i<s; i++) {
	for(int j = 0; j<L;j++){
            tmpmat(i,j) = Rcpp::as<double>(eigenfun(j+1,thetaobs[i]));
	}
}


calc1 = 2.0*M_PI*std::exp(-2.0*lambdaini(0));
for(int i=0; i<(i_J+1); i++){
	cnini(i) = calc1(0)*Rcpp::as<double>(besselIs(2.0*lambdaini,i));
}

eigencnini(0) = cnini(0);
for(int i=1; i<(i_J+1); i++){
	eigencnini(2*i-1) = cnini(i);
	eigencnini(2*i) = cnini(i);
} 

for(int i = 0; i<s; i++){
	if(diffini(i)<0){
		tmpp(i) = 1;
	}else{tmpp(i)=0;}
}

nin = sum(tmpp);
nout = s - nin;

mu1 = tmpp.t()*intensityobs/nin; mu2 = intensityobs.t()*(1.0-tmpp)/nout;
ssin = trans(intensityobs % tmpp - mu1(0)* tmpp)*(intensityobs % tmpp - mu1(0)* tmpp);
ssout = trans(intensityobs % (1.0-tmpp) - mu2(0)* (1.0-tmpp))*(intensityobs % (1.0-tmpp) - mu2(0)* (1.0-tmpp));
sd1 = pow(ssin/(nin-1), 0.5);
sd2 = pow(ssout/(nout-1), 0.5);

sdcomp = (intensityobs - (tmpp*mu1(0) + (1.0-tmpp)*mu2(0)))/(tmpp*sd1(0) + (1.0-tmpp)*sd2(0));
std_sum = -sdcomp.t()*sdcomp*0.5;
gxZ = nin * (-log(sd1(0)))-nin*(-log(sd2(0))) + std_sum - anini.t()*diagmat(1.0/eigencnini)*anini*tauinirg*0.5;

//main loop for posterior sampling
for(int i=0; i<(i_nrun+i_nburn); i++){
for(int k=0; k<L; k++){
	if(slice_[0] == FALSE){
	tmpp.fill(0.0);											// begin MH sampling of kth element of z
	anok[0] = anini(k);
	annew = anini; 
	ankk = Rcpp::rnorm(1, anini(k), 0.05); 
	annew(k) = ankk[0];
	rmat = tmpmat*annew + mu;
	for(int i = 0; i<s; i++){
		if(demean(i)<rmat(i)){
			tmpp(i) = 1.0;}
	}
	sdcomp = (intensityobs - (tmpp*mu1 + (1.0-tmpp)*mu2))/(tmpp*sd1 + (1.0-tmpp)*sd2);
	std_sum = -sdcomp.t()*sdcomp*0.5;
	ninnew = sum(tmpp);
	lLnew = ninnew * (-log(sd1))-ninnew*(-log(sd2)) + std_sum - annew.t()*(annew/eigencnini)*tauinirg*0.5;
	r = Rcpp::dnorm(ankk, anok[0],.05) - Rcpp::dnorm(anok,ankk[0],.05);
	r = r + lLnew - gxZ;
      r(0) = fmin(std::exp(r(0)), 1.0);
	uu = Rcpp::runif(1);
      if(uu(0) <= r(0)) {
		anini(k) = ankk[0];									// update z, nin, and loglikelihood if sample is accepted
		nin = ninnew;
		gxZ = lLnew;
      }

	// end of MH sampling of all elements of z parameter anini for this sweep
	}

	else{
	// reset values of variables used in slice sampling of z
	rndsmp = Rcpp::rexp(1);
	logy(0) = gxZ(0) - rndsmp[0];
	rndsmp = Rcpp::runif(1, 0, 1.0);
	LL(0) = anini(k) - rndsmp[0];
	R(0) = LL(0) + 1.0;
	x1(0)=0.0;
	g(0)=0.0;
	gxF(0)=0.0;
	count(0)=-1.0;
	ninnew = 0; 
	annew = anini; 
	res1 = TRUE;
	res2 = FALSE;
	res3 = FALSE;
	tmpp.fill(0.0);

	// slice sampling of kth element of z parameter anini
	while( res1 ){
            count = count+1.0;
		res1 =  count(0) < 200.0 ;
	      res2 = LL(0)<=-10.0;
            if (res2){ 
			LL(0) = LL(0)+1.0;
			res1= FALSE;}
		else {
			annew(k) = LL(0);
			rmat = tmpmat*annew + mu;
			for(int i = 0; i<s; i++){
				if(demean(i)<rmat(i)){
					tmpp(i) = 1;}
				else{tmpp(i)=0;}
			}			
			sdcomp = (intensityobs - (tmpp*mu1(0) + (1.0-tmpp)*mu2(0)))/(tmpp*sd1(0) + (1.0-tmpp)*sd2(0));
			std_sum = -sdcomp.t()*sdcomp/2.0;
			ninnew = sum(tmpp);
			g = ninnew * (-log(sd1))-ninnew*(-log(sd2)) + std_sum - annew.t()*diagmat(1.0/eigencnini)*annew*tauinirg*0.5;
			res3 = g(0) <= logy(0);
		}
            if (res3){ 
			LL(0) = LL(0)+1.0;
			res1=FALSE;}
            LL(0) = LL(0) - 1.0;
            }
		count = -1.0;
		res1 = TRUE;res2 = FALSE;res3 = FALSE;
       while( res1 ){
            count = count+1.0;
		res1 =  count(0) < 200.0 ;
		res2 = R(0)>=10.0;
            if (res2){ 
			R(0) = R(0)-1.0;
			res1= FALSE;}
		else {
			annew(k) = R(0);
			rmat = tmpmat*annew + mu;
			for(int i = 0; i<s; i++){
				if(demean(i)<rmat(i)){
					tmpp(i) = 1;}else{tmpp(i)=0;}
			}
			sdcomp = (intensityobs - (tmpp*mu1(0) + (1.0-tmpp)*mu2(0)))/(tmpp*sd1(0) + (1.0-tmpp)*sd2(0));
			std_sum = -sdcomp.t()*sdcomp/2.0;
			ninnew = sum(tmpp);
			g = ninnew * (-log(sd1))-ninnew*(-log(sd2)) + std_sum - annew.t()*diagmat(1.0/eigencnini)*annew*tauinirg*0.5;
			res3 = g(0) <= logy(0);
		}
            if (res3){ 
			R(0) = R(0)-1.0;
			res1=FALSE;}
            R(0) = R(0) + 1.0;
            }
		
	res1 = LL(0)<=-10.0;
      if (res1) {
           LL(0) = -10.0;
      }
	res1 = R(0)>10.0;
      if (res1) {
           R(0) = 10.0;
      }
	count = -1.0;
	res1 = TRUE;res2 = FALSE;res3 = FALSE;
      while(res1) {
		count=count+1;
		if (count(0) > 200.0){
		res1 = FALSE;
		gxF = gxZ;
		x1 = anini(k);}
		else{
       	if (res2){ 
			gxF = g;
			res1 = FALSE;
			}
		else {
			x1 = Rcpp::runif(1, LL(0), R(0));
			annew(k) = x1(0);
			rmat = tmpmat*annew + mu;
			for(int i = 0; i<s; i++){
				if(demean(i)<rmat(i)){
					tmpp(i) = 1;}else{tmpp(i)=0;}
			}
			sdcomp = (intensityobs - (tmpp*mu1(0) + (1.0-tmpp)*mu2(0)))/(tmpp*sd1(0) + (1.0-tmpp)*sd2(0));
			std_sum = -sdcomp.t()*sdcomp/2.0;
			ninnew = sum(tmpp);
			g = ninnew * (-log(sd1))-ninnew*(-log(sd2)) + std_sum - annew.t()*diagmat(1.0/eigencnini)*annew*tauinirg*0.5;
			res3 = (x1(0) > anini(k));
		}
	  	  	res2 = (g(0) >= logy(0));
        if (res3) {
            R(0) = x1(0);
        }
        else {
            LL(0) = x1(0);
        }}
      }  
 	anini(k) = x1(0);											// new sample of kth element of z parameter (anini)
	nin = ninnew;											// new value of number of datapoints inside image boundary
	gxZ = gxF;												// new loglikelihood value
	
	} 													// end of slice sampling of all elements of z parameter anini for this sweep
}
astar = alphatau + L/2;
interim = anini.t()*(anini / eigencnini);
bstar = betatau + interim*0.5;
tauinirg = Rcpp::rgamma(1,astar(0),bstar(0));								// sampling the tau parameter
rmat = tmpmat*anini + mu;										// calculating new values of nin and nout for new z (anini) sample.  Also, sampling new values of mu1, mu2, sd1, sd2.
for(int j = 0; j< s; j++){
	if(demean(j)<=rmat(j)){tmpp(j)=1.0;}
	else{tmpp(j)=0.0;}
}
nin = sum(tmpp); nout = s - nin;
xbarin = tmpp.t()*intensityobs/nin; xbarout = intensityobs.t()*(1.0-tmpp)/nout;
ssin = trans(intensityobs % tmpp - xbarin(0)* tmpp)*(intensityobs % tmpp - xbarin(0)* tmpp);
ssout = trans(intensityobs % (1.0-tmpp) - xbarout(0)* (1.0-tmpp))*(intensityobs % (1.0-tmpp) - xbarout(0)* (1.0-tmpp));
sd1 = pow(1/Rcpp::rgamma(1, alpha_2 + nin*0.5, 1/(beta_2+0.5*ssin(0) + (pow(sd0(0), -2)*nin/(nin+pow(sd0(0), -2)))*(pow(xbarin(0) - mu0(0),2.0)*0.5))),0.5);
sd2 = pow(1/Rcpp::rgamma(1, alpha_2 + nout*0.5, 1/(beta_2+0.5*ssout(0) + (pow(sd0(0), -2)*nout/(nout+pow(sd0(0), -2)))*(pow(xbarout(0) - mu0(0),2.0)*0.5))),0.5);
mu1 = Rcpp::rnorm(1,pow(sd0(0), -2)*mu0(0)/(nin+pow(sd0(0), -2)) + nin*xbarin(0)/(nin+pow(sd0(0), -2)), pow(nin+pow(sd0(0), -2),-0.5));
mu2 = Rcpp::rnorm(1,pow(sd0(0), -2)*mu0(0)/(nout+pow(sd0(0), -2)) + nout*xbarout(0)/(nout+pow(sd0(0), -2)),pow(nout+pow(sd0(0), -2),-0.5));
sdcomp = (intensityobs - (tmpp*mu1(0) + (1.0-tmpp)*mu2(0)))/(tmpp*sd1(0) + (1.0-tmpp)*sd2(0));
std_sum = -sdcomp.t()*sdcomp*0.5;
tempmu1(0) = mu1(0);tempmu2(0) = mu2(0);tempsd1(0) = sd1(0);tempsd2(0) = sd2(0);
if(order_mu(0) == I(0)){
mu1 = arma::max(tempmu1, tempmu2);
mu2 = arma::min(tempmu1, tempmu2);
}else if(order_mu(0) == O(0)){
mu1 = arma::min(tempmu1, tempmu2);
mu2 = arma::max(tempmu1, tempmu2);
} else {
mu1=mu1;mu2=mu2;
}
if(order_sd(0) == I(0)){
sd1 = arma::max(tempsd1, tempsd2);
sd2 = arma::min(tempsd1, tempsd2);
}else if(order_sd(0) == O(0)){
sd1 = arma::min(tempsd1, tempsd2);
sd2 = arma::max(tempsd1, tempsd2);
} else {
sd1=sd1;sd2=sd2;
}

gx0 = -0.5 * (sum(log(eigencnini)) + tauinirg * interim) + (alpha_a - 1) * log(lambdaini) - beta_a;
rndsmp = Rcpp::rexp(1);												// resetting values of variables to be used in slice sampling of a parameters (lambdaini) 
logy(0) = gx0(0) - rndsmp[0];
rndsmp = Rcpp::runif(1, 0, 1.0);
LL(0) = lambdaini(0) - rndsmp[0];
R(0) = LL(0) + 1.0;
x1(0)=0.0;
g(0)=0.0;
count(0)=0.0;
interim(0) = 0.0;
calc1(0) = 0.0;
res1 = TRUE;
res2 = FALSE;
res3 = FALSE;

while( res1 ){												// begin slice sampling of a
		interim(0)=0.0;
            count = count+1.0;
		res1 =  count(0) < 100.0 ;
	      res2 = LL(0)<=0.0;
            if (res2){ 
			LL(0) = LL(0)+1.0;
			res1= FALSE;}
		else {
			calc1 = 2*M_PI*std::exp(-2.0*LL(0));
			for(int i=0; i<(i_J+1); i++){
				cnini(i) = calc1(0)*Rcpp::as<double>(besselIs(2*LL,i));
			}
			eigencnini(0) = cnini(0);
			for(int i=1; i<(i_J+1); i++){
				eigencnini(2*i-1) = cnini(i);
				eigencnini(2*i) = cnini(i);
			} 
			interim = anini.t()*(anini / eigencnini);
			g = -0.5 * (sum(log(eigencnini)) + tauinirg * interim) + (alpha_a - 1) * log(LL) - beta_a;
			res3 = g(0) <= logy(0);
		}
            if (res3){ 
			LL(0) = LL(0)+1.0;
			res1=FALSE;}
            LL(0) = LL(0) - 1.0;
            }
		count = 0.0;
		res1 = TRUE;res2 = FALSE;res3 = FALSE;
       while( res1 ){
		interim(0)=0.0;
            count = count+1.0;
		res1 =  count(0) < 100.0 ;
		res2 = R(0)>=500.0;
            if (res2){ 
			R(0) = R(0)-1.0;
			res1= FALSE;}
		else {
			calc1 = 2*M_PI*std::exp(-2.0*R(0));
			for(int i=0; i<(i_J+1); i++){
				cnini(i) = calc1(0)*Rcpp::as<double>(besselIs(2*R,i));
			}
			eigencnini(0) = cnini(0);
			for(int i=1; i<(i_J+1); i++){
				eigencnini(2*i-1) = cnini(i);
				eigencnini(2*i) = cnini(i);
			}
			interim = anini.t()*(anini / eigencnini);
			g = -0.5 * (sum(log(eigencnini)) + tauinirg * interim) + (alpha_a - 1) * log(R) - beta_a;
			res3 = (g(0) <= logy(0));
		}
            if (res3){ 
			R(0) = R(0)-1.0;
			res1=FALSE;}
            R(0) = R(0) + 1.0;
            }
		
	res1 = LL(0)<=0.0;
    if (res1) {
        LL(0) = 0.0;
    }
	res1 = R(0)>500.0;
    if (res1) {
        R(0) = 500.0;
    }
	count = 0.0;
	res1 = TRUE;res2 = FALSE;res3 = FALSE;
    while(res1) {
		interim(0)=0.0;
		count=count+1;
		res1 =  count(0) < 200.0  ;
       	if (res2){ 
			res1 = FALSE;
			}
		else {
			x1 = Rcpp::runif(1, LL(0), R(0));
			calc1 = 2*M_PI*std::exp(-2.0*x1(0));
			for(int i=0; i<(i_J+1); i++){
				cnini(i) = calc1(0)*Rcpp::as<double>(besselIs(2*x1,i));
			}
			eigencnini(0) = cnini(0);
			for(int i=1; i<(i_J+1); i++){
				eigencnini(2*i-1) = cnini(i);
				eigencnini(2*i) = cnini(i);
			}
			interim = anini.t()*(anini / eigencnini);
			g = -0.5 * (sum(log(eigencnini)) + tauinirg * interim) + (alpha_a - 1) * log(x1) - beta_a;
			res3 = (x1(0) > lambdaini(0));
		}
	  	  	res2 = (g(0) >= logy(0));
        if (res3) {
            R(0) = x1(0);
        }
        else {
            LL(0) = x1(0);
        }
    }  

lambdaini = x1;												// end of slice sampling of a and taking on new sample



														// loglikelihood calculation for z using new sample of a
gxZ = nin * (-log(sd1(0)))-nin*(-log(sd2(0))) + std_sum - anini.t()*(anini / eigencnini)*tauinirg*0.5;

if(i > i_nburn){
ansum = ansum + anini;
ansmp.col(i-i_nburn) = anini; 
musigvec(0) = mu1(0); musigvec(1) = mu2(0); musigvec(2) = sd1(0); musigvec(3) = sd2(0);
musigsmp.col(i-i_nburn) = musigvec;
}

}

ansum = ansum/i_nrun;

arma::colvec thetaplot = arma::colvec(200);
for(int k = 0; k<200; k++){
	thetaplot(k) = 2.0*M_PI*(k/200.0); 
}
arma::mat estfunc = arma::mat(200,L);
arma::colvec esttheta = arma::colvec(200);
double est = 0.0;

for(int k = 0; k<200; k++){
		est=0.0;
  		for(int j = 0; j<L; j++){
			estfunc(k,j) = Rcpp::as<double>(eigenfun(j+1,thetaplot(k)));
			est = est + Rcpp::as<double>(eigenfun(j+1,thetaplot(k)))*ansum(j)+(d_inimean/L);
		}
		esttheta(k) =est;
}

arma::mat estthetapts = arma::mat(200,i_nrun);
arma::mat estpts = arma::mat(200,i_nrun);
arma::colvec mutheta = arma::colvec(200);mutheta.fill(d_inimean);

for(int i = 0; i<i_nrun; i++){
		estpts.col(i) = estfunc*ansmp.col(i)+mutheta;
		estthetapts.col(i) = esttheta;
}
arma::colvec variance = arma::colvec(200);

variance = sum(pow(estpts-estthetapts,2),1)/(i_nrun-1.0);
arma::colvec maxed = arma::colvec(1);
arma::colvec diffed = arma::colvec(1);
arma::colvec maxes = arma::colvec(i_nrun);
arma::colvec sortval = arma::colvec(200);
arma::colvec sorted = arma::colvec(200);
for(int j = 0; j<i_nrun; j++){
maxed(0)=0.0;
for(int i = 0; i<200; i++){
diffed = abs(estpts(i,j)-estthetapts(i,j))/pow(variance(i),0.5);
maxed(0) = fmax(maxed(0),diffed(0));
}
maxes(j) = maxed(0);
}
sorted = sort(maxes);
sortval.fill(sorted(floor(0.95*i_nrun)));
arma::colvec lower = arma::colvec(200);
arma::colvec upper = arma::colvec(200);
lower = esttheta - sortval%pow(variance,0.5);
upper = esttheta + sortval%pow(variance,0.5);

if(output_All[0] == TRUE){result = Rcpp::List::create(Rcpp::Named("estimate") = esttheta,Rcpp::Named("theta") = thetaplot,Rcpp::Named("lower") = lower,Rcpp::Named("upper") = upper, Rcpp::Named("musig.smp") = musigsmp, Rcpp::Named("coef.smp") = ansmp);}
else{result = Rcpp::List::create(Rcpp::Named("estimate") = esttheta,Rcpp::Named("theta") = thetaplot,Rcpp::Named("lower") = lower,Rcpp::Named("upper") = upper);}

return result;
}
