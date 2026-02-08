// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

const double xieps = 0.0;

// //' Discrete Extended generalized Pareto distribution of type 1 (deGPD1) negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for each deGPD parameter
// //' @param X1 a design matrix for the deGPD log scale parameter
// //' @param X2 a design matrix for the deGPD log shape parameter
// //' @param X3 a design matrix for the deGPD log kappa
// //' @param yvec a vector
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return degpd1d0 a scalar, the negative log-liklihood
// //' @return degpd1d12 a matrix, first then second derivatives w.r.t. deGPD1 parameters
// //' @return degpd1d34 a matrix, third then fourth derivatives w.r.t. deGPD1 parameters (Not given)
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double degpd1d0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2, const arma::mat& X3, arma::vec yvec, const arma::uvec& dupid, int dcate, const Rcpp::List& offsets)
{
  
  arma::vec lsigmavec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lxivec = X2 * Rcpp::as<arma::vec>(pars[1]);
  arma::vec lkappavec = X3 * Rcpp::as<arma::vec>(pars[2]);
  int nobs = yvec.size();
   
   if (dcate == 1) {
    lsigmavec = lsigmavec.elem(dupid);
    lxivec = lxivec.elem(dupid);
    lkappavec = lkappavec.elem(dupid);
  }

  {
  arma::vec off0 = Rcpp::as<arma::vec>(offsets[0]);
  if (off0.n_elem > 0) lsigmavec += off0;
  arma::vec off1 = Rcpp::as<arma::vec>(offsets[1]);
  if (off1.n_elem > 0) lxivec += off1;
  arma::vec off2 = Rcpp::as<arma::vec>(offsets[2]);
  if (off2.n_elem > 0) lkappavec += off2;
  }

  double y, lsigma, lxi, lkappa;
  double e1,e2,e3, e4, e5;
  double hi; 
  double lo;
  double nllh=0.0;
  
  for (int j=0; j < nobs; j++) {
    
    y = yvec[j];
    lsigma = lsigmavec[j];
    lxi = lxivec[j];
    lkappa = lkappavec[j];
    
    e1=1/exp(lxi);
    e2 =  exp(lxi) / exp(lsigma);
    e3= 1+ (y+1)*e2;
	e4 = 1+  y*e2;
	e5= exp(lkappa);
    hi=  R_pow(1- R_pow(1/e3, e1) , e5);
    lo= R_pow(1-R_pow(1/e4, e1) , e5);
    
    if (e3 <= 0) {
      nllh = 1e20;
      break;
    }
    
    if (e4 <= 0) {
      nllh = 1e20;
      break;
    }
    
    nllh += -log(hi-lo);
    //if (!ISNA(nllh)){
    //nllh = 1e20;
    // break;
//}
 //Rprintf("hi %f ",hi); 
   // Rprintf("lo %f \n",lo);
  }
  return(nllh);
}


// //' @rdname degpd1d0
// [[Rcpp::export]]
arma::mat degpd1d12(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::vec yvec , const arma::uvec dupid, int dcate, const Rcpp::List& offsets)
{
  
  arma::vec lsigmavec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lxivec = X2 * Rcpp::as<arma::vec>(pars[1]);
  arma::vec lkappavec = X3 * Rcpp::as<arma::vec>(pars[2]);
  int nobs = yvec.size();
  arma::mat out = arma::mat(nobs, 9);
  
  if (dcate == 1) {
    lsigmavec = lsigmavec.elem(dupid);
    lxivec = lxivec.elem(dupid);
    lkappavec = lkappavec.elem(dupid);
  }

  {
  arma::vec off0 = Rcpp::as<arma::vec>(offsets[0]);
  if (off0.n_elem > 0) lsigmavec += off0;
  arma::vec off1 = Rcpp::as<arma::vec>(offsets[1]);
  if (off1.n_elem > 0) lxivec += off1;
  arma::vec off2 = Rcpp::as<arma::vec>(offsets[2]);
  if (off2.n_elem > 0) lkappavec += off2;
  }

  double y, lsigma, lxi, lkappa;
  double ee1, ee2, ee3, ee4, ee6, ee8, ee9;
  double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
  double ee20, ee21,ee22, ee23, ee24, ee25, ee26, ee27, ee28, ee29,ee30;
  double ee31, ee32, ee33, ee34, ee35, ee36, ee37, ee38, ee39, ee40;
  double ee42, ee44, ee46, ee47;
  double ee49, ee51, ee52, ee53, ee55, ee56;
  double ee58, ee59, ee60, ee61, ee62;
  double ee63, ee64, ee65, ee66, ee67, ee68;
  double ee69, ee70, ee73, ee74;
  double ee77, ee78, ee81, ee82, ee83;
  double ee84, ee85, ee86, ee87,ee91;
  double ee92, ee94, ee95, ee97, ee100;
  double ee102, ee103, ee105;
  
  double eee1, eee2, eee3, eee5, eee6, eee7, eee8;
  double eee9, eee10;
  double eee11, eee12, eee13, eee14, eee15, eee16;
  double eee17, eee18, eee20, eee22;
  double eee25, eee26, eee29, eee31;
  double eee33, eee34,eee35, eee39, eee41;
  
  
  for (int j=0; j < nobs; j++) {
    
    y  = yvec[j];
      lsigma  = lsigmavec[j];
      lxi  = lxivec[j];
      lkappa  = lkappavec[j];
   
    if(y>0){   
    ee1  = exp(lxi);
    ee2  = exp(lsigma);
    ee3  = 1/ee1;
    ee4  = 1 + y;
    ee6  = (ee4 * ee1)/ee2;
    ee8  = (y * ee1)/ee2;
    ee9  = exp(lkappa);
    ee10  = ee6 + 1;
    ee11  = 1 + ee8;
    ee12  = R_pow(ee10,ee3);
    ee13  = R_pow(ee11,ee3);
    ee14  = 1 - 1/ee12;
    ee15  = 1 - 1/ee13;
    ee16  = 1 + ee3;
    ee17 = ee9 - 1;
    ee18 = R_pow(ee14,ee9);
    ee19 = R_pow(ee15,ee9);
    ee20 = R_pow(ee10,ee16);
    ee21 = R_pow(ee11,ee16);
    ee22 = R_pow(ee14,ee17);
    ee23 = R_pow(ee15,ee17);
    ee24 = ee18 - ee19;
    ee25 = log1p(ee6);
    ee26 = log1p(ee8);
    ee27 = ee3 - 1;
    ee28 = ee20 * ee2;
    ee29 = ee21 * ee2;
    ee30 = 2/ee1;
    ee31 = log(ee14);
    ee32 = log(ee15);
    ee33 = ee12 * ee1;
    ee34 = ee13 * ee1;
    ee35 = ee24 * ee2;
    ee36 = ee22 * ee4;
    ee37 = y * ee23;
    ee38 = R_pow(ee10,ee27);
    ee39 = ee36/ee20;
    ee40 = R_pow(ee11,ee27);
    ee42 = ee9 - 2;
    ee44 = ee25/ee33 - ee4/ee28;
    ee46 = ee26/ee34 - y/ee29;
    ee47 = ee37/ee21;
    ee49 = (ee38 * ee4)/ee2;
    ee51 = (ee12 * ee25)/ee1;
    ee52 = ee18 * ee31;
    ee53 = ee19 * ee32;
    ee55 = (ee13 * ee26)/ee1;
    ee56 = ee47 - ee39;
    ee58 = (y * ee40)/ee2;
    ee59 = ee49 - ee51;
    ee60 = R_pow(ee10,ee30);
    ee61 = ee52 - ee53;
    ee62 = R_pow(ee11,ee30);
    ee63 = 2 * ee16;
    ee64 = ee58 - ee55;
    ee65 = ee22 * ee44;
    ee66 = ee23 * ee46;
    ee67 = R_pow(ee14,ee42);
    ee68 = ee66 - ee65;
    ee69 = R_pow(ee15,ee42);
    ee70 = R_pow(ee35,2);
    ee73 = (ee59 * ee22)/ee60 - (ee23 * ee64)/ee62;
    ee74 = ee61 * ee9;
    ee77 = ee12 * ee16 * ee4 * ee1; 
    ee78 = ee69 * ee17;
    ee81 = y * ee16 * ee13 * ee1;
    ee82 = R_pow(ee28,2);
    ee83 = ee59 * ee67;
    ee84 = R_pow(ee33,2);
    ee85 = ee68 * ee9;
    ee86 = R_pow(ee29,2);
    ee87 = R_pow(ee34,2);
    ee91 = ee77/ee2 - (ee20 * ee25)/ee1;
    ee92 = R_pow(ee10,ee63);
    ee94 = ee22 * ee9 * ee31;
    ee95 = ee67 * ee17;
    ee97 = ee23 * ee9 * ee32;
    ee100 = R_pow(ee11,ee63);
    //ee101 = 1 + 3/ee1;
    ee102 = ee3 - ee63;
    ee103 = ee9 * ee56;
    ee105 = ee81/ee2 - ee21 * ee26/ee1; 

out(j, 0) = -(ee103/ee35);
    out(j, 1) = -(ee85/ee24); 
    out(j, 2) = -(ee74/ee24);
    out(j, 3) = -(((R_pow(y,2) * 
                      (ee23 * ee16 * R_pow(ee11,ee102) * ee1 - ee78/ee100) - (R_pow(ee10,ee102) * 
                                                                         ee22 * ee16 * ee1 - ee95/ee92) * R_pow(ee4,2))/(ee24 * R_pow(ee2,2)) - 
                     (ee35 + ee103) * ee56/ee70) * ee9);
    out(j, 4) = -(ee9 * (y * 
                           (((ee29 - ee81)/ee86 + (ee40 * ee1 * ee26/ee87 - 1/ee21)/ee2) * 
                              ee23 - ee78 * ee46/ee29) - ((((ee28 - ee77)/ee82 + 
                                                              (ee38 * ee1 * ee25/ee84 - 1/ee20)/ee2) * ee22 - ee95 * 
                                                             ee44/ee28) * ee4 + ee85 * ee56/ee35))/ee24);
    out(j, 5) = -(ee9 * 
                    (y * (ee97/ee21 + ee23/ee21) - ((ee94/ee20 + ee22/ee20) * 
                                                      ee4 + ee74 * ee56/ee24))/ee35); 
    out(j, 6) = -(((ee69 * ee46 * 
                      ee64/ee62 - ee83 * ee44/ee60) * ee17 + ee23 * (y * (1/ee29 + 
                                                                            ee2 * ee105/ee86) - (ee13 + ee58 - ee55) * ee1 * ee26/ee87) - 
                     (((ee91 * ee2/ee82 + 1/ee28) * ee4 - (ee49 + ee12 - ee51) * 
                         ee1 * ee25/ee84) * ee22 + ee73 * ee68 * ee9/ee24)) * ee9/ee24);
    out(j, 7) =  -((ee59 * (ee94/ee60 + ee22/ee60) - 
                      (ee73 * ee61 * ee9/ee24 + (ee97/ee62 + ee23/ee62) * ee64)) * ee9/ee24);
    out(j, 8) = -(((ee18 * R_pow(ee31,2) - (R_pow(ee61,2)/ee24 + 
                                       ee19 * R_pow(ee32,2)))*ee9 + ee52 - ee53) * ee9/ee24);
  }

else {
     	   
    eee1 = exp(lxi);
    eee2 = exp(lsigma);
    eee3 = 1 + y;
    eee5 = eee3 * eee1/eee2;
    eee6 = 1/eee1;
    eee7 = eee5 + 1;
    eee8 = R_pow(eee7,eee6);
    eee9 = 1 + eee6;
    eee10 = R_pow(eee7,eee9);
    eee11 = 1 - 1/eee8;
    eee12 = exp(lkappa);
    eee13 = log1p(eee5);
    eee14 = eee10 * eee2;
    eee15 = eee10 * eee11;
    eee16 = eee15 * eee2;
    eee17 = eee8 * eee1;
    eee18 = R_pow(eee7,eee6 - 1);
    eee20 = eee18 * eee3/eee2;
    eee22 = eee8 * eee13/eee1;
    eee25 = eee13/eee17 - eee3/eee14;
    eee26 = eee20 - eee22;
    eee29 = eee8 * eee9 * eee3 * eee1;
    eee31 = eee3 * eee12/eee16;
    eee33 = R_pow(eee16,2);
    eee34 = R_pow(eee14,2);
    eee35 = R_pow(eee17,2);
    eee39 = eee29/eee2 - eee10 * eee13/eee1;
    eee41 = R_pow(eee7,2/eee1) * eee11;
    
	
	out(j, 0) = eee31; 
    out(j, 1) = eee12 * eee25/eee11;
    out(j, 2) = -(eee12 *log(eee11)); 
    out(j, 3) = -((eee16 - (eee8*eee11 * eee9 * eee1 + 1) * eee3) * eee3 * eee12/eee33); 
    out(j, 4) = ((eee14 - 
                    eee29)/eee34 + (eee18 * eee1 * eee13/eee35 + eee25/eee15 - 1/eee10)/eee2) * eee3 * eee12/eee11; 
    out(j, 5) = eee31;  
    out(j, 6) = ((eee39 * eee2/eee34 + 1/eee14) * eee3 - (eee26 * eee25/eee41 + 
                                                       (eee20 + eee8 - eee22) * eee1 * eee13/eee35)) * eee12/eee11;
    out(j, 7) =  -(eee26 * eee12/eee41); 
    out(j, 8) =  -(eee12 *log(eee11));                                                                            
  }
}
     
   return out;
}
  
  

// //' Discrete Extended generalized Pareto distribution of type 2 (DeGPD2) negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for each DeGPD parameter
// //' @param X1 a design matrix for the eGPD log scale parameter
// //' @param X2 a design matrix for the eGPD log shape parameter
// //' @param X3 a design matrix for the eGPD log kappa1 parameter
// //' @param X4 a design matrix for the eGPD log kappa2 parameter
// //' @param X5 a design matrix for the eGPD logit p parameter
// //' @param yvec a vector
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return degpd2d0 a scalar, the negative log-liklihood
// //' @return degpd2d12 a matrix, first then second derivatives w.r.t. deGPD2 parameters
// //' @return degpd2d34 a matrix, third then fourth derivatives w.r.t. deGPDD2 parameters (Not available)
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]

double degpd2d0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2, const arma::mat& X3, const arma::mat& X4, const arma::mat& X5, arma::vec yvec, const arma::uvec& dupid, int dcate, const Rcpp::List& offsets)
{
  arma::vec lsigmavec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lxivec = X2 * Rcpp::as<arma::vec>(pars[1]);
  arma::vec lkappa1vec = X3 * Rcpp::as<arma::vec>(pars[2]);
  arma::vec ldkappavec = X4 * Rcpp::as<arma::vec>(pars[3]);
  arma::vec logitpvec = X5 * Rcpp::as<arma::vec>(pars[4]);
  int nobs = yvec.size();

   if (dcate == 1) {
  lsigmavec = lsigmavec.elem(dupid);
  lxivec = lxivec.elem(dupid);
  lkappa1vec = lkappa1vec.elem(dupid);
  ldkappavec = ldkappavec.elem(dupid);
  logitpvec = logitpvec.elem(dupid);
}

  {
  arma::vec off0 = Rcpp::as<arma::vec>(offsets[0]);
  if (off0.n_elem > 0) lsigmavec += off0;
  arma::vec off1 = Rcpp::as<arma::vec>(offsets[1]);
  if (off1.n_elem > 0) lxivec += off1;
  arma::vec off2 = Rcpp::as<arma::vec>(offsets[2]);
  if (off2.n_elem > 0) lkappa1vec += off2;
  arma::vec off3 = Rcpp::as<arma::vec>(offsets[3]);
  if (off3.n_elem > 0) ldkappavec += off3;
  arma::vec off4 = Rcpp::as<arma::vec>(offsets[4]);
  if (off4.n_elem > 0) logitpvec += off4;
  }

  double y, lsigma, lxi, lkappa1, lkappa2, logitp;
  double e1,e2,e3, e4,e5,e6,e7, e8, e9,e10,e11, e12;
  double lo, hi;
  double nllh=0.0;

  for (int j=0; j < nobs; j++) {

    y = yvec[j];
    lsigma = lsigmavec[j];
    lxi = lxivec[j];
    lkappa1 = lkappa1vec[j];
    // Reparameterization: ldkappa = log(kappa2 - kappa1), compute lkappa2 via log-sum-exp
    {
      double ldk = ldkappavec[j];
      double mx = std::max(lkappa1, ldk);
      lkappa2 = mx + log1p(exp(-std::abs(lkappa1 - ldk)));
    }
    logitp = logitpvec[j];
    
    e1=1.0/exp(lxi);
    e2 =  exp(lxi) / exp(lsigma);
    e3= 1+ (y+1)*e2;
    e4= R_pow(1/e3, e1);
    e5= R_pow(1-e4, exp(lkappa1)); //(H(y+1))^lkappa1
    e6= 1+ y*e2;
    e7= R_pow(1/e6, e1);
    e8= R_pow(1-e7, exp(lkappa1)); //(H(y))^lkappa1
    e9= exp(logitp)/ (1+exp(logitp));
	hi= e9*(e5-e8);
	e10= R_pow(1-e4, exp(lkappa2)); //(H(y+1))^lkappa2
	e11= R_pow(1-e7, exp(lkappa2)); //(H(y))^lkappa2
	e12= 1/(1+exp(logitp));
	lo= e12*(e10-e11);
    
    nllh += -log(hi+lo);
    //if (!ISNA(nllh)){
    //nllh = 1e20;
    // break;
//}
// Rprintf("hi %f ",hi); 
  //  Rprintf("lo %f \n",lo);
 }
  return(nllh);
}


// //' @rdname degpd2d0
// [[Rcpp::export]]
arma::mat degpd2d12(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::mat X3,  arma::mat X4, arma::mat X5, arma::vec yvec, const arma::uvec dupid, int dcate, const Rcpp::List& offsets)
{
  
  arma::vec lsigmavec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lxivec = X2 * Rcpp::as<arma::vec>(pars[1]);
  arma::vec lkappa1vec = X3 * Rcpp::as<arma::vec>(pars[2]);
  arma::vec ldkappavec = X4 * Rcpp::as<arma::vec>(pars[3]);
  arma::vec logitpvec = X5 * Rcpp::as<arma::vec>(pars[4]);
  int nobs = yvec.size();
  arma::mat out = arma::mat(nobs, 20);

 if (dcate == 1) {
  lsigmavec = lsigmavec.elem(dupid);
  lxivec = lxivec.elem(dupid);
  lkappa1vec = lkappa1vec.elem(dupid);
  ldkappavec = ldkappavec.elem(dupid);
  logitpvec = logitpvec.elem(dupid);
}

  {
  arma::vec off0 = Rcpp::as<arma::vec>(offsets[0]);
  if (off0.n_elem > 0) lsigmavec += off0;
  arma::vec off1 = Rcpp::as<arma::vec>(offsets[1]);
  if (off1.n_elem > 0) lxivec += off1;
  arma::vec off2 = Rcpp::as<arma::vec>(offsets[2]);
  if (off2.n_elem > 0) lkappa1vec += off2;
  arma::vec off3 = Rcpp::as<arma::vec>(offsets[3]);
  if (off3.n_elem > 0) ldkappavec += off3;
  arma::vec off4 = Rcpp::as<arma::vec>(offsets[4]);
  if (off4.n_elem > 0) logitpvec += off4;
  }

  double y, lsigma, lxi, lkappa1, lkappa2, logitp;
  double ee1, ee2, ee3, ee4,  ee6, ee8,ee9, ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
  double ee20, ee21,ee22, ee23, ee24, ee25, ee26,ee27, ee28, ee30, ee31,  ee32, ee33, ee34, ee35, ee36, ee37, ee38, ee39, ee40, ee41, ee42,ee43, ee44, ee45, ee46, ee47;
  double ee50,ee52, ee54, ee56,ee58, ee60, ee61, ee62, ee63;
  double ee64,ee65, ee66, ee67,ee68,ee69,ee70,ee71, ee72,ee73, ee74,ee75, ee76, ee77, ee78,ee79,ee80,ee81, ee82,ee83,ee84,ee85,ee86;
  double ee87,ee88,ee89,ee90,ee91, ee92,ee93,ee94,ee95, ee96,ee99, ee100;
  double ee102, ee105,ee108, ee109, ee112, ee113,ee114,ee115,ee116, ee120,ee121,ee122,ee123, ee124;
  double ee125,ee128,ee130,ee132, ee134;
  double ee135,ee136, ee137,ee140,ee142,ee143,ee144 ,ee149,ee150,ee151,ee153,ee158,ee160,ee161,ee163,ee164,ee166;
  double ee168,ee170,ee171,ee172,ee175,ee176,ee178,ee179;
  
  double eee1, eee2, eee3, eee5, eee6, eee7, eee8, eee9, eee10, eee11, eee12, eee13, eee14,eee15;
  double eee17, eee18,eee19, eee20,eee21,  eee22, eee23, eee24, eee25, eee26,eee27, eee28, eee29, eee30, eee31;
  double eee32,eee34,eee35, eee36, eee37,eee38,eee40, eee42,eee43, eee44,eee45, eee48, eee49,eee50,eee51,eee52,  eee54;
  double  eee57,eee58, eee59, eee60,eee61, eee64,eee68, eee69,eee70,eee71,eee75,  eee76;
  double eee77,eee80, eee82, eee84,eee88,eee93, eee95;
  
  
  
  for (int j=0; j < nobs; j++) {
    
    y  = yvec[j];
      lsigma  = lsigmavec[j];
      lxi  = lxivec[j];
      lkappa1  = lkappa1vec[j];
      // Reparameterization: ldkappa = log(kappa2 - kappa1), compute lkappa2 via log-sum-exp
      {
        double ldk = ldkappavec[j];
        double mx = std::max(lkappa1, ldk);
        lkappa2 = mx + log1p(exp(-std::abs(lkappa1 - ldk)));
      }
      logitp  = logitpvec[j];

    if(y>0){
      ee1 = exp(lxi);
    ee2 = exp(lsigma);
    ee3 = 1/ee1;
    ee4 = 1 + y;
    ee6 = ee4 * ee1/ee2;
    ee8 = y * ee1/ee2;
    ee9 = ee6 + 1;
    ee10 = 1 + ee8;
    ee11 = R_pow(ee9,ee3);
    ee12 = R_pow(ee10,ee3);
    ee13 = exp(lkappa1);
    ee14 = exp(lkappa2);
    ee15 = 1 - 1/ee11;
    ee16 = 1 - 1/ee12;
    ee17 = 1 + ee3;
    ee18 = exp(logitp);
    ee19 = ee13 - 1;
    ee20 = ee14 - 1;
    ee21 = R_pow(ee15,ee13);
    ee22 = R_pow(ee16,ee13);
    ee23 = R_pow(ee15,ee14);
    ee24 = R_pow(ee16,ee14);
    ee25 = R_pow(ee9,ee17);
    ee26 = R_pow(ee10,ee17);
    ee27 = ee21 - ee22;
    ee28 = ee27 * ee18;
    ee30 = ee28 + ee23 - ee24;
    ee31 = log1p(ee6);
    ee32 = log1p(ee8);
    ee33 = ee3 - 1;
    ee34 = R_pow(ee15,ee19);
    ee35 = R_pow(ee16,ee19);
    ee36 = R_pow(ee15,ee20);
    ee37 = R_pow(ee16,ee20);
    ee38 = 2/ee1;
    ee39 = log(ee15);
    ee40 = log(ee16);
    ee41 = ee25 * ee2;
    ee42 = ee26 * ee2;
    ee43 = ee11 * ee1;
    ee44 = ee12 * ee1;
    ee45 = R_pow(ee9,ee33);
    ee46 = R_pow(ee10,ee33);
    ee47 = 1 + ee18;
    ee50 = ee31/ee43 - ee4/ee41;
    ee52 = ee32/ee44 - y/ee42;
    ee54 = ee45 * ee4/ee2;
    ee56 = ee11 * ee31/ee1;
    ee58 = ee12 * ee32/ee1;
    ee60 = y * ee46/ee2;
    ee61 = ee54 - ee56;
    ee62 = R_pow(ee9,ee38);
    ee63 = R_pow(ee10,ee38);
    ee64 = ee60 - ee58;
    ee65 = ee30 * ee2;
    ee66 = ee34 * ee4;
    ee67 = ee13 * ee18;
    ee68 = y * ee35;
    ee69 = ee66/ee25;
    ee70 = ee36 * ee4;
    ee71 = ee68/ee26;
    ee72 = y * ee37;
    ee73 = ee70/ee25;
    ee74 = ee71 - ee69;
    ee75 = ee72/ee26;
    ee76 = 2 * ee17;
    ee77 = ee14 * (ee75 - ee73);
    ee78 = ee21 * ee39;
    ee79 = ee23 * ee39;
    ee80 = ee22 * ee40;
    ee81 = ee24 * ee40;
    ee82 = ee67 * ee74;
    ee83 = ee78 - ee80;
    ee84 = ee79 - ee81;
    ee85 = ee82 + ee77;
    ee86 = 1 - ee18/ee47;
    ee87 = ee34 * ee50;
    ee88 = ee35 * ee52;
    ee89 = ee13 - 2;
    ee90 = ee14 - 2;
    ee91 = ee27 * ee86;
    ee92 = (ee23 - ee24)/ee47;
    ee93 = (ee88 - ee87) * ee13;
    ee94 = ee36 * ee50;
    ee95 = ee37 * ee52;
    ee96 = (ee61 * ee36/ee62 - ee37 * ee64/ee63) * ee14;
    ee99 = ee61 * ee34/ee62 - ee35 * ee64/ee63;
    ee100 = ee91 - ee92;
    ee102 = ee93 * ee18 + (ee95 - ee94) * ee14;
    ee105 = ee11 * ee17 * ee4 * ee1;
    ee108 = y * ee17 * ee12 * ee1;
    ee109 = R_pow(ee65,2);
    ee112 = ee99 * ee13 * ee18 + ee96;
    ee113 = R_pow(ee41,2);
    ee114 = R_pow(ee43,2);
    ee115 = R_pow(ee42,2);
    ee116 = R_pow(ee44,2);
    ee120 = ee105/ee2 - ee25 * ee31/ee1;
    ee121 = R_pow(ee9,ee76);
    ee122 = R_pow(ee15,ee89);
    ee123 = R_pow(ee15,ee90);
    ee124 = R_pow(ee16,ee89);
    ee125 = R_pow(ee16,ee90);
    ee128 = R_pow(ee10,ee76);
    //ee129 = 1 + 3/ee1;
    ee130 = ee3 - ee76;
    ee132 = ee108/ee2 - ee26 * ee32/ee1;
    ee134 = ee100/ee30 + 1/ee47;
    ee135 = R_pow(ee30,2);
    ee136 = ee124 * ee19;
    ee137 = ee125 * ee20;
    ee140 = (ee120 * ee2/ee113 + 1/ee41) * ee4 - (ee54 + ee11 - 
                                                     ee56) * ee1 * ee31/ee114;
    ee142 = (ee41 - ee105)/ee113 + (ee45 * ee1 * ee31/ee114 - 
                                       1/ee25)/ee2;
    ee143 = ee61 * ee122;
    ee144 = ee61 * ee123;
    ee149 = ee83 * ee84 * ee13 * ee14 * ee18/ee135;
    ee150 = ee83 * ee85;
    ee151 = ee84 * ee85;
    ee153 = (ee42 - ee108)/ee115 + (ee46 * ee1 * ee32/ee116 - 
                                       1/ee26)/ee2;
    //ee157 = R_pow(ee9,ee129);
    ee158 = R_pow(ee9,ee130);
    ee160 = ee34 * ee13 * ee39;
    ee161 = ee122 * ee19;
    ee163 = ee36 * ee14 * ee39;
    ee164 = ee123 * ee20;
    ee166 = ee35 * ee13 * ee40;
    ee168 = ee37 * ee14 * ee40;
    //ee169 = R_pow(ee10,ee129);
    ee170 = R_pow(ee10,ee130);
    ee171 = R_pow(ee4,2);
    ee172 = 1 - ee134 * ee18;
    ee175 = R_pow(ee39,2);
    ee176 = R_pow(ee40,2);
    ee178 = y * (1/ee42 + ee2 * ee132/ee115) - (ee12 + ee60 - 
                                                   ee58) * ee1 * ee32/ee116;
    ee179 = R_pow(y,2);

    out(j, 0) = -(ee85/ee65); //  #w.r.t to lsigma
    out(j, 1) = -(ee102/ee30); // #w.r.t to lxi
    out(j, 2) = -(ee83 * ee13 * ee18/ee30); // #w.r.t to lkappa1
    out(j, 3) = -(ee84 * 
                  ee14/ee30); // #w.r.t to lkappa2
	out(j, 4) = -(ee100 * ee18/ee30); // #w.r.t to logitp
	
    out(j, 5) = -((ee67 * 
                    (ee179 * (ee35 * ee17 * ee170 * ee1 - ee136/ee128) - 
                       (ee158 * ee34 * ee17 * ee1 - ee161/ee121) * ee171) + 
                    ee14 * (ee179 * (ee37 * ee17 * ee170 * ee1 - ee137/ee128) - 
                              (ee158 * ee36 * ee17 * ee1 - ee164/ee121) * ee171))/(ee30 * 
                                                                                     R_pow(ee2,2)) - (ee65 + ee82 + ee77) * ee85/ee109); // #w.r.t to lsigma,lsigma
    out(j, 6) = -((ee67 * 
                   (y * (ee153 * ee35 - ee136 * ee52/ee42) - (ee142 * ee34 - 
                                                                ee161 * ee50/ee41) * ee4) + ee14 * (y * (ee153 * 
                                                                                                           ee37 - ee137 * ee52/ee42) - (ee142 * ee36 - ee164 * ee50/ee41) * 
                                                                                                      ee4) - ee102 * ee85/ee65)/ee30);// #w.r.t to lsigma,lxi
    out(j, 7) = -(ee67 * (y * 
                          (ee166/ee26 + ee35/ee26) - ((ee160/ee25 + ee34/ee25) * 
                                                        ee4 + ee150/ee30))/ee65);  //#w.r.t to lsigma, lkappa1
    out(j, 8) = -(ee14 * (y * (ee168/ee26 + 
                               ee37/ee26) - ((ee163/ee25 + ee36/ee25) * ee4 + ee151/ee30))/ee65); // #w.r.t to lsigma,lkappa2
                               
                               
out(j, 9) =  -((ee86 * ee13 * ee74 - (ee100 * ee85/ee30 + 
                                          ee77/ee47)) * ee18/ee65); // #w.r.t to lsigma,logitp
out(j, 10) = -((((ee124 * ee52 * ee64/ee63 - 
                     ee143 * ee50/ee62) * ee19 + ee35 * ee178 - ee140 * ee34) * 
                   ee13 * ee18 + ((ee125 * ee52 * ee64/ee63 - ee144 * ee50/ee62) * 
                                    ee20 + ee37 * ee178 - ee140 * ee36) * ee14 - ee112 * 
                   ee102/ee30)/ee30); //  #w.r.t to lxi,lxi
out(j, 11) = -((ee61 * (ee160/ee62 + 
                            ee34/ee62) - (ee112 * ee83/ee30 + (ee166/ee63 + ee35/ee63) * 
                                            ee64)) * ee13 * ee18/ee30); //  #w.r.t to lxi,lkappa1
out(j, 12) = -((ee61 * (ee163/ee62 + 
                            ee36/ee62) - (ee112 * ee84/ee30 + (ee168/ee63 + ee37/ee63) * 
                                            ee64)) * ee14/ee30); //  #w.r.t to lxi,lkappa2
out(j, 13) = -((ee99 * ee86 * ee13 - 
                    (ee112 * ee100/ee30 + ee96/ee47)) * ee18/ee30); // #w.r.t to lxi,logitp

out(j, 14) = -(((ee21 * ee175 - 
                     (R_pow(ee83,2)* ee18/ee30 + ee22 * ee176)) * ee13 + ee78 - 
                    ee80) * ee13 * ee18/ee30); //  #w.r.t to lkappa1,lkappa1




out(j, 15) = ee149; //  #w.r.t to lkappa1,lkappa2
out(j, 16) = -(ee83 * 
                   ee172 * ee13 * ee18/ee30); // #w.r.t to lkappa1,logitp
    
out(j, 17) = -(((ee23 * 
                     ee175 - (R_pow(ee84,2)/ee30 + ee24 * ee176)) * ee14 + ee79 - 
                    ee81) * ee14/ee30); //  #w.r.t to lkappa2,lkappa2
out(j, 18) =  ee134 * ee84 * ee14 * ee18/ee30; // #lkappa2,logitp
out(j, 19) = -((((ee92 - ee91) * ee18 + ee24 -
                      ee23)/ee47 + ee27 * ee172) * ee18/ee30); //#w.r.t to logitp,logitp

// Chain rule: transform derivatives from (lkappa1, lkappa2) to (lkappa1, ldkappa)
{
  double r1 = exp(lkappa1 - lkappa2);
  double r2 = 1.0 - r1;
  double g3 = out(j, 2);  // grad lkappa1
  double g4 = out(j, 3);  // grad lkappa2

  // Gradient transformation
  out(j, 2) = g3 + g4 * r1;
  out(j, 3) = g4 * r2;

  // Hessian cross terms with param 0 (lsigma): H(0,3)=col7, H(0,4)=col8
  double H03 = out(j, 7), H04 = out(j, 8);
  out(j, 7) = H03 + H04 * r1;
  out(j, 8) = H04 * r2;

  // Cross terms with param 1 (lxi): H(1,3)=col11, H(1,4)=col12
  double H13 = out(j, 11), H14 = out(j, 12);
  out(j, 11) = H13 + H14 * r1;
  out(j, 12) = H14 * r2;

  // Diagonal block: H(3,3)=col14, H(3,4)=col15, H(4,4)=col17
  double H33 = out(j, 14), H34 = out(j, 15), H44 = out(j, 17);
  out(j, 14) = H33 + 2.0 * r1 * H34 + r1 * r1 * H44 + g4 * r1 * r2;
  out(j, 15) = r2 * H34 + r1 * r2 * H44 - g4 * r1 * r2;
  out(j, 17) = r2 * r2 * H44 + g4 * r1 * r2;

  // Cross terms with param 4 (logitp): H(4,3)=col16, H(4,4)=col18
  double H43 = out(j, 16), H44p = out(j, 18);
  out(j, 16) = H43 + H44p * r1;
  out(j, 18) = H44p * r2;
}

  }

else {
     	   
    eee1 = exp(lxi);
    eee2 = exp(lsigma);
    eee3 = 1 + y;
    eee5 = eee3 * eee1/eee2;
    eee6 = eee5 + 1;
    eee7 = 1/eee1;
    eee8 = R_pow(eee6,eee7);
    eee9 = 1 - 1/eee8;
    eee10 = exp(lkappa1);
    eee11 = exp(lkappa2);
    eee12 = exp(logitp);
    eee13 = R_pow(eee9,eee10);
    eee14 = R_pow(eee9,eee11);
    eee15 = 1 + eee7;
    eee17 = eee13 * eee12 + eee14;
    eee18 = R_pow(eee6,eee15);
    eee19 = eee10 - 1;
    eee20 = eee11 - 1;
    eee21 = R_pow(eee9,eee19);
    eee22 = R_pow(eee9,eee20);
    eee23 = eee21 * eee10;
    eee24 = eee22 * eee11;
    eee25 = eee23 * eee12;
    eee26 = 1 + eee12;
    eee27 = log(eee9);
    eee28 = log1p(eee5);
    eee29 = R_pow(eee6,(2/eee1));
    eee30 = eee24/eee18;
    eee31 = eee17 * eee2;
    eee32 = eee18 * eee2;
    eee34 = eee25/eee18 + eee30;
    eee35 = 1 - eee12/eee26;
    eee36 = R_pow(eee6,(eee7 - 1));
    eee37 = eee8 * eee1;
    eee38 = eee13 * eee35;
    eee40 = eee36 * eee3/eee2;
    eee42 = eee8 * eee28/eee1;
    eee43 = eee14/eee26;
    eee44 = eee40 - eee42;
    eee45 = eee25 + eee24;
    eee48 = eee28/eee37 - eee3/eee32;
    eee49 = eee24/eee29;
    eee50 = eee38 - eee43;
    eee51 = 2 * eee15;
    eee52 = R_pow(eee31,2);
    eee54 = eee25/eee29 + eee49;
    eee57 = eee8 * eee15 * eee3 * eee1;
    eee58 = R_pow(eee6,eee51);
    eee59 = R_pow(eee9,(eee10 - 2));
    eee60 = R_pow(eee9,(eee11 - 2));
    eee61 = eee34 * eee13;
    //eee63 = eee17 * eee18 * eee2;
    eee64 = R_pow(eee17,2);
    eee68 = eee57/eee2 - eee18 * eee28/eee1;
    eee69 = R_pow(eee9,(eee10 + eee11));
    eee70 = R_pow(eee32, 2);
    eee71 = R_pow(eee37, 2);
    //eee74 = eee61 * eee2;
    eee75 = eee34 * eee14;
    eee76 = eee34 * eee3;
    eee77 = (eee13 - R_pow(eee9, (2 * eee10)) * eee12/eee17) * eee10;
    eee80 = eee50 * eee13 * eee12/eee17;
   // eee81 = R_pow(eee6, (1 + 3/eee1));
    eee82 = R_pow(eee6, (eee7 - eee51));
    eee84 = eee21 * eee35 * eee10;
    //eee85 = eee23 - eee45 * eee13/eee17;
    eee88 = eee59 * eee10 * eee19 * eee12;
    eee93 = eee69 * eee10 * eee11 * eee12 * R_pow(eee27,2)/eee64;
    eee95 = eee60 * eee11 * eee20;


    
	
	out(j, 0) = eee76/eee31; // # w.r.t. lsigma
    out(j, 1) = eee45 * eee48/eee17; // # w.r.t. lxi
    out(j, 2) = -(eee13 * eee10 * eee12 * eee27/eee17); // # w.r.t. lkappa1
    out(j, 3) = -(eee14 * 
                    eee11 * eee27/eee17); // # w.r.t. lkappa2
    out(j, 4) =  -(eee50 * eee12/eee17); //  # w.r.t. logitp
    
    out(j, 5) = (((eee82 * eee21 * eee15 * eee1 - 
                      eee59 * eee19/eee58) * eee10 * eee12 + (eee82 * eee22 * 
                                                           eee15 * eee1 - eee60 * eee20/eee58) * eee11) * eee3/(eee17 * 
                                                                                                           R_pow(eee2,2)) - (eee31 - eee76) * eee34/eee52) * eee3;  // # w.r.t (lsigma, lsigma)
    out(j, 6) = (((eee34 * 
                     eee48/eee17 + eee36 * eee1 * eee28/eee71 - 1/eee18)/eee2 + 
                    (eee32 - eee57)/eee70) * eee45 - (eee88/eee18 + eee95/eee18) * 
                   eee48/eee2) * eee3/eee17 ; // # weereet (lsigma, lxi)
    out(j, 7) =    -(((eee61/eee17 - eee23/eee18) * 
         eee27 - eee21/eee18) * eee3 * eee10 * eee12/eee31); //# w.r.t (lsigma, lkappa1)
    out(j, 8) = -(((eee75/eee17 - 
                      eee30) * eee27 - eee22/eee18) * eee3 * eee11/eee31); // # w.r.t (lsigma, lkappa2)  
	out(j, 9) = -((eee34 * 
                       eee50/eee17 + eee24/(eee18 * eee26) - eee84/eee18) * eee3 * 
                      eee12/eee31); //# w.r.t (lsigma, logitp)
	out(j, 10) = (((eee68 * eee2/eee70 + 
                      1/eee32) * eee3 - (eee40 + eee8 - eee42) * eee1 * eee28/eee71) * 
                    eee45 + eee44 * (eee88/eee29 + eee95/eee29 - eee45 * eee54/eee17) * 
                    eee48)/eee17; // # w.r.t (lxi, lxi)
	out(j, 11) =  -(((eee23/eee29 - eee54 * eee13/eee17) * 
                      eee27 + eee21/eee29) * eee44 * eee10 * eee12/eee17); // # w.r.t (lxi, lkappa1)
	out(j, 12) = -(((eee49 - 
                        eee54 * eee14/eee17) * eee27 + eee22/eee29) * eee44 * eee11/eee17); //# w.r.t (lxi, lkappa2)
	out(j, 13) = -(eee44 * (eee84/eee29 - (eee54 * eee50/eee17 + 
                                          eee24/(eee29 * eee26))) * eee12/eee17); // # w.r.t (lxi, logitp)
    
    out(j, 14) = -((eee77 * eee27 + eee13) * eee10 * 
                     eee12 * eee27/eee17); // # w.r.t (lkappa1, lkappa1)
    
    out(j, 15) = eee93; //  # w.r.t (lkappa1, lkappa2)
    
    out(j, 16) = -((eee38 - 
                       eee80) * eee10 * eee12 * eee27/eee17); //# w.r.t (lkappa2, logitp)
    
    out(j, 17) = -(((eee14 - 
                       R_pow(eee9,(2 * eee11))/eee17) * eee11 * eee27 + eee14) * eee11 * 
                     eee27/eee17); //  # w.r.t (lkappa2, lkappa2)
    
    out(j, 18) = (eee50 * eee14/eee17 + eee43) * eee11 * eee12 * eee27/eee17; //# w.r.t (lkappa2, logitp)
    
    out(j, 19) = -((((eee43 -
                         eee38) * eee12 - eee14)/eee26 + eee38 - eee80) * eee12/eee17); //# w.r.t (logitp, logitp)

// Chain rule: transform derivatives from (lkappa1, lkappa2) to (lkappa1, ldkappa)
{
  double r1 = exp(lkappa1 - lkappa2);
  double r2 = 1.0 - r1;
  double g3 = out(j, 2);  // grad lkappa1
  double g4 = out(j, 3);  // grad lkappa2

  // Gradient transformation
  out(j, 2) = g3 + g4 * r1;
  out(j, 3) = g4 * r2;

  // Hessian cross terms with param 0 (lsigma): H(0,3)=col7, H(0,4)=col8
  double H03 = out(j, 7), H04 = out(j, 8);
  out(j, 7) = H03 + H04 * r1;
  out(j, 8) = H04 * r2;

  // Cross terms with param 1 (lxi): H(1,3)=col11, H(1,4)=col12
  double H13 = out(j, 11), H14 = out(j, 12);
  out(j, 11) = H13 + H14 * r1;
  out(j, 12) = H14 * r2;

  // Diagonal block: H(3,3)=col14, H(3,4)=col15, H(4,4)=col17
  double H33 = out(j, 14), H34 = out(j, 15), H44 = out(j, 17);
  out(j, 14) = H33 + 2.0 * r1 * H34 + r1 * r1 * H44 + g4 * r1 * r2;
  out(j, 15) = r2 * H34 + r1 * r2 * H44 - g4 * r1 * r2;
  out(j, 17) = r2 * r2 * H44 + g4 * r1 * r2;

  // Cross terms with param 4 (logitp): H(4,3)=col16, H(4,4)=col18
  double H43 = out(j, 16), H44p = out(j, 18);
  out(j, 16) = H43 + H44p * r1;
  out(j, 18) = H44p * r2;
}
  }
}

   return out;
}



// //' Discrete Extended generalized Pareto distribution of type 3 (deGPD3) negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for each deGPD parameter
// //' @param X1 a design matrix for the deGPD log scale parameter
// //' @param X2 a design matrix for the deGPD log shape parameter
// //' @param X3 a design matrix for the deGPD log delta
// //' @param yvec a vector
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return degpd3d0 a scalar, the negative log-liklihood
// //' @return degpd3d12 a matrix, first then second derivatives w.r.t. deGPD3 parameters
// //' @return degpd3d34 a matrix, third then fourth derivatives w.r.t. deGPD3 parameters (Not given)
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double degpd3d0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2, const arma::mat& X3, arma::vec yvec, const arma::uvec& dupid, int dcate, const Rcpp::List& offsets)
{
  
  arma::vec lsigmavec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lxivec = X2 * Rcpp::as<arma::vec>(pars[1]);
  arma::vec ldeltavec = X3 * Rcpp::as<arma::vec>(pars[2]);
  int nobs = yvec.size();
  
  if (dcate == 1) {
    lsigmavec = lsigmavec.elem(dupid);
    lxivec = lxivec.elem(dupid);
    ldeltavec = ldeltavec.elem(dupid);
  }

  {
  arma::vec off0 = Rcpp::as<arma::vec>(offsets[0]);
  if (off0.n_elem > 0) lsigmavec += off0;
  arma::vec off1 = Rcpp::as<arma::vec>(offsets[1]);
  if (off1.n_elem > 0) lxivec += off1;
  arma::vec off2 = Rcpp::as<arma::vec>(offsets[2]);
  if (off2.n_elem > 0) ldeltavec += off2;
  }

  double y, lsigma, lxi, ldelta;
  double e1,e2,e3, e4,e5,e6,e7, e8, e9,e10,e11, e12, e13;
  double lo, hi;
  double nllh=0.0;
  
  for (int j=0; j < nobs; j++) {
    
    y = yvec[j];
    lsigma = lsigmavec[j];
    lxi = lxivec[j];
    ldelta = ldeltavec[j];
    
    e1=1.0/exp(lxi);
    e2 =  exp(lxi) / exp(lsigma);
    e3= 1+ (y+1)*e2;
    e4= R_pow(1/e3, e1);
    e5= 1-e4; //H(y+1)
    e6= exp(ldelta)+1; 
    e7= R_pow(e4, e6)/ exp(ldelta); //[bar H(y+1)]^(ldelta+1)//ldelta
    e8=e4/ exp(ldelta);  //[bar H(y+1)]//ldelta
	hi= e5+e7-e8;
	e9= 1+ y*e2;
	e10= R_pow(1/e9, e1);
    e11= 1-e10; //H(y)
    e12= R_pow(e10, e6)/ exp(ldelta); //[bar H(y)]^(ldelta+1)//ldelta
	e13=e10/ exp(ldelta);  //[bar H(y+1)]//ldelta
	lo= e11+e12-e13;
    
    nllh += -log(hi-lo);
    //if (!ISNA(nllh)){
    //nllh = 1e20;
    //break;
//}
// Rprintf("hi %f ",hi); 
  //  Rprintf("lo %f \n",lo);
 }
  return(nllh);
}


// //' @rdname degpd3d0
// [[Rcpp::export]]
arma::mat degpd3d12(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::vec yvec, const arma::uvec dupid, int dcate, const Rcpp::List& offsets)
{
  
  arma::vec lsigmavec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lxivec = X2 * Rcpp::as<arma::vec>(pars[1]);
  arma::vec ldeltavec = X3 * Rcpp::as<arma::vec>(pars[2]);
  int nobs = yvec.size();
  arma::mat out = arma::mat(nobs, 9);
  
  if (dcate == 1) {
    lsigmavec = lsigmavec.elem(dupid);
    lxivec = lxivec.elem(dupid);
    ldeltavec = ldeltavec.elem(dupid);
  }

  {
  arma::vec off0 = Rcpp::as<arma::vec>(offsets[0]);
  if (off0.n_elem > 0) lsigmavec += off0;
  arma::vec off1 = Rcpp::as<arma::vec>(offsets[1]);
  if (off1.n_elem > 0) lxivec += off1;
  arma::vec off2 = Rcpp::as<arma::vec>(offsets[2]);
  if (off2.n_elem > 0) ldeltavec += off2;
  }

  double y, lsigma, lxi, ldelta;
  double ee1, ee2, ee3, ee4, ee6, ee8, ee9;
  double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
  double ee20, ee21,ee22, ee23, ee24, ee25, ee26, ee27, ee28, ee29;
  double ee32, ee33, ee34, ee35, ee36, ee37, ee38, ee39, ee40, ee41, ee42,ee43, ee44, ee45, ee46, ee47, ee48;
  double  ee50, ee53, ee55, ee56, ee58, ee59, ee60, ee61, ee62;
  double ee63, ee64, ee65, ee67,ee69, ee72,ee73, ee74, ee76, ee77, ee81, ee82, ee83;
  double ee87,ee89,ee90,ee91, ee92,ee93, ee96, ee97, ee100;
  double ee103,ee106,ee107,ee108,ee111,ee112,ee115, ee119,ee123;
  double ee125,ee126,ee127,ee128,ee129,ee130,ee131,ee132, ee133;
  double ee135,ee136,ee138,ee139,ee140,ee141, ee143,ee144,ee145;
  
  double eee1, eee2, eee3, eee5, eee6, eee7, eee8, eee9, eee10, eee11, eee12, eee13, eee14, eee16;
  double eee17, eee19, eee20, eee21, eee22, eee23, eee24, eee25, eee26,eee27, eee28, eee29, eee30, eee31;
  double eee32, eee33, eee36, eee38, eee41, eee43, eee45, eee46,eee47, eee48, eee49, eee53, eee54;
  double eee56, eee58, eee61, eee62, eee63, eee64, eee65, eee66,eee67, eee71, eee73, eee74;
  
  
  
  
  for (int j=0; j < nobs; j++) {
    
    y  = yvec[j];
      lsigma  = lsigmavec[j];
      lxi  = lxivec[j];
      ldelta  = ldeltavec[j];
   
    if(y>0){   
    ee1 = exp(lxi);
    ee2 = exp(lsigma);
    ee3 = 1 + y;
    ee4 = exp(ldelta);
    ee6 = ee3 * ee1/ee2;
    ee8 = y * ee1/ee2;
    ee9 = 1/ee1;
    ee10 = 1 + ee4;
    ee11 = ee6 + 1;
    ee12 = 1 + ee8;
    ee13 = ee10/ee1;
    ee14 = 1 + ee9;
    ee15 = R_pow(ee11,ee9);
    ee16 = R_pow(ee12,ee9);
    ee17 = log1p(ee6);
    ee18 = log1p(ee8);
    ee19 = R_pow(ee11,ee13);
    ee20 = R_pow(ee12,ee13);
    ee21 = ee13 + 1;
    ee22 = R_pow(ee11,ee14);
    ee23 = R_pow(ee12,ee14);
    ee24 = 1/ee15;
    ee25 = 1/ee16;
    ee26 = 1/ee19;
    ee27 = 1/ee20;
    ee28 = ee25 - ee24;
    ee29 = ee4/ee1;
    ee32 = ee10 * ee28 + ee26 - ee27;
    ee33 = ee9 - 1;
    ee34 = R_pow(ee11,ee21);
    ee35 = R_pow(ee12,ee21);
    ee36 = ee22 * ee2;
    ee37 = ee23 * ee2;
    ee38 = ee15 * ee1;
    ee39 = ee13 - 1;
    ee40 = ee16 * ee1;
    ee41 = 2 * ee13;
    ee42 = 2/ee1;
    ee43 = R_pow(ee11,ee29);
    ee44 = R_pow(ee12,ee29);
    ee45 = 1/ee22;
    ee46 = 1/ee23;
    ee47 = ee32 * ee2;
    ee48 = ee3 * (1/ee34 - ee45);
    ee50 = ee3/ee22;
    ee53 = 2 * ee14;
    ee55 = ee17/ee38 - ee3/ee36;
    ee56 = ee17/ee19;
    ee58 = ee18/ee40 - y/ee37;
    ee59 = ee18/ee20;
    ee60 = y * (ee46 - 1/ee35);
    ee61 = y/ee23;
    ee62 = R_pow(ee11,ee33);
    ee63 = R_pow(ee12,ee33);
    ee64 = ee48 + ee60;
    ee65 = (ee59 - ee56)/ee1;
    ee67 = ee62 * ee3/ee2;
    ee69 = ee15 * ee17/ee1;
    ee72 = ee16 * ee18/ee1;
    ee73 = 1 - ee10/ee4;
    ee74 = 2 * ee21;
    ee76 = y * ee63/ee2;
    ee77 = (ee67 - ee69)/R_pow(ee11,ee42);
    ee81 = R_pow(ee11,ee39) * ee3/ee2 - ee19 * ee17/ee1;
    ee82 = R_pow(ee11,ee41);
    ee83 = (ee50 - ee61)/ee2;
    ee87 = R_pow(ee12,ee41);
    ee89 = ee55/ee43;
    ee90 = ee58/ee44;
    ee91 = (ee76 - ee72)/R_pow(ee12,ee42);
    ee92 = ee17/ee15;
    ee93 = ee18/ee16;
    ee96 = y * R_pow(ee12,ee39)/ee2 - ee20 * ee18/ee1;
    ee97 = ee81/ee82;
    ee100 = ee15 * ee14 * ee3 * ee1;
    ee103 = ee83 + ee89 + (ee93 - ee92)/ee1 - ee90;
    ee106 = ee73 * ee28 + (ee27 - ee26)/ee4 + ee65;
    ee107 = ee96/ee87;
    ee108 = ee9 - ee53;
    ee111 = y * ee14 * ee16 * ee1;
    ee112 = R_pow(ee47,2);
    ee115 = ee77 + ee107 - (ee97 + ee91);
    ee119 = ee100/ee2 - ee22 * ee17/ee1;
    ee123 = ee65 + ee25 - ee24;
    ee125 = ee111/ee2 - ee23 * ee18/ee1;
    ee126 = R_pow(ee36,2);
    ee127 = R_pow(ee38,2);
    ee128 = R_pow(ee37,2);
    ee129 = R_pow(ee40,2);
    ee130 = ee64 * ee10;
    ee131 = R_pow(ee11, ee108);
    ee132 = R_pow(ee11,ee53);
    ee133 = ee13 - ee74;
    //ee134 = ee21 - ee74;
    ee135 = R_pow(ee12,ee108);
    ee136 = R_pow(ee12,ee53);
    ee138 = R_pow(ee3,2);
    ee139 = 1 + ee29;
    ee140 = 2 * ee29;
    ee141 = ee29 - 1;
    ee143 = y/ee35 - ee3/ee34;
    ee144 = ee61 - ee50;
    ee145 = R_pow(y,2);

    out(j, 0) = -(ee130/ee47);
    out(j, 1) = -(ee103 * ee10/ee32); 
    out(j, 2) = -(ee106 * ee4/ee32);
    out(j, 3) = -((((ee21 * 
                       R_pow(ee11,ee133) - ee131 * ee14) * ee138 + ee145 * (ee14 * 
                                                                       ee135 - ee21 * R_pow(ee12,ee133))) * ee1/(ee32 * R_pow(ee2,2)) - (ee47 + 
                                                                                                                             ee130) * ee64/ee112) * ee10);   // weret (lsigma, lsigma);
    out(j, 4) = -(((((ee36 - ee100)/ee126 + 
                       (ee62 * ee1 * ee17/ee127 - ee45)/ee2)/ee43 + ee4 * ee55/(R_pow(ee11,ee139) * 
                                                                                  ee2)) * ee3 + ((ee131 * ee138 - ee145 * ee135) * ee14 * 
                                                                                                   ee1/ee2 + (y * (ee18/ee23 - ee1/ee23) - ee3 * (ee17/ee22 - 
                                                                                                                                                    ee1/ee22))/ee1 + ee61 - (ee103 * ee64 * ee10/ee32 + ee50))/ee2 - 
                     y * (((ee37 - ee111)/ee128 + (ee63 * ee1 * ee18/ee129 - 
                                                     ee46)/ee2)/ee44 + ee4 * ee58/(R_pow(ee12, ee139) * ee2))) * 
                    ee10/ee32); // weret (lsigma, xi)
    out(j, 5) = -(((ee143/ee4 - ee106 * ee64/ee32) * 
                     ee10 + ee73 * ee144 + (y * (ee10 * ee18/ee35 - ee1/ee35) - 
                                              (ee10 * ee17/ee34 - ee1/ee34) * ee3)/ee1) * ee4/ee47); // weret (lsigma, ldelta) 
    out(j, 6) = -((((ee77 + ee24) * ee17 + ee1 * ee144/ee2 - 
                      (ee91 + ee25) * ee18)/ee1 + ((ee119 * ee2/ee126 + 
                                                      1/ee36) * ee3 - (ee67 + ee15 - ee69) * ee1 * ee17/ee127)/ee43 + 
                     (ee58 * (y * R_pow(ee12, ee141)/ee2 - ee44 * ee18/ee1)/R_pow(ee12, ee140) - 
                        (R_pow(ee11,ee141) * ee3/ee2 - ee43 * ee17/ee1) * ee55/R_pow(ee11,ee140)) * 
                     ee4 + (y * ee125/ee136 - ee119 * ee3/ee132)/ee2 - 
                     (ee115 * ee103 * ee10/ee32 + (y * (1/ee37 + ee2 * 
                                                          ee125/ee128) - (ee16 + ee76 - ee72) * ee1 * ee18/ee129)/ee44)) * ee10/ee32); // weret (lxi, lxi)
    out(j, 7) =  -((((ee81 * ee17/ee82 - ee18 * 
                        ee96/ee87) * ee10 + ee1 * ee143/ee2 + ee56 - ee59)/ee1 + 
                      ((ee97 - ee107)/ee4 - ee115 * ee106/ee32) * ee10 + 
                      (ee77 - ee91) * ee73) * ee4/ee32); // weret (lxi, ldelta)
    out(j, 8) = -(((R_pow(ee17,2)/ee19 - 
                      R_pow(ee18,2)/ee20)/R_pow(ee1,2) - ee106 * ee123/ee32) * R_pow(ee4, 2)/ee32); // w.r.t (ldelta, ldelta)
  }

else {
     	   
    eee1 = exp(lxi);
    eee2 = exp(lsigma);
    eee3 = 1 + y;
    eee5 = eee3 * eee1/eee2;
    eee6 = exp(ldelta);
    eee7 = eee5 + 1;
    eee8 = 1 + eee6;
    eee9 = 1/eee1;
    eee10 = eee8/eee1;
    eee11 = R_pow(eee7,eee9);
    eee12 = R_pow(eee7,eee10);
    eee13 = log1p(eee5);
    eee14 = 1 + eee9;
    eee16 = 1/eee12 + eee6;
    eee17 = R_pow(eee7,eee14);
    eee19 = eee16 - eee8/eee11;
    eee20 = eee10 + 1;
    eee21 = R_pow(eee7,eee20);
    eee22 = eee12 * eee1;
    eee23 = eee17 * eee2;
    eee24 = eee6/eee1;
    eee25 = 1/eee21;
    eee26 = 1/eee17;
    eee27 = eee11 * eee1;
    eee28 = eee19 * eee2;
    eee29 = R_pow(eee7,eee24);
    eee30 = eee25 - eee26;
    eee31 = eee13/eee22;
    eee32 = R_pow(eee7,(eee10 - 1));
    eee33 = R_pow(eee7,(eee9 - 1));
    eee36 = 1 - eee8/eee6;
    eee38 = eee13/eee27 - eee3/eee23;
    eee41 = eee32 * eee3/eee2 - eee12 * eee13/eee1;
    eee43 = eee33 * eee3/eee2;
    eee45 = eee11 * eee13/eee1;
    eee46 = 1/eee29;
    eee47 = eee43 - eee45;
    eee48 = R_pow(eee7,(2 * eee10));
    eee49 = R_pow(eee7,(2/eee1));
    eee53 = 1 - (eee36/eee11 + eee16/eee6 + eee31);
    eee54 = eee46 - 1;
    eee56 = R_pow(eee22,2);
    eee58 = eee47/eee49 - eee41/eee48;
    eee61 = eee11 * eee14 * eee3 * eee1;
    eee62 = R_pow(eee28,2);
    eee63 = 1 - (1/eee11 + eee31);
    eee64 = 2 * eee20;
    eee65 = eee58 * eee8;
    eee66 = R_pow(eee23,2);
    eee67 = R_pow(eee27,2);
    eee71 = eee61/eee2 - eee17 * eee13/eee1;
    eee73 = eee8 * eee3 * eee30;
    eee74 = 2 * eee14;
    
	
	out(j, 0) = -(eee73/eee28);  // weereetee lsigma
    out(j, 1) = -(eee8 * eee54 * 
                    eee38/eee19); //# weereetee lxi
    out(j, 2) = -(eee53 * eee6/eee19); //# weereetee ldelta 
    out(j, 3) = -(((eee20 * 
                      R_pow(eee7,(eee10 - eee64)) - R_pow(eee7,(eee9 - eee74)) * eee14) * eee3 * 
                     eee1/(eee19 * R_pow(eee2,2)) - (eee73 + eee28) * eee30/eee62) * eee8 * 
                    eee3);  //# weereet (lsigma, lsigma)
    out(j, 4) = -((((eee23 - eee61)/eee66 + (eee33 * eee1 * eee13/eee67 - 
                                            (eee8 * eee30 * eee38/eee19 + eee26))/eee2) * eee54 + eee6 * 
                     eee38/(R_pow(eee7,(1 + eee24)) * eee2)) * eee8 * eee3/eee19); //# weereet (lsigma, xi)
    out(j, 5) = ((eee32 * 
                    eee1 * eee13/eee56 + eee53 * eee30/eee19 + 1/(eee21 * eee6)) * 
                   eee8 + eee36/eee17 - eee25) * eee3 * eee6/eee28;   // # weereet (lsigma, ldelta) 
    out(j, 6) = -((((eee71 * 
                       eee2/eee66 + 1/eee23) * eee3 - (eee65 * eee38/eee19 + (eee43 + 
                                                                         eee11 - eee45) * eee1 * eee13/eee67)) * eee54 - (R_pow(eee7,(eee24 - 
                                                                                                                           1)) * eee3/eee2 - eee29 * eee13/eee1) * eee6 * eee38/R_pow(eee7,(2 *eee24))) * eee8/eee19); //# weereet (lxi, lxi)
    out(j, 7) =  ((eee58 * eee53/eee19 - eee41/(eee48 * 
                                             eee6)) * eee8 + eee3/(eee21 * eee2) - ((eee41 * eee8 + eee12) * 
                                                                                 eee1 * eee13/eee56 + eee47 * eee36/eee49)) * eee6/eee19; //# weereet (lxi, ldelta) 
    out(j, 8) =  -((eee12 * R_pow(eee13,2)/eee56 - eee53 * 
                     eee63/eee19) * R_pow(eee6,2)/eee19);   //# weereet (ldelta, ldelta)                                                                           
  }
}
     
   return out;
}
    

// //' Discrete Extended generalized Pareto distribution of type 4 (deGPD4) negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for each deGPD parameter
// //' @param X1 a design matrix for the deGPD log scale parameter
// //' @param X2 a design matrix for the deGPD log shape parameter
// //' @param X3 a design matrix for the deGPD log delta
// //' @param X3 a design matrix for the deGPD log kappa
// //' @param yvec a vector
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return degpd4d0 a scalar, the negative log-liklihood
// //' @return degpd4d12 a matrix, first then second derivatives w.r.t. deGPD4 parameters
// //' @return degpd4d34 a matrix, third then fourth derivatives w.r.t. deGPD4 parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]

double degpd4d0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2, const arma::mat& X3, const arma::mat& X4, arma::vec yvec, const arma::uvec& dupid, int dcate, const Rcpp::List& offsets)
{
  arma::vec lsigmavec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lxivec = X2 * Rcpp::as<arma::vec>(pars[1]);
  arma::vec ldeltavec = X3 * Rcpp::as<arma::vec>(pars[2]);
  arma::vec lkappavec = X4 * Rcpp::as<arma::vec>(pars[3]);
  int nobs = yvec.size();
  
 if (dcate == 1) {
  lsigmavec = lsigmavec.elem(dupid);
  lxivec = lxivec.elem(dupid);
  ldeltavec = ldeltavec.elem(dupid);
  lkappavec = lkappavec.elem(dupid);
}

  {
  arma::vec off0 = Rcpp::as<arma::vec>(offsets[0]);
  if (off0.n_elem > 0) lsigmavec += off0;
  arma::vec off1 = Rcpp::as<arma::vec>(offsets[1]);
  if (off1.n_elem > 0) lxivec += off1;
  arma::vec off2 = Rcpp::as<arma::vec>(offsets[2]);
  if (off2.n_elem > 0) ldeltavec += off2;
  arma::vec off3 = Rcpp::as<arma::vec>(offsets[3]);
  if (off3.n_elem > 0) lkappavec += off3;
  }

  double y, lsigma, lxi, ldelta, lkappa;
  double e1,e2,e3, e4,e5,e6,e7, e8, e9,e10,e11, e12, e13;
  double lo, hi;
  double nllh=0.0;

  for (int j=0; j < nobs; j++) {

    y = yvec[j];
    lsigma = lsigmavec[j];
    lxi = lxivec[j];
    ldelta = ldeltavec[j];
    lkappa = lkappavec[j];
    
    e1=1.0/exp(lxi);
    e2 =  exp(lxi) / exp(lsigma);
    e3= 1+ (y+1)*e2;
    e4= R_pow(1/e3, e1);
    e5= 1-e4; //H(y+1)
    e6= exp(ldelta)+1; 
    e7= R_pow(e4, e6)/ exp(ldelta); //[bar H(y+1)]^(ldelta+1)//ldelta
    e8=e4/ exp(ldelta);  //[bar H(y+1)]//ldelta
	hi= R_pow((e5+e7-e8),exp(lkappa)/2);
	e9= 1+ y*e2;
	e10= R_pow(1/e9, e1);
    e11= 1-e10; //H(y)
    e12= R_pow(e10, e6)/ exp(ldelta); //[bar H(y)]^(ldelta+1)//ldelta
	e13=e10/ exp(ldelta);  //[bar H(y+1)]//ldelta
	lo= R_pow((e11+e12-e13),exp(lkappa)/2);
    
    nllh += -log(hi-lo);
    //if (!ISNA(nllh)){
    //nllh = 1e20;
    // break;
//}
// Rprintf("hi %f ",hi); 
  //  Rprintf("lo %f \n",lo);
 }
  return(nllh);
}


// //' @rdname degpd4d0
// [[Rcpp::export]]
arma::mat degpd4d12(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::mat X3,  arma::mat X4, arma::vec yvec, const arma::uvec dupid, int dcate, const Rcpp::List& offsets)
{
  
  arma::vec lsigmavec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lxivec = X2 * Rcpp::as<arma::vec>(pars[1]);
  arma::vec ldeltavec = X3 * Rcpp::as<arma::vec>(pars[2]);
  arma::vec lkappavec = X4 * Rcpp::as<arma::vec>(pars[3]);
  int nobs = yvec.size();
  arma::mat out = arma::mat(nobs, 14);
  
 if (dcate == 1) {
  lsigmavec = lsigmavec.elem(dupid);
  lxivec = lxivec.elem(dupid);
  ldeltavec = ldeltavec.elem(dupid);
  lkappavec = lkappavec.elem(dupid);
}

  {
  arma::vec off0 = Rcpp::as<arma::vec>(offsets[0]);
  if (off0.n_elem > 0) lsigmavec += off0;
  arma::vec off1 = Rcpp::as<arma::vec>(offsets[1]);
  if (off1.n_elem > 0) lxivec += off1;
  arma::vec off2 = Rcpp::as<arma::vec>(offsets[2]);
  if (off2.n_elem > 0) ldeltavec += off2;
  arma::vec off3 = Rcpp::as<arma::vec>(offsets[3]);
  if (off3.n_elem > 0) lkappavec += off3;
  }

  double y, lsigma, lxi, ldelta, lkappa;
  double ee1, ee2, ee3, ee4, ee5, ee7, ee9, ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
  double ee20, ee21,ee22, ee23, ee24, ee25, ee26, ee28, ee30, ee31,  ee32, ee33, ee34, ee35, ee36, ee37, ee38, ee39, ee40, ee41, ee42,ee43, ee44, ee45, ee46, ee47, ee48;
  double ee49, ee50, ee51,ee52,ee53,ee54, ee55, ee56, ee57,ee58, ee59, ee60, ee61, ee62;
  double ee64, ee66, ee67,ee68,ee70, ee72,ee73, ee74,ee75, ee76, ee77, ee78,ee79,ee80,ee81, ee83,ee84,ee85,ee86;
  double ee87,ee88,ee89,ee90,ee91, ee92,ee93,ee95, ee96, ee97,ee98,ee99, ee100;
  double ee101,ee102, ee103, ee105,ee107,ee108, ee109, ee110,ee111,ee112, ee113,ee115, ee117, ee118,ee119,ee121,ee123, ee124;
  double ee125,ee126,ee127,ee128,ee129,ee130,ee131,ee132, ee133, ee134;
  double ee135,ee138,ee139,ee140,ee141,ee144,ee145,ee146 ,ee149,ee150,ee151,ee152,ee153,ee154,ee158,ee163,ee166;
  double ee168,ee169,ee172,ee173,ee177,ee181,ee182,ee183,ee186,ee188,ee189,ee192,ee193,ee194,ee195,ee197,ee200;
  double ee204,ee206,ee208,ee211,ee215,ee217,ee218,ee222,ee226,ee230,ee231,ee235,ee236;
  
  double eee1, eee2, eee3, eee5, eee6, eee7, eee8, eee9, eee10, eee11, eee12, eee13, eee14,eee15, eee16;
  double eee17, eee18,eee19, eee20,  eee22, eee23, eee24, eee25, eee26,eee27, eee28, eee29, eee30, eee31;
  double eee33, eee35, eee36, eee37,eee38,eee39,eee40, eee41, eee42,eee43, eee44,eee45, eee46,eee47, eee48, eee49,eee52, eee53, eee54;
  double eee55, eee57,eee58, eee60,eee61, eee62, eee63,eee67,eee68, eee69,eee70,eee71,eee74,eee75,  eee76;
  double eee77,eee84,eee85,eee87,eee90,eee93, eee97,eee99, eee100;
  
  
  
  for (int j=0; j < nobs; j++) {
    
    y  = yvec[j];
      lsigma  = lsigmavec[j];
      lxi  = lxivec[j];
      ldelta  = ldeltavec[j];
      lkappa  = lkappavec[j];
   
    if(y>0){   
      ee1 = exp(lxi);
      ee2 = exp(lsigma);
      ee3 = exp(ldelta);
      ee4 = 1 + y;
      ee5 = 1/ee1;
      ee7 = ee4 * ee1/ee2;
      ee9 = y * ee1/ee2;
      ee10 = ee7 + 1;
      ee11 = 1 + ee9;
      ee12 = 1 + ee3;
      ee13 = ee12/ee1;
      ee14 = R_pow(ee10,ee5);
      ee15 = R_pow(ee11,ee5);
      ee16 = exp(lkappa);
      ee17 = 1/ee14;
      ee18 = 1/ee15;
      ee19 = ee16/2;
      ee20 = R_pow(ee10,ee13);
      ee21 = R_pow(ee11,ee13);
      ee22 = 1 + ee5;
      ee23 = 1/ee20;
      ee24 = 1/ee21;
      ee25 = (ee23 -ee17)/ee3;
      ee26 = (ee24 - ee18)/ee3;
      ee28 = ee25 + 1 - ee17;
      ee30 = ee26 + 1 - ee18;
      ee31 = ee19 - 1;
      ee32 = ee13 + 1;
      ee33 = log1p(ee7);
      ee34 = log1p(ee9);
      ee35 = R_pow(ee10,ee22);
      ee36 = R_pow(ee11,ee22);
      ee37 = R_pow(ee28,ee31);
      ee38 = R_pow(ee30,ee31);
      ee39 = R_pow(ee28,ee19);
      ee40 = R_pow(ee30,ee19);
      ee41 = R_pow(ee10,ee32);
      ee42 = R_pow(ee11,ee32);
      ee43 = ee5 - 1;
      ee44 = 1/ee35;
      ee45 = 1/ee36;
      ee46 = ee20 * ee1;
      ee47 = ee39 - ee40;
      ee48 = ee21 * ee1;
      ee49 = 2/ee1;
      ee50 = ee35 * ee2;
      ee51 = ee36 * ee2;
      ee52 = ee33/ee46;
      ee53 = ee34/ee48;
      ee54 = ee13 - 1;
      ee55 = ee14 * ee1;
      ee56 = ee15 * ee1;
      ee57 = R_pow(ee10,ee43);
      ee58 = R_pow(ee11,ee43);
      ee59 = ee4/ee50;
      ee60 = ee33/ee55;
      ee61 = ee34/ee56;
      ee62 = y/ee51;
      ee64 = ee57 * ee4/ee2;
      ee66 = ee14 * ee33/ee1;
      ee67 = ee12/ee41;
      ee68 = ee12/ee42;
      ee70 = ee15 * ee34/ee1;
      ee72 = y * ee58/ee2;
      ee73 = (ee67 - ee44)/ee3;
      ee74 = (ee68 - ee45)/ee3;
      ee75 = ee41 * ee2;
      ee76 = ee42 * ee2;
      ee77 = 2 * ee13;
      ee78 = ee19 - 2;
      ee79 = (ee64 - ee66)/R_pow(ee10,ee49);
      ee80 = ee73 - ee44;
      ee81 = ee74 - ee45;
      ee83 = (ee72 - ee70)/R_pow(ee11,ee49);
      ee84 = log(ee28);
      ee85 = log(ee30);
      ee86 = R_pow(ee10,ee54);
      ee87 = R_pow(ee11,ee54);
      ee88 = ee4/ee75;
      ee89 = y/ee76;
      ee90 = ee47 * ee2;
      ee91 = ee20 * ee33;
      ee92 = ee21 * ee34;
      ee93 = y * ee81;
      ee95 = ee80 * ee37 * ee4;
      ee96 = (ee86 * ee4/ee2 - ee91/ee1) * ee12;
      ee97 = (ee12 * (ee52 - ee88) + ee59 - ee60)/ee3;
      ee98 = (ee12 * (ee53 - ee89) + ee62 - ee61)/ee3;
      ee99 = ee12 * (y * ee87/ee2 - ee92/ee1);
      ee100 = (ee17 - ee23)/ee3;
      ee101 = (ee18 - ee24)/ee3;
      ee102 = ee93 * ee38;
      ee103 = ee95/2;
      ee105 = ee97 + ee59 - ee60;
      ee107 = ee98 + ee62 - ee61;
      ee108 = ee100 - ee52;
      ee109 = ee101 - ee53;
      ee110 = ee102/2;
      ee111 = ee103 - ee110;
      ee112 = ee96/R_pow(ee10,ee77);
      ee113 = R_pow(ee28,ee78);
      ee115 = R_pow(ee30,ee78);
      ee117 = ee99/R_pow(ee11,ee77);
      ee118 = 0.5 * (ee39 * ee84);
      ee119 = 0.5 * (ee40 * ee85);
      ee121 = (ee79 - ee112)/ee3 + ee79;
      ee123 = (ee83 - ee117)/ee3 + ee83;
      ee124 = ee25 + ee52;
      ee125 = ee26 + ee53;
      ee126 = ee118 - ee119;
      ee127 = 2 * ee22;
      ee128 = ee105 * ee37;
      ee129 = ee107 * ee38;
      ee130 = ee37 * ee108;
      ee131 = ee38 * ee109;
      ee132 = ee128/2;
      ee133 = ee129/2;
      ee134 = R_pow(ee46,2);
      ee135 = R_pow(ee48,2);
      ee138 = ee14 * ee22 * ee4 * ee1;
      ee139 = ee130/2;
      ee140 = ee131/2;
      ee141 = 2 * ee32;
      ee144 = y * ee22 * ee15 * ee1;
      ee145 = ee132 - ee133;
      ee146 = ee139 - ee140;
      ee149 = ee121 * ee37/2 - ee123 * ee38/2;
      ee150 = R_pow(ee90,2);
      ee151 = R_pow(ee50,2);
      ee152 = R_pow(ee55,2);
      ee153 = R_pow(ee51,2);
      ee154 = R_pow(ee56,2);
      ee158 = ee138/ee2 - ee35 * ee33/ee1;
      ee163 = ee38 * ee125/2 - ee37 * ee124/2;
      ee166 = ee5 - ee127;
      ee168 = ee144/ee2 - ee36 * ee34/ee1;
      ee169 = R_pow(ee75,2);
      ee172 = ee32 * ee20 * ee4 * ee1;
      ee173 = R_pow(ee76,2);
      ee177 = ee37 + ee37 * ee16 * ee84/2;
      ee181 = ee38 + ee38 * ee16 * ee85/2;
      ee182 = 1/ee41;
      ee183 = 1/ee42;
      ee186 = y * ee32 * ee21 * ee1;
      ee188 = ee111 * ee126 * ee16;
      ee189 = ee111 * ee16;
      ee192 = (ee96 + ee20) * ee1 * ee33/ee134;
      ee193 = (ee158 * ee2/ee151 + 1/ee50) * ee4;
      ee194 = ee80 * ee113;
      ee195 = ee81 * ee115;
      ee197 = (ee50 - ee138)/ee151 + (ee57 * ee1 * ee33/ee152 -ee44)/ee2;
      ee200 = (ee64 + ee14 - ee66) * ee1 * ee33/ee152;
      //ee201 = ee158/R_pow(ee10,ee127);
      ee204 = (ee99 + ee21) * ee1 * ee34/ee135;
      ee206 = ee172/ee2 - ee41 * ee12 * ee33/ee1;
      ee208 = (ee51 - ee144)/ee153 + (ee58 * ee1 * ee34/ee154 - ee45)/ee2;
      ee211 = (ee15 + ee72 - ee70) * ee1 * ee34/ee154;
      ee215 = ee86 * ee12 * ee1 * ee33/ee134;
      ee217 = R_pow(ee10,ee166) * ee22;
      ee218 = ee22 * R_pow(ee11,ee166);
      ee222 = ee12 * ee87 * ee1 * ee34/ee135;
      ee226 = ee13 - ee141;
      //ee227 = ee32 - ee141;
     // ee229 = ee168/R_pow(ee11,ee127);
      ee230 = 1/ee46;
      ee231 = 1/ee48;
      ee235 = ee186/ee2 - ee12 * ee42 * ee34/ee1;
      ee236 = y * (1/ee51 + ee2 * ee168/ee153);

    out(j, 0) = -(ee189/ee90); // # werete lsigma
    out(j, 1) = -(ee145 * ee16/ee47); // # werete lxi 
    out(j, 2) = -(ee146 * ee16/ee47); //# werete ldelta
    out(j, 3) = -(ee126 * ee16/ee47); //# werete lkappa
    
	
	out(j, 4) = -(((0.5 * ((((ee32 * R_pow(ee10,ee226) * 
                                ee12 - ee217)/ee3 - ee217) * ee37 * ee1 + R_pow(ee80,2) * 
                              ee113 * ee31) * R_pow(ee4,2)) - 0.5 * (R_pow(y,2) * (((ee32 * ee12 * 
                                                                         R_pow(ee11,ee226) - ee218)/ee3 - ee218) * ee38 * ee1 + R_pow(ee81,2) * 
                                                                       ee115 * ee31)))/(ee47 * R_pow(ee2,2)) - (ee189 + ee90) * 
                     ee111/ee150) * ee16);   //# weret (lsigma, lsigma)
    out(j, 5) = -((0.5 * ((((((ee215 - 
                                 ee182)/ee2 + (ee75 - ee172)/ee169) * ee12 - ee197)/ee3 - 
                              ee197) * ee37 + ee105 * ee80 * ee113 * ee31/ee2) * 
                            ee4) - (ee145 * ee111 * ee16/ee90 + 0.5 * (y * (((((ee222 - 
                                                                                  ee183)/ee2 + (ee76 - ee186)/ee173) * ee12 - ee208)/ee3 - 
                                                                               ee208) * ee38 + ee107 * ee81 * ee115 * ee31/ee2)))) * 
                    ee16/ee47); //# weret (lsigma, xi)
    out(j, 6) = -((0.5 * ((ee194 * ee108 * ee31 + 
                             ee37 * ((ee44 - ee67)/ee3 + ee182 - ee215)) * ee4) - 
                     (ee111 * ee146 * ee16/ee47 + 0.5 * (y * (ee195 * 
                                                                ee109 * ee31 + ee38 * ((ee45 - ee68)/ee3 + ee183 - 
                                                                                         ee222))))) * ee16/ee90); //# weret (lsigma, ldelta)
    out(j, 7) = -((0.5 * (ee80 * 
                             ee177 * ee4) - (ee188/ee47 + 0.5 * (ee93 * ee181))) * 
                     ee16/ee90); //# weret (lsigma, lkappa)
    out(j, 8) = -((0.5 * (((((ee206 * ee2/ee169 + 1/ee75) * 
                               ee4 - ee192) * ee12 + ee200 - ee193)/ee3 + ee200 - 
                             ee193) * ee37 + ee121 * ee105 * ee113 * ee31) - (ee149 * 
                                                                                ee145 * ee16/ee47 + 0.5 * (((ee211 + ee12 * (y * 
                                                                                                                               (1/ee76 + ee2 * ee235/ee173) - ee204) - ee236)/ee3 + 
                                                                                                              ee211 - ee236) * ee38 + ee107 * ee123 * ee115 * ee31))) * 
                    ee16/ee47); // # weret (lxi, lxi)
out(j, 9) =  -((0.5 * ((ee192 + (ee112 - 
                                        ee79)/ee3 - ee88) * ee37 + ee121 * ee113 * ee108 * 
                              ee31) - (ee149 * ee146 * ee16/ee47 + 0.5 * ((ee204 + 
                                                                             (ee117 - ee83)/ee3 - ee89) * ee38 + ee123 * ee115 * 
                                                                            ee109 * ee31))) * ee16/ee47); //# weret (lxi, ldelta)
out(j, 10) = -((0.5 * (ee121 * 
                             ee177) - (ee149 * ee126 * ee16/ee47 + 0.5 * (ee181 * 
                                                                            ee123))) * ee16/ee47); //# weret (lxi, lkappa)
out(j, 11) = -((0.5 * (((ee20 * 
                               ee3 * ee33/ee134 + ee230) * ee33 - ee100) * ee37 - 
                             ee113 * ee124 * ee108 * ee31) - (ee146 * ee163 * 
                                                                ee16/ee47 + 0.5 * (((ee21 * ee3 * ee34/ee135 + ee231) * 
                                                                                      ee34 - ee101) * ee38 - ee115 * ee125 * ee109 * ee31))) * 
                     ee16/ee47); //# weret (ldelta, ldelta) 
out(j, 12) = -((0.5 * (ee181 * ee125) - (ee163 * 
                                               ee126 * ee16/ee47 + 0.5 * (ee177 * ee124))) * ee16/ee47); //# weret (ldelta, lkappa)
out(j, 13) = -(((0.25 * 
                        (ee39 * R_pow(ee84,2)) - (R_pow(ee126,2)/ee47 + 0.25 * (ee40 * 
                                                                    R_pow(ee85,2)))) * ee16 + ee118 - ee119) * ee16/ee47); //# weret (lkappa, lkappa)
 

  }

else {
     	   
    eee1 = exp(lxi);
    eee2 = exp(lsigma);
    eee3 = 1 + y;
    eee5 = eee3 * eee1/eee2;
    eee6 = eee5 + 1;
    eee7 = exp(ldelta);
    eee8 = 1/eee1;
    eee9 = 1 + eee7;
    eee10 = eee9/eee1;
    eee11 = R_pow(eee6,eee8);
    eee12 = R_pow(eee6,eee10);
    eee13 = 1/eee11;
    eee14 = log1p(eee5);
    eee15 = 1 + eee8;
    eee16 = R_pow(eee6,eee15);
    eee17 = 1/eee12;
    eee18 = eee10 + 1;
    eee19 = (eee17 - eee13)/eee7;
    eee20 = R_pow(eee6,eee18);
    eee22 = eee19 + 1 - eee13;
    eee23 = eee12 * eee1;
    eee24 = 1/eee16;
    eee25 = exp(lkappa);
    eee26 = eee16 * eee2;
    eee27 = 2 * eee22;
    eee28 = eee11 * eee1;
    eee29 = eee14/eee23;
    eee30 = eee20 * eee2;
    eee31 = R_pow(eee6,(eee8 - 1));
    eee33 = eee31 * eee3/eee2;
    eee35 = eee11 * eee14/eee1;
    eee36 = eee9/eee20;
    eee37 = eee3/eee26;
    eee38 = eee14/eee28;
    eee39 = (eee36 - eee24)/eee7;
    eee40 = (eee33 - eee35)/R_pow(eee6,(2/eee1));
    eee41 = eee39 - eee24;
    eee42 = R_pow(eee6,(eee10 - 1));
    eee43 = eee22 * eee2;
    eee44 = eee12 * eee14;
    eee45 = eee3/eee30;
    eee46 = (eee42 * eee3/eee2 - eee44/eee1) * eee9;
    eee47 = 2 * eee43;
    eee48 = R_pow(eee23,2);
    eee49 = (eee9 * (eee29 - eee45) + eee37 - eee38)/eee7;
    eee52 = eee11 * eee15 * eee3 * eee1;
    eee53 = (eee13 - eee17)/eee7;
    eee54 = R_pow(eee27,2);
    eee55 = eee46/R_pow(eee6,(2 * eee10));
    eee57 = eee49 + eee37 - eee38;
    eee58 = eee53 - eee29;
    eee60 = (eee40 - eee55)/eee7 + eee40;
    eee61 = eee41 * eee3;
    eee62 = R_pow(eee26,2);
    eee63 = R_pow(eee28,2);
    eee67 = eee52/eee2 - eee16 * eee14/eee1;
    eee68 = eee19 + eee29;
    eee69 = 2 * eee15;
    eee70 = eee61 * eee25;
    eee71 = R_pow(eee30,2);
    eee74 = eee18 * eee12 * eee3 * eee1;
    eee75 = R_pow(eee47,2);
    eee76 = 1/eee20;
    eee77 = 2 * eee18;
    eee84 = (eee46 + eee12) * eee1 * eee14/eee48;
    eee85 = (eee67 * eee2/eee62 + 1/eee26) * eee3;
    eee87 = (eee26 - eee52)/eee62 + (eee31 * eee1 * eee14/eee63 - eee24)/eee2;
    eee90 = (eee33 + eee11 - eee35) * eee1 * eee14/eee63;
   // eee91 = eee67/R_pow(eee6,eee69);
    eee93 = eee74/eee2 - eee20 * eee9 * eee14/eee1;
    eee97 = eee42 * eee9 * eee1 * eee14/eee48;
    eee99 = R_pow(eee6,(eee8 - eee69)) * eee15;
    eee100 = 1/eee23;

    
	
	out(j, 0) = -(eee70/eee47); //# weereetee lsigma
    out(j, 1) = -(eee57 * eee25/eee27); //# weereetee lxi
    out(j, 2) = -(eee58 * eee25/eee27); //# weereetee ldelta
    out(j, 3) = -(0.5 * (eee25 * 
                            log(eee22))); //# weereetee lkappa
    out(j, 4) = -((((eee18 * 
                       R_pow(eee6,(eee10 - eee77)) * eee9 - eee99)/eee7 - eee99) * eee3 * eee1/(2 * 
                                                                                  (eee22 * R_pow(eee2,2))) - 2 * ((eee61 + eee43) * eee41/eee75)) * 
                    eee3 * eee25); //# weereet (lsigma, lsigma)
    out(j, 5) = -((((((eee97 - eee76)/eee2 + (eee30 - 
                                             eee74)/eee71) * eee9 - eee87)/eee7 - eee87)/eee27 - 2 * (eee57 * 
                                                                                                 eee41/(eee54 * eee2))) * eee3 * eee25); //# weereet (lsigma, xi)
    out(j, 6) = -((((eee24 - 
                       eee36)/eee7 + eee76 - eee97)/eee27 - 2 * (eee41 * eee58/eee54)) * 
                    eee3 * eee25/eee2); // # weereet (lsigma, ldelta)
    out(j, 7) =    -(0.5 * (eee70/eee43)); //# weereet (lsigma, lkappa)
    out(j, 8) = -((((((eee93 * eee2/eee71 + 
                         1/eee30) * eee3 - eee84) * eee9 + eee90 - eee85)/eee7 + eee90 - 
                      eee85)/eee27 - 2 * (eee60 * eee57/eee54)) * eee25); // # weereet (lxi, lxi)  
	out(j, 9) = -(((eee84 + 
                       (eee55 - eee40)/eee7 - eee45)/eee27 - 2 * (eee60 * eee58/eee54)) * 
                     eee25); //# weereet (lxi, ldelta)
	out(j, 10) = -(0.5 * (eee60 * eee25/eee22)); // # weereet (lxi, lkappa)
	out(j, 11) = -((((eee12 * eee7 * eee14/eee48 + eee100) * 
                       eee14 - eee53)/eee27 + 2 * (eee68 * eee58/eee54)) * eee25); // # weereet (ldelta, ldelta)
	out(j, 12) = 0.5 * (eee68 * eee25/eee22); //  # weereet (ldelta, lkappa)
	out(j, 13) = -(0.5 * (eee25 *
                               log(eee22))); //# weer.t (lkappa, lkappa)

  }
}

   return out;
}


// //' Discrete Extended generalized Pareto distribution of type 5 (deGPD5) negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for each deGPD parameter
// //' @param X1 a design matrix for the deGPD log scale parameter
// //' @param X2 a design matrix for the deGPD log shape parameter
// //' @param X3 a design matrix for the deGPD log kappa
// //' @param yvec a vector
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return degpd5d0 a scalar, the negative log-liklihood
// //' @return degpd5d12 a matrix, first then second derivatives w.r.t. deGPD5 parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double degpd5d0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2, const arma::mat& X3, arma::vec yvec, const arma::uvec& dupid, int dcate, const Rcpp::List& offsets)
{

  arma::vec lsigmavec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lxivec = X2 * Rcpp::as<arma::vec>(pars[1]);
  arma::vec lkappavec = X3 * Rcpp::as<arma::vec>(pars[2]);
  int nobs = yvec.size();

  if (dcate == 1) {
    lsigmavec = lsigmavec.elem(dupid);
    lxivec = lxivec.elem(dupid);
    lkappavec = lkappavec.elem(dupid);
  }

  {
  arma::vec off0 = Rcpp::as<arma::vec>(offsets[0]);
  if (off0.n_elem > 0) lsigmavec += off0;
  arma::vec off1 = Rcpp::as<arma::vec>(offsets[1]);
  if (off1.n_elem > 0) lxivec += off1;
  arma::vec off2 = Rcpp::as<arma::vec>(offsets[2]);
  if (off2.n_elem > 0) lkappavec += off2;
  }

  double y, lsigma, lxi, lkappa;
  double nllh = 0.0;

  for (int j = 0; j < nobs; j++) {

    y = yvec[j];
    lsigma = lsigmavec[j];
    lxi = lxivec[j];
    lkappa = lkappavec[j];

    double sigma = exp(lsigma);
    double xi = exp(lxi);
    double kappa = exp(lkappa);
    double sk = sqrt(kappa);

    // GPD CDF at y and y+1
    double v_lo = xi * y / sigma;
    double t_lo = 1.0 + v_lo;
    double v_hi = xi * (y + 1.0) / sigma;
    double t_hi = 1.0 + v_hi;

    if (t_lo <= 0.0 || t_hi <= 0.0) {
      nllh = 1e20;
      break;
    }

    double F_lo = 1.0 - R_pow(t_lo, -1.0 / xi);
    double F_hi = 1.0 - R_pow(t_hi, -1.0 / xi);

    // Truncated normal G-function
    double Fmin = R::pnorm(-sk, 0.0, 1.0, 1, 0);
    double denom = 0.5 - Fmin;
    if (denom < 1e-300) denom = 1e-300;
    double G_lo = (R::pnorm(sk * (F_lo - 1.0), 0.0, 1.0, 1, 0) - Fmin) / denom;
    double G_hi = (R::pnorm(sk * (F_hi - 1.0), 0.0, 1.0, 1, 0) - Fmin) / denom;

    double pmf = G_hi - G_lo;
    if (pmf <= 0.0) pmf = 1e-20;

    nllh += -log(pmf);
  }
  return nllh;
}


// //' @rdname degpd5d0
// [[Rcpp::export]]
arma::mat degpd5d12(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::vec yvec, const arma::uvec dupid, int dcate, const Rcpp::List& offsets)
{

  arma::vec lsigmavec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lxivec = X2 * Rcpp::as<arma::vec>(pars[1]);
  arma::vec lkappavec = X3 * Rcpp::as<arma::vec>(pars[2]);
  int nobs = yvec.size();
  arma::mat out = arma::mat(nobs, 9);

  if (dcate == 1) {
    lsigmavec = lsigmavec.elem(dupid);
    lxivec = lxivec.elem(dupid);
    lkappavec = lkappavec.elem(dupid);
  }

  {
  arma::vec off0 = Rcpp::as<arma::vec>(offsets[0]);
  if (off0.n_elem > 0) lsigmavec += off0;
  arma::vec off1 = Rcpp::as<arma::vec>(offsets[1]);
  if (off1.n_elem > 0) lxivec += off1;
  arma::vec off2 = Rcpp::as<arma::vec>(offsets[2]);
  if (off2.n_elem > 0) lkappavec += off2;
  }

  double y, lsigma, lxi, lkappa;

  // Per-observation NLL helper for model 5 (truncated normal)
  auto nll5 = [](double y, double lsigma, double lxi, double lkappa) -> double {
    double sigma = exp(lsigma);
    double xi = exp(lxi);
    double kappa = exp(lkappa);
    double sk = sqrt(kappa);
    double v_lo = xi * y / sigma;
    double t_lo = 1.0 + v_lo;
    double v_hi = xi * (y + 1.0) / sigma;
    double t_hi = 1.0 + v_hi;
    if (t_lo <= 0.0 || t_hi <= 0.0) return 1e20;
    double F_lo = 1.0 - R_pow(t_lo, -1.0 / xi);
    double F_hi = 1.0 - R_pow(t_hi, -1.0 / xi);
    double Fmin = R::pnorm(-sk, 0.0, 1.0, 1, 0);
    double denom = 0.5 - Fmin;
    if (denom < 1e-300) denom = 1e-300;
    double G_lo = (R::pnorm(sk * (F_lo - 1.0), 0.0, 1.0, 1, 0) - Fmin) / denom;
    double G_hi = (R::pnorm(sk * (F_hi - 1.0), 0.0, 1.0, 1, 0) - Fmin) / denom;
    double pmf = G_hi - G_lo;
    if (pmf <= 0.0) pmf = 1e-20;
    return -log(pmf);
  };

  for (int j = 0; j < nobs; j++) {

    y = yvec[j];
    lsigma = lsigmavec[j];
    lxi = lxivec[j];
    lkappa = lkappavec[j];

    double params[3] = {lsigma, lxi, lkappa};
    double grads[3];

    // Central difference gradients
    for (int k = 0; k < 3; k++) {
      double h = std::max(std::abs(params[k]) * 1e-5, 1e-8);
      double pp[3] = {params[0], params[1], params[2]};
      double pm[3] = {params[0], params[1], params[2]};
      pp[k] += h;
      pm[k] -= h;
      grads[k] = (nll5(y, pp[0], pp[1], pp[2]) - nll5(y, pm[0], pm[1], pm[2])) / (2.0 * h);
    }

    out(j, 0) = grads[0];
    out(j, 1) = grads[1];
    out(j, 2) = grads[2];

    // Central difference Hessian (upper triangle)
    int col = 3;
    for (int i = 0; i < 3; i++) {
      for (int k = i; k < 3; k++) {
        double hi = std::max(std::abs(params[i]) * 1e-4, 1e-7);
        double hk = std::max(std::abs(params[k]) * 1e-4, 1e-7);
        double pp[3], pm[3], mp[3], mm[3];
        for (int l = 0; l < 3; l++) {
          pp[l] = params[l];
          pm[l] = params[l];
          mp[l] = params[l];
          mm[l] = params[l];
        }
        pp[i] += hi; pp[k] += hk;
        pm[i] += hi; pm[k] -= hk;
        mp[i] -= hi; mp[k] += hk;
        mm[i] -= hi; mm[k] -= hk;
        out(j, col) = (nll5(y, pp[0], pp[1], pp[2])
                      - nll5(y, pm[0], pm[1], pm[2])
                      - nll5(y, mp[0], mp[1], mp[2])
                      + nll5(y, mm[0], mm[1], mm[2])) / (4.0 * hi * hk);
        col++;
      }
    }
  }
  return out;
}


// //' Discrete Extended generalized Pareto distribution of type 6 (deGPD6) negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for each deGPD parameter
// //' @param X1 a design matrix for the deGPD log scale parameter
// //' @param X2 a design matrix for the deGPD log shape parameter
// //' @param X3 a design matrix for the deGPD log kappa
// //' @param yvec a vector
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return degpd6d0 a scalar, the negative log-liklihood
// //' @return degpd6d12 a matrix, first then second derivatives w.r.t. deGPD6 parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double degpd6d0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2, const arma::mat& X3, arma::vec yvec, const arma::uvec& dupid, int dcate, const Rcpp::List& offsets)
{

  arma::vec lsigmavec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lxivec = X2 * Rcpp::as<arma::vec>(pars[1]);
  arma::vec lkappavec = X3 * Rcpp::as<arma::vec>(pars[2]);
  int nobs = yvec.size();

  if (dcate == 1) {
    lsigmavec = lsigmavec.elem(dupid);
    lxivec = lxivec.elem(dupid);
    lkappavec = lkappavec.elem(dupid);
  }

  {
  arma::vec off0 = Rcpp::as<arma::vec>(offsets[0]);
  if (off0.n_elem > 0) lsigmavec += off0;
  arma::vec off1 = Rcpp::as<arma::vec>(offsets[1]);
  if (off1.n_elem > 0) lxivec += off1;
  arma::vec off2 = Rcpp::as<arma::vec>(offsets[2]);
  if (off2.n_elem > 0) lkappavec += off2;
  }

  double y, lsigma, lxi, lkappa;
  double nllh = 0.0;

  for (int j = 0; j < nobs; j++) {

    y = yvec[j];
    lsigma = lsigmavec[j];
    lxi = lxivec[j];
    lkappa = lkappavec[j];

    double sigma = exp(lsigma);
    double xi = exp(lxi);
    double kappa = exp(lkappa);

    // GPD CDF at y and y+1
    double v_lo = xi * y / sigma;
    double t_lo = 1.0 + v_lo;
    double v_hi = xi * (y + 1.0) / sigma;
    double t_hi = 1.0 + v_hi;

    if (t_lo <= 0.0 || t_hi <= 0.0) {
      nllh = 1e20;
      break;
    }

    double F_lo = 1.0 - R_pow(t_lo, -1.0 / xi);
    double F_hi = 1.0 - R_pow(t_hi, -1.0 / xi);

    // Truncated Beta G-function
    // G(u; kappa) = [pbeta((0.5 - 1/32)*u + 1/32, kappa, kappa) - pbeta(1/32, kappa, kappa)]
    //             / [pbeta(0.5, kappa, kappa) - pbeta(1/32, kappa, kappa)]
    double c1 = 0.5 - 1.0 / 32.0;  // = 15/32
    double c2 = 1.0 / 32.0;         // = 1/32
    double pb_min = R::pbeta(c2, kappa, kappa, 1, 0);
    double pb_half = R::pbeta(0.5, kappa, kappa, 1, 0);
    double denom_b = pb_half - pb_min;
    if (denom_b < 1e-300) denom_b = 1e-300;

    double u_lo = c1 * F_lo + c2;
    double u_hi = c1 * F_hi + c2;
    double G_lo = (R::pbeta(u_lo, kappa, kappa, 1, 0) - pb_min) / denom_b;
    double G_hi = (R::pbeta(u_hi, kappa, kappa, 1, 0) - pb_min) / denom_b;

    double pmf = G_hi - G_lo;
    if (pmf <= 0.0) pmf = 1e-20;

    nllh += -log(pmf);
  }
  return nllh;
}


// //' @rdname degpd6d0
// [[Rcpp::export]]
arma::mat degpd6d12(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::vec yvec, const arma::uvec dupid, int dcate, const Rcpp::List& offsets)
{

  arma::vec lsigmavec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lxivec = X2 * Rcpp::as<arma::vec>(pars[1]);
  arma::vec lkappavec = X3 * Rcpp::as<arma::vec>(pars[2]);
  int nobs = yvec.size();
  arma::mat out = arma::mat(nobs, 9);

  if (dcate == 1) {
    lsigmavec = lsigmavec.elem(dupid);
    lxivec = lxivec.elem(dupid);
    lkappavec = lkappavec.elem(dupid);
  }

  {
  arma::vec off0 = Rcpp::as<arma::vec>(offsets[0]);
  if (off0.n_elem > 0) lsigmavec += off0;
  arma::vec off1 = Rcpp::as<arma::vec>(offsets[1]);
  if (off1.n_elem > 0) lxivec += off1;
  arma::vec off2 = Rcpp::as<arma::vec>(offsets[2]);
  if (off2.n_elem > 0) lkappavec += off2;
  }

  double y, lsigma, lxi, lkappa;

  // Per-observation NLL helper for model 6 (truncated beta)
  auto nll6 = [](double y, double lsigma, double lxi, double lkappa) -> double {
    double sigma = exp(lsigma);
    double xi = exp(lxi);
    double kappa = exp(lkappa);
    double v_lo = xi * y / sigma;
    double t_lo = 1.0 + v_lo;
    double v_hi = xi * (y + 1.0) / sigma;
    double t_hi = 1.0 + v_hi;
    if (t_lo <= 0.0 || t_hi <= 0.0) return 1e20;
    double F_lo = 1.0 - R_pow(t_lo, -1.0 / xi);
    double F_hi = 1.0 - R_pow(t_hi, -1.0 / xi);
    double c1 = 0.5 - 1.0 / 32.0;
    double c2 = 1.0 / 32.0;
    double pb_min = R::pbeta(c2, kappa, kappa, 1, 0);
    double pb_half = R::pbeta(0.5, kappa, kappa, 1, 0);
    double denom_b = pb_half - pb_min;
    if (denom_b < 1e-300) denom_b = 1e-300;
    double u_lo = c1 * F_lo + c2;
    double u_hi = c1 * F_hi + c2;
    double G_lo = (R::pbeta(u_lo, kappa, kappa, 1, 0) - pb_min) / denom_b;
    double G_hi = (R::pbeta(u_hi, kappa, kappa, 1, 0) - pb_min) / denom_b;
    double pmf = G_hi - G_lo;
    if (pmf <= 0.0) pmf = 1e-20;
    return -log(pmf);
  };

  for (int j = 0; j < nobs; j++) {

    y = yvec[j];
    lsigma = lsigmavec[j];
    lxi = lxivec[j];
    lkappa = lkappavec[j];

    double params[3] = {lsigma, lxi, lkappa};
    double grads[3];

    // Central difference gradients
    for (int k = 0; k < 3; k++) {
      double h = std::max(std::abs(params[k]) * 1e-5, 1e-8);
      double pp[3] = {params[0], params[1], params[2]};
      double pm[3] = {params[0], params[1], params[2]};
      pp[k] += h;
      pm[k] -= h;
      grads[k] = (nll6(y, pp[0], pp[1], pp[2]) - nll6(y, pm[0], pm[1], pm[2])) / (2.0 * h);
    }

    out(j, 0) = grads[0];
    out(j, 1) = grads[1];
    out(j, 2) = grads[2];

    // Central difference Hessian (upper triangle)
    int col = 3;
    for (int i = 0; i < 3; i++) {
      for (int k = i; k < 3; k++) {
        double hi = std::max(std::abs(params[i]) * 1e-4, 1e-7);
        double hk = std::max(std::abs(params[k]) * 1e-4, 1e-7);
        double pp[3], pm[3], mp[3], mm[3];
        for (int l = 0; l < 3; l++) {
          pp[l] = params[l];
          pm[l] = params[l];
          mp[l] = params[l];
          mm[l] = params[l];
        }
        pp[i] += hi; pp[k] += hk;
        pm[i] += hi; pm[k] -= hk;
        mp[i] -= hi; mp[k] += hk;
        mm[i] -= hi; mm[k] -= hk;
        out(j, col) = (nll6(y, pp[0], pp[1], pp[2])
                      - nll6(y, pm[0], pm[1], pm[2])
                      - nll6(y, mp[0], mp[1], mp[2])
                      + nll6(y, mm[0], mm[1], mm[2])) / (4.0 * hi * hk);
        col++;
      }
    }
  }
  return out;
}
