#include <Rcpp.h>

using namespace Rcpp;

// This file contains the following functions:

// General functions for most models
// tempupdate - updates the geometrically spaced temperatures
// matcomp - computes proposal values for the beta's
// offsetcompute - computes the part of the offset containing the trends information
// matN - takes a vector and creates a matrix across number of time points/number of chains
// poissondevfit - computes deviance for poisson
// linpredcompute - computing the linear predictor for covariates
// poissonbetablockupdate - for updating covariate effects based on a poisson likelihood
// gammaproposal - for proposing new candidate values for the gamma's
// poissongammaupdate - updates the gamma parameters for the poisson
// poissonwupdate - updates the trend indicator variables (omega's) for the poisson
// lambdaupdate - updates the lambda parameters
// poissonphiupdate - updates the spatial random effects phi for the poisson
// tau2quadform - for computing the sum of quadratic forms for updating tau2
// tau2compute - computes the full conditional of tau2 and updates tau2
// quadform - for computing the sum of quadratic forms
// rhoquadformcompute - for computing the sum of quadratic forms for updating rho
// Qdet - computes determinant of precision matrix Q
// poissoncouplingAllupdate - performs the coupling step for mixing chains for the poisson


// [[Rcpp::export]]
NumericVector tempupdate(const int Nchains, double dt)
{
  //Create new objects
  NumericVector newtemps(Nchains);
  
  newtemps[0] = 1;
  
  for(int i = 1; i < Nchains; i++)
  {
    newtemps[i] = dt * newtemps[(i-1)];
  }
  // Return the result
  return newtemps;
}


// [[Rcpp::export]]
NumericMatrix matcomp(NumericMatrix X, NumericMatrix beta, NumericVector prop, const int p, const int Nchains)
{
  //Create new objects
  NumericVector newprop = clone(prop);
  NumericMatrix proposal(p, Nchains), newbeta = clone(beta), newX = clone(X);
  NumericVector gen(p);
  NumericVector matmult(p);
  
  Environment base("package:stats");
  Function sample = base["rnorm"];
  
  for(int i = 0; i < Nchains; i++)
  {
    gen = as<NumericVector>( rnorm(p, 0, 1) );
    
    for(int j = 0; j < p; j++)
    {
      matmult[j] = sum((sqrt(newprop[i]) * newX(j, _)) * gen);
    }
    proposal(_, i) = newbeta(_, i) +  matmult;
  }
  // Return the result
  return proposal;
}


// [[Rcpp::export]]
NumericMatrix offsetcompute(NumericMatrix w, NumericMatrix gamma, NumericMatrix time, const int Nchains, const int nsites, const int Ntrends, NumericVector begin)
{
  //Create new objects
  NumericMatrix wnew = clone(w), gammanew = clone(gamma), timenew = clone(time);
  NumericMatrix offset(nsites, Nchains);
  int rowstart;
  
  for(int i = 0; i < nsites; i++)
  {
    for(int j = 0; j < Ntrends; j++)
    {
      rowstart = begin[j] - 1;
      
      offset(i, _) = offset(i, _) + wnew((rowstart + i), _) * (gammanew((rowstart + i), _) * timenew((rowstart + i), _));
    }
    
  }
  
  return offset;
}


// [[Rcpp::export]]
NumericMatrix matN(NumericVector x, const int nsites, const int Nchains)
{
  //Create new objects
  NumericMatrix mat(nsites, Nchains);
  NumericVector newx = clone(x);
  
  for(int j = 0; j < nsites; j++)
  {
    mat(j, _) = newx; 
  }
  // Return the result
  return mat;
}  


// [[Rcpp::export]]
List poissondevfit(NumericVector y, NumericMatrix fitted, const int nsites, const int Nchains)
{
  NumericVector newy = clone(y);
  NumericMatrix newfitted = clone(fitted), like_all(nsites, Nchains);
  NumericVector deviance(Nchains), deviance_all(nsites), fit(nsites);
  
  Environment base2("package:stats");
  Function dpois = base2["dpois"];
  
  for(int j = 0; j < Nchains; j++)
  {
    fit = newfitted(_, j);
    deviance_all = as<NumericVector>( dpois(newy, fit) );
    like_all(_, j) = deviance_all;
    deviance[j] = -2 * sum(log(deviance_all));
  }
  
  List out(2);
  out[0] = deviance;
  out[1] = like_all;
  return out;
}


// [[Rcpp::export]]
NumericMatrix linpredcompute(NumericMatrix X, const int nsites, const int p, NumericMatrix beta, const int Nchains)
{
  // Create new objects
  // Compute the linear predictor
  NumericMatrix linpred(nsites, Nchains);
  double temp;
  
  for(int j = 0; j < Nchains; j++)
  {
    // Compute the linear predictor via a double for loop
    for(int k = 0; k < nsites; k++)
    {
      temp = 0;
      
      for(int l = 0; l < p; l++) temp = temp + X(k, l) * beta(l, j);
      
      linpred(k, j) = temp;
    }
  }
  
  // Return the result
  return linpred;
}


// [[Rcpp::export]]
NumericVector poissonbetablockupdate(const int nsites, NumericMatrix beta, NumericMatrix betaprop, NumericMatrix lp_beta, NumericMatrix lp_betaprop,
                                     NumericMatrix offset, NumericVector y, NumericVector prior_meanbeta, NumericVector prior_varbeta,
                                     const int Nchains, NumericVector temps, const int p)
{
  // Compute the acceptance probability for beta
  //Create new objects
  double oldlikebit=0, newlikebit=0, likebit, priorbit=0;
  NumericVector lp_current(nsites), lp_proposal(nsites), p_current(nsites), p_proposal(nsites), acceptance(Nchains);
  NumericMatrix betanew = clone(beta), betapropnew = clone(betaprop);
  NumericMatrix newlp_beta = clone(lp_beta), newlp_betaprop = clone(lp_betaprop), newoffset = clone(offset);
  
  for(int j = 0; j < Nchains; j++)
  {
    oldlikebit = 0;
    newlikebit = 0;
    priorbit = 0;
    
    for(int i = 0; i < nsites; i++)
    {
      lp_current[i] = newlp_beta(i, j) + newoffset(i, j);
      lp_proposal[i] = newlp_betaprop(i, j) + newoffset(i, j);
      
      p_current[i] = exp(lp_current[i]);
      p_proposal[i] = exp(lp_proposal[i]);
      
      oldlikebit += y[i] * lp_current[i] - p_current[i];
      newlikebit += y[i] * lp_proposal[i] - p_proposal[i];
    }
    likebit = newlikebit - oldlikebit;
    
    // Create the prior acceptance component
    for(int k = 0; k < p; k++)
    {
      priorbit += 0.5 * pow((betanew(k, j) - prior_meanbeta[k]), 2) / prior_varbeta[k] - 0.5 * pow((betapropnew(k, j) - prior_meanbeta[k]), 2) / prior_varbeta[k];
    }
    
    // Compute the acceptance probability and return the value
    acceptance[j] = exp((likebit + priorbit) * temps[j]);
  }
  
  return acceptance;
}


// [[Rcpp::export]]
NumericVector gammaproposal(const int Nchains, NumericVector gamma, NumericVector gamma_tune, const int prior_vargamma, NumericVector Wareas, const int trend, const int knots)
{
  NumericVector proposal(Nchains), gammanew = clone(gamma);
  NumericVector Wareasnew = clone(Wareas);
  double inf = std::numeric_limits<double>::infinity();
  
  Environment base("package:truncdist");
  Function rtrunc = base["rtrunc"];
  
  for(int j = 0; j < Nchains; j++)
  {
    if(Wareasnew[j] == 0)
    {
      
      if(trend == 2 | trend == 5 | trend == 6 | (trend >= 8 & trend <= (8 + knots)))
      {
        proposal[j] = as<double>( rtrunc(1, "norm", -inf, 0, 0, 0.01) );
      }else if(trend == 3 | trend == 4 | trend == 7 | (trend >= (8 + knots + 1) & trend <= (8 + knots + 1 + knots)))
      {
        proposal[j] = as<double>( rtrunc(1, "norm", 0, inf, 0, 0.01) );
      }else
      {}
    }
    else
    {
      if(trend == 2 | trend == 5 | trend == 6 | (trend >= 8 & trend <= (8 + knots)))
      {
        proposal[j] = as<double>( rtrunc(1, "norm", -inf, 0, gammanew[j], gamma_tune[j]) );
      }else if(trend == 3 | trend == 4 | trend == 7 | (trend >= (8 + knots + 1) & trend <= (8 + knots + 1 + knots)))
      {
        proposal[j] = as<double>( rtrunc(1, "norm", 0, inf, gammanew[j], gamma_tune[j]) );
      }else
      {}
    }
  }
  return proposal;
}


// [[Rcpp::export]]
List poissongammaupdate(const int nsites, NumericVector gamma, NumericVector proposal, NumericMatrix offset, NumericMatrix offset_proposal, NumericVector y,
                        double prior_meangamma, double prior_vargamma, const int Nchains, NumericVector temps)
{
  // Compute the acceptance probability for gamma
  //Create new objects
  double acceptance, oldlikebit=0, newlikebit=0, likebit, priorbit;
  NumericVector lp_current(nsites), lp_proposal(nsites), p_current(nsites), p_proposal(nsites);
  NumericVector accept(Nchains);
  NumericVector gammanew = clone(gamma), proposalnew = clone(proposal);
  NumericMatrix newoffset = clone(offset), newproposal = clone(offset_proposal);
  
  for(int j = 0; j < Nchains; j++)
  {
    oldlikebit = 0;
    newlikebit = 0;
    
    for(int i = 0; i < nsites; i++)
    {
      
      lp_current[i] = newoffset(i, j);
      lp_proposal[i] = newproposal(i, j);
      
      p_current[i] = exp(lp_current[i]);
      p_proposal[i] = exp(lp_proposal[i]);
      
      oldlikebit += y[i] * lp_current[i] - p_current[i];
      newlikebit += y[i] * lp_proposal[i] - p_proposal[i];
    }
    likebit = newlikebit - oldlikebit;
    
    priorbit = 0.5 * pow((gammanew[j] - prior_meangamma), 2) / prior_vargamma - 0.5 * pow((proposalnew[j] - prior_meangamma), 2) / prior_vargamma;
    
    // Compute the acceptance probability and return the value
    acceptance = exp((likebit + priorbit) * temps[j]);
    if(runif(1)[0] <= acceptance)
    {
      gammanew[j] = proposalnew[j];
      accept[j] = accept[j] + 1;
    }
    else
    {
    }
  }
  
  List out(2);
  out[0] = gammanew;
  out[1] = accept;
  return out;
}


// [[Rcpp::export]]
List poissonwupdate(const int nsites, const int ntimes, NumericMatrix w, NumericMatrix offset, NumericMatrix offset_proposal, NumericMatrix w_proposal,
                    NumericMatrix y, NumericMatrix lambda, const int Nchains, NumericVector temps, NumericVector begin, NumericVector regbegin, const int Ntrends)
{
  // Compute the acceptance probability for beta
  //Create new objects
  double acceptance, oldlikebit=0, newlikebit=0, likebit, priorbit;
  NumericVector lp_current(ntimes), lp_proposal(ntimes), p_current(ntimes), p_proposal(ntimes);
  NumericMatrix accept(nsites, Nchains);
  NumericMatrix wnew = clone(w), wproposal = clone(w_proposal), lambdanew = clone(lambda), newoffset = clone(offset), newproposal = clone(offset_proposal);
  int rowstart, regstart;
  NumericVector proposal;
  
  for(int j = 0; j < Nchains; j++)
  {
    rowstart = begin[j] - 1;
    
    for(int k = 0; k < nsites; k++)
    {
      proposal = wproposal((rowstart + k), _);
      
      oldlikebit = 0;
      newlikebit = 0;
      priorbit = 0;
      
      for(int t = 0; t < ntimes; t++)
      {
        regstart = regbegin[t] - 1;
        
        lp_current[t] = newoffset((regstart + k), j);
        lp_proposal[t] = newproposal((regstart + k), j);
        
        p_current[t] = exp(lp_current[t]);
        p_proposal[t] = exp(lp_proposal[t]);
        
        oldlikebit += y(k, t) * lp_current[t] - p_current[t];
        newlikebit += y(k, t) * lp_proposal[t] - p_proposal[t];
      }
      
      likebit = newlikebit - oldlikebit;
      
      for(int i = 0; i < Ntrends; i++)
      {
        priorbit += proposal[i] * log(lambdanew(i, j)) - wnew((rowstart + k), i) * log(lambdanew(i, j));
      }
      
      // Compute the acceptance probability and return the value
      acceptance = exp((likebit + priorbit) * temps[j]);
      if(runif(1)[0] <= acceptance) 
      {
        wnew((rowstart + k), _) = proposal;
        accept(k, j) = accept(k, j) + 1;
      }
      else
      { 
      }
    }
  }
  List out(2);
  out[0] = wnew;
  out[1] = accept;
  return out;
}


// [[Rcpp::export]]
NumericMatrix lambdaupdate(const int Nchains, NumericMatrix temp)
{
  Environment base("package:LearnBayes");
  Function rdirichlet = base["rdirichlet"];
  //Create new objects
  NumericMatrix tempnew = clone(temp), lambdanew = clone(temp);
  
  for(int j = 0; j < Nchains; j++)
  {
    lambdanew(_, j) = as<NumericVector>( rdirichlet(1, tempnew(_, j)) );
  }
  return lambdanew;
}


// [[Rcpp::export]]
List poissonphiupdate(NumericMatrix Wtriplet, NumericMatrix Wbegfin, NumericVector Wtripletsum, const int nsites, const int ntimes,
                      NumericMatrix phi, NumericMatrix offset, NumericMatrix y, NumericVector tau2, NumericVector rho, const int Nchains,
                      NumericVector temps, NumericMatrix phi_tune, NumericVector regbegin)
{
  // Compute the acceptance probability
  //Create new objects
  double acceptance, oldlikebit=0, newlikebit=0, likebit, priorbit;
  NumericVector lp_current(ntimes), lp_proposal(ntimes), p_current(ntimes), p_proposal(ntimes);
  NumericVector rhonew = clone(rho), tau2new = clone(tau2);
  NumericMatrix accept(nsites, Nchains);
  NumericMatrix phinew = clone(phi), newoffset = clone(offset);
  double proposal;
  double priorvardenom, priormean, priorvar, sumphi;
  int rowstart=0, rowend=0;
  int regstart;
  
  for(int j = 0; j < Nchains; j++)
  {
    
    for(int k = 0; k < nsites; k++)
    {
      
      oldlikebit = 0;
      newlikebit = 0;
      priorbit = 0;
      
      // Calculate prior variance
      priorvardenom = rhonew[j] * Wtripletsum[k] + 1 - rhonew[j];
      priorvar = tau2new[j] / priorvardenom;
      
      // Calculate the prior mean
      rowstart = Wbegfin(k, 0) - 1;
      rowend = Wbegfin(k, 1);
      sumphi = 0;
      for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phinew((Wtriplet(l, 1) - 1), j);
      priormean = rhonew[j] * sumphi / priorvardenom; 
      
      // propose a value 
      proposal = rnorm(1, phinew(k, j), sqrt(priorvar*phi_tune(k, j)))[0];
      
      for(int t = 0; t < ntimes; t++)
      {
        regstart = regbegin[t] - 1;
        
        lp_current[t] = newoffset((regstart + k), j) + phinew(k, j);
        lp_proposal[t] = newoffset((regstart + k), j) + proposal;
        
        p_current[t] = exp(lp_current[t]);
        p_proposal[t] = exp(lp_proposal[t]);
        
        oldlikebit += y(k, t) * lp_current[t] - p_current[t];
        newlikebit += y(k, t) * lp_proposal[t] - p_proposal[t];
      }
      
      likebit = newlikebit - oldlikebit;
      priorbit = (0.5/priorvar) * pow((phinew(k, j) - priormean), 2) - (0.5/priorvar) * pow((proposal - priormean), 2);
      
      // Compute the acceptance probability and return the value
      acceptance = exp((likebit + priorbit) * temps[j]);
      if(runif(1)[0] <= acceptance) 
      {
        phinew(k, j) = proposal;
        accept(k, j) = accept(k, j) + 1;
      }
      else
      { 
      }
    }
  }
  List out(2);
  out[0] = phinew;
  out[1] = accept;
  return out;
}


// [[Rcpp::export]]
NumericVector tau2quadform(NumericMatrix Wtriplet, NumericVector Wtripletsum, const int n_triplet, const int nsites, 
                           NumericMatrix phi, NumericMatrix theta, NumericVector rho, const int Nchains)
{
  // Compute a quadratic form for the random effects
  // Create new objects 
  NumericVector tau2_posteriorscale(Nchains), rhonew = clone(rho);
  NumericMatrix phinew = clone(phi), thetanew = clone(theta);
  double tau2_quadform = 0, tau2_phisq = 0;
  int row, col;
  
  for(int j = 0; j < Nchains; j++)
  {
    tau2_quadform = 0;
    tau2_phisq = 0;
    // Compute the off diagonal elements of the quadratic form
    for(int i = 0; i < n_triplet; i++)
    {
      row = Wtriplet(i, 0) - 1;
      col = Wtriplet(i, 1) - 1;
      tau2_quadform += phinew((Wtriplet(i, 0) - 1), j) * thetanew((Wtriplet(i, 1) - 1), j) * Wtriplet(i, 2); 
    }
    // Compute the diagonal elements of the quadratic form          
    for(int l = 0; l < nsites; l++)
    {
      tau2_phisq += phinew(l, j) * thetanew(l, j) * (rhonew[j] * Wtripletsum[l] + 1 - rhonew[j]);    
    }
    // Compute the quadratic form
    tau2_posteriorscale[j] = 0.5 * (tau2_phisq - rhonew[j] * tau2_quadform);
  }
  // Return the simulated values
  return tau2_posteriorscale;
}


// [[Rcpp::export]]
NumericVector tau2compute(NumericVector temp, const double tau2_shape, const double prior_tau2, const int Nchains)
{
  NumericVector tau2new(Nchains);
  double tau2_scale;
  
  for (int j = 0; j < Nchains; j++)
  {
    tau2_scale = temp[j] + prior_tau2;
    tau2new[j] = 1 / rgamma(1, tau2_shape, (1/tau2_scale))[0];
  }
  return tau2new;
}


// [[Rcpp::export]]
double quadform(NumericMatrix Wtriplet, NumericVector Wtripletsum, const int n_triplet, const int nsites, 
                NumericVector phi, NumericVector delta, double rho)
{
  // Compute a quadratic form for the random effects
  // Create new objects 
  double rho_quad;
  double rho_quadform = 0, rho_phisq = 0;
  int row, col;
  
  // Compute the off diagonal elements of the quadratic form
  for(int l = 0; l < n_triplet; l++)
  {
    row = Wtriplet(l, 0) - 1;
    col = Wtriplet(l, 1) - 1;
    rho_quadform += phi[(Wtriplet(l, 0) - 1)] * delta[(Wtriplet(l, 1) - 1)] * Wtriplet(l, 2); 
  }
  
  // Compute the diagonal elements of the quadratic form          
  for(int l = 0; l < nsites; l++)
  {
    rho_phisq += phi[l] * delta[l] * (rho * Wtripletsum[l] + 1 - rho);    
  }
  
  // Compute the quadratic form
  rho_quad = 0.5 * (rho_phisq - rho * rho_quadform);
  
  // Return the simulated value
  return rho_quad;
}


// [[Rcpp::export]]
NumericVector rhoquadformcompute(NumericMatrix Wtriplet, NumericVector Wtripletsum, const int n_triplet, 
                                 const int nsites, const int Nchains, NumericMatrix phi, NumericVector rho, NumericVector tau2)
{    
  NumericVector temp(nsites);
  NumericVector rho_quadform(Nchains);
  NumericVector rhonew = clone(rho), tau2new = clone(tau2);
  NumericMatrix phinew = clone(phi);
  
  for(int j = 0; j < Nchains; j++)
  {
    temp = phinew(_, j);  
    rho_quadform[j] = quadform(Wtriplet, Wtripletsum, n_triplet, nsites, temp, temp, rhonew[j]) / tau2new[j];
  }
  return rho_quadform;
}


// [[Rcpp::export]]
NumericVector Qdet(const int Nchains, NumericVector rho, NumericVector Wstar_val)
{
  NumericVector detQ(Nchains);
  NumericVector rhonew = clone(rho);
  
  for(int j = 0; j < Nchains; j++)
  {
    detQ[j] = 0.5 * sum(log((rhonew[j] * Wstar_val + (1 - rhonew[j]))));
  }
  return detQ;
}


// [[Rcpp::export]]
int poissoncouplingAllupdate(const int nsites, const int K, const int p, NumericMatrix w, NumericMatrix offset, NumericMatrix beta, NumericMatrix gamma,
                             NumericMatrix lambda, NumericMatrix phi, NumericVector rho, NumericVector tau2, NumericVector Wtripletsum,
                             NumericMatrix Wtriplet, NumericMatrix Wbegfin, NumericVector y, NumericVector prior_meanbeta, NumericVector prior_varbeta,
                             NumericVector prior_meantrends, NumericVector prior_vartrends, NumericVector prior_lambda, NumericVector prior_tau2,
                             NumericVector swap, NumericVector temps, NumericVector begin, const int Ntrends, const int TrendSel)
{
  // Compute the acceptance probability for swapping
  //Create new objects
  double acceptance, chain1oldlikebit=0, chain2oldlikebit=0, chain1newlikebit=0, chain2newlikebit=0, chain1likebit, chain2likebit, chain1priorbitbeta=0, chain2priorbitbeta=0;
  double wpriorchain1old=0, wpriorchain1new=0, wpriorchain2old=0, wpriorchain2new=0, chain1priorbitw, chain2priorbitw, chain1priorbitlambda, chain2priorbitlambda;
  double phichain1=0, phichain2=0;
  double chain1priorbitgamma=0, chain2priorbitgamma=0;
  double chain1priorbittau2, chain2priorbittau2, chain1priorbitphi, chain2priorbitphi;
  NumericVector lpchain1_current(nsites), lpchain1_proposal(nsites), pchain1_current(nsites), pchain1_proposal(nsites);
  NumericVector lpchain2_current(nsites), lpchain2_proposal(nsites), pchain2_current(nsites), pchain2_proposal(nsites);
  int accept = 0;
  NumericVector swapnew;
  swapnew = swap - 1;
  double swap1, swap2;
  swap1 = swapnew[0];
  swap2 = swapnew[1];
  NumericVector rhonew = clone(rho), tau2new = clone(tau2);
  NumericMatrix gammanew = clone(gamma), wnew = clone(w), newoffset = clone(offset), betanew = clone(beta), lambdanew = clone(lambda), phinew = clone(phi);
  int beginchain1 = begin[swap1] - 1;
  int beginchain2 = begin[swap2] - 1;
  double temp1, temp2;
  temp1 = temps[swap1];
  temp2 = temps[swap2];
  double priorvardenomchain1, priorvardenomchain2, phipriorvarchain1, phipriorvarchain2, sumphichain1=0, sumphichain2=0, phipriormeanchain1, phipriormeanchain2;
  int rowstart=0, rowend=0;
  
  Environment base("package:gtools");
  Function ddirichlet = base["ddirichlet"];
  
  for(int i = 0; i < nsites; i++)
  {
    
    lpchain1_current[i] = newoffset(i, swap1);
    pchain1_current[i] = exp(lpchain1_current[i]);
    
    lpchain2_current[i] = newoffset(i, swap2);
    pchain2_current[i] = exp(lpchain2_current[i]);
    
    lpchain1_proposal[i] = newoffset(i, swap2);
    pchain1_proposal[i] = exp(lpchain1_proposal[i]);
    
    lpchain2_proposal[i] = newoffset(i, swap1);
    pchain2_proposal[i] = exp(lpchain2_proposal[i]);
    
    chain1oldlikebit += y[i] * lpchain1_current[i] - pchain1_current[i];
    chain2oldlikebit += y[i] * lpchain2_current[i] - pchain2_current[i];
    
    chain1newlikebit += y[i] * lpchain1_proposal[i] - pchain1_proposal[i];
    chain2newlikebit += y[i] * lpchain2_proposal[i] - pchain2_proposal[i];
  }
  
  chain1likebit = chain1newlikebit - chain1oldlikebit;
  chain2likebit = chain2newlikebit - chain2oldlikebit;
  
  for(int j = 0; j < K; j++)
  {
    
    for(int l = 0; l < Ntrends; l++)
    {
      wpriorchain1old += wnew((beginchain1 + j), l) * log(lambdanew(l, swap1));
      wpriorchain1new += wnew((beginchain2 + j), l) * log(lambdanew(l, swap1));
      
      wpriorchain2old += wnew((beginchain2 + j), l) * log(lambdanew(l, swap2));
      wpriorchain2new += wnew((beginchain1 + j), l) * log(lambdanew(l, swap2));
    }
    
    priorvardenomchain1 = rhonew[swap1] * Wtripletsum[j] + 1 - rhonew[swap1];
    phipriorvarchain1 = tau2new[swap1] / priorvardenomchain1;
    
    priorvardenomchain2 = rhonew[swap2] * Wtripletsum[j] + 1 - rhonew[swap2];
    phipriorvarchain2 = tau2new[swap2] / priorvardenomchain2;
    
    rowstart = Wbegfin(j, 0) - 1;
    rowend = Wbegfin(j, 1);
    for(int l = rowstart; l < rowend; l++)
    {
      sumphichain1 += Wtriplet(l, 2) * phinew((Wtriplet(l, 1) - 1), swap1);
      sumphichain2 += Wtriplet(l, 2) * phinew((Wtriplet(l, 1) - 1), swap2);
    }
    phipriormeanchain1 = rhonew[swap1] * sumphichain1 / priorvardenomchain1; 
    phipriormeanchain2 = rhonew[swap2] * sumphichain2 / priorvardenomchain2; 
    
    phichain1 += (0.5/phipriorvarchain1) * pow((phinew(j, swap1) - phipriormeanchain1), 2) - (0.5/phipriorvarchain1) * pow((phinew(j, swap2) - phipriormeanchain1), 2);
    phichain2 += (0.5/phipriormeanchain2) * pow((phinew(j, swap2) - phipriormeanchain2), 2) - (0.5/phipriormeanchain2) * pow((phinew(j, swap1) - phipriormeanchain2), 2);
  }
  
  for(int k = 0; k < p; k++)     
  {
    chain1priorbitbeta += 0.5 * pow((betanew(k, swap1) - prior_meanbeta[k]), 2) / prior_varbeta[k] - 0.5 * pow((betanew(k, swap2) - prior_meanbeta[k]), 2) / prior_varbeta[k];
    chain2priorbitbeta += 0.5 * pow((betanew(k, swap2) - prior_meanbeta[k]), 2) / prior_varbeta[k] - 0.5 * pow((betanew(k, swap1) - prior_meanbeta[k]), 2) / prior_varbeta[k];
  }
  
  for(int l = 1; l < TrendSel; l++)
  {
    chain1priorbitgamma +=  0.5 * pow((gammanew(l, swap1) - prior_meantrends[l]), 2) / prior_vartrends[l] - 0.5 * pow((gammanew(l, swap2) - prior_meantrends[l]), 2) / prior_vartrends[l];
    chain2priorbitgamma +=  0.5 * pow((gammanew(l, swap2) - prior_meantrends[l]), 2) / prior_vartrends[l] - 0.5 * pow((gammanew(l, swap1) - prior_meantrends[l]), 2) / prior_vartrends[l];
  }
  
  chain1priorbitw = wpriorchain1new - wpriorchain1old;  
  chain2priorbitw = wpriorchain2new - wpriorchain2old;  
  
  chain1priorbitphi = phichain1;  
  chain2priorbitphi = phichain2; 
  
  chain1priorbitlambda = log(as<double>( ddirichlet(lambdanew(_, swap1), prior_lambda) )) - log(as<double>( ddirichlet(lambdanew(_, swap2), prior_lambda) ));
  chain2priorbitlambda = log(as<double>( ddirichlet(lambdanew(_, swap2), prior_lambda) )) - log(as<double>( ddirichlet(lambdanew(_, swap1), prior_lambda) ));
  
  chain1priorbittau2 = (-prior_tau2[0] - 1) * log(tau2[swap1]) - (prior_tau2[1] / tau2[swap1]) - (-prior_tau2[0] - 1) * log(tau2[swap2]) - (prior_tau2[1] / tau2[swap2]);
  chain2priorbittau2 = (-prior_tau2[0] - 1) * log(tau2[swap2]) - (prior_tau2[1] / tau2[swap2]) - (-prior_tau2[0] - 1) * log(tau2[swap1]) - (prior_tau2[1] / tau2[swap1]);
  
  // Compute the acceptance probability and return the value
  acceptance = exp(((chain1likebit + chain1priorbitbeta + chain1priorbitgamma + chain1priorbitw + chain1priorbitlambda + chain1priorbitphi + chain1priorbittau2) * temps[swap1]) + ((chain2likebit + chain2priorbitbeta + chain2priorbitgamma + chain2priorbitw + chain2priorbitlambda + chain2priorbitphi + chain2priorbittau2) * temps[swap2]));
  
  if(runif(1)[0] <= acceptance) 
  {
    accept = accept + 1;
  }
  else
  { 
  }
  return accept;
}
