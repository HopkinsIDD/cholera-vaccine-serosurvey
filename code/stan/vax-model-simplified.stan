
data {
  int  S0Q0;
  int  S1Q0;
  int  S0Q1;
  int  S1Q1;

  // real gamma_p1;
  // real gamma_p2;
  
  real alpha_p1;
  real alpha_p2;
  real beta_p1;
  real beta_p2;
  real epsilon_p1;
  real epsilon_p2;
  
}

parameters {
  real<lower=0,upper=1> lambda; // seroincidence
  real<lower=0,upper=1> gamma; //coverage

  real<lower=0,upper=1> alpha; //sens
  real<lower=0,upper=1> beta; //spec among unvaxxed
  real<lower=0,upper=1> epsilon; //spec among vaxxed

}




model {
        
  // target += beta_lpdf(gamma|gamma_p1,gamma_p2) ; //coverage prior
// 
  target += beta_lpdf(alpha|alpha_p1,alpha_p2) ; //sensitivity prior
  target += beta_lpdf(beta|beta_p1,beta_p2); //specificity prior on unvaxxed
  target += beta_lpdf(epsilon| epsilon_p1, epsilon_p2) ; //specificity prior on vaxxed
//   
  // target += normal_lpdf(alpha|alpha_p1,alpha_p2) ; //sensitivity prior
  // target += normal_lpdf(beta|beta_p1,beta_p2); //specificity prior on unvaxxed
  // target += normal_lpdf(epsilon| epsilon_p1, epsilon_p2) ; //specificity prior on vaxxed
  
  target += S1Q1 *log(gamma*(lambda*alpha+(1-lambda)*(1-epsilon))) ;
  target += S1Q0 *log((1-gamma)*(lambda*alpha+(1-lambda)*(1-beta)));
  target += S0Q1 *log(gamma*(lambda*(1-alpha)+(1-lambda)*epsilon)) ;
  target += S0Q0 *log((1-gamma)*(lambda*(1-alpha)+ (1-lambda)*beta));

}


