data{
    int<lower=0> NSUB;                      // number of patients
    int<lower=0> NOBS[NSUB];                // number of observations for each patient
    int<lower=0> NDOSE[NSUB];               // number of doses for each patient
    vector[sum(NOBS)] conc;                 // observations of concentration
    vector<lower=0>[sum(NOBS)] obs_time;    // observation time for all patients
    vector<lower=0>[sum(NDOSE)] dose_time;  // dose time for all patients
    vector<lower=0>[sum(NDOSE)] dose_amt;   // dose amounts for all patients 
}
 
parameters{
    vector<lower=-5, upper=5>[6] theta;       //fixed effect
    vector[6] eta[NSUB];              //inter-individual random effect
    vector<lower=0, upper=2>[6] sigma_eta;     //variance of inter-individual random effect
    real<lower=0> sigma;              //variance of intra-individual random effect
}
 
transformed parameters{
    vector[NSUB] CL;
    vector[NSUB] V;
    vector[NSUB] Q;
    vector[NSUB] V2;
    vector[NSUB] ka;
    vector[NSUB] tlag;

    vector[sum(NOBS)] y_pred;

    for(n in 1:NSUB){
        CL[n] <- exp( theta[1] + sigma_eta[1] * eta[n,1] );
        V[n] <- exp( theta[2] + sigma_eta[2] * eta[n,2] );
        Q[n] <- exp( theta[3] + sigma_eta[3] * eta[n,3] );
        V2[n] <- exp( theta[4] + sigma_eta[4] * eta[n,4] );
        ka[n] <- exp( theta[5] + sigma_eta[5] * eta[n,5] );
        tlag[n] <- exp( theta[6] + sigma_eta[6] * eta[n,6] );

    }
    
    {
        int y_index;
        int dose_index;
        y_index <- 1;
        dose_index <- 1;

        for(i in 1:NSUB){
            vector[NOBS[i]] g;
            vector[6] params;           
            params[1] <- CL[i];
            params[2] <- V[i];
            params[3] <- Q[i];
            params[4] <- V2[i];
            params[5] <- ka[i];
            params[6] <- tlag[i];

            g <- linear_cmpt_1order_absor(
                     segment(obs_time, y_index, NOBS[i]),
                     segment(dose_time, dose_index, NDOSE[i]),
                     segment(dose_amt, dose_index, NDOSE[i]), 
                     params, 
                     2,    // number of compartment(s) 
                     1);   // parameterization option: 1(CL_V), 2(micro_rate)

            for(j in 1:NOBS[i])
                y_pred[y_index + j - 1] <- g[j];
            y_index <- y_index + NOBS[i];
            dose_index <- dose_index + NDOSE[i];

        }// end of for loop
    } //end of local variable
}//end of transformed parameters block
 
model{
    for(k in 1:6){
        for(i in 1:NSUB)
            eta[i,k] ~ normal(0.,1.);
        theta[k] ~ normal(0.,1000.);
        sigma_eta[k] ~ normal(0.,1000.);
    }
    sigma ~ normal(0.,1000.);
    conc ~ normal(y_pred, sigma);
}

generated quantities{
    vector[sum(NOBS)] log_lik;
    for(n in 1:sum(NOBS)){
        log_lik[n] <- normal_log(conc[n], y_pred[n], sigma);
    }
}
