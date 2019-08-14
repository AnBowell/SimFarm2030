functions {
    real d_yield(real d_t, 
                 real d_p, 
                 real mu_t, 
                 real mu_p, 
                 real sigma_t, 
                 real sigma_p, 
                 real rho){
        // return d_t; //test function can run 
        return (1/(2 * pi() * sigma_t * sigma_p * sqrt(1 - square(rho)))) * 
            exp( 
                - (1/(2* (1- square(rho) ))) * 
                (
                    square((d_t - mu_t)/(sigma_t))
                    + square((d_p - mu_p)/(sigma_p))
                    - (2* rho * (d_t - mu_t) * (d_p - mu_p) / (sigma_p * sigma_p) ) 
                )
            );
 
    }
}

data {
    int<lower=0> n_regions;
    int<lower=0> n_years;
    real d_temp[n_regions,n_years,12];
    real d_precip[n_regions,n_years,12];
    real d_yields[n_regions,n_years];
}

parameters {
    real mu_t[n_regions,12];
    real mu_p[n_regions,12];
    //Wikipedia notation
    real y_norm[n_regions,12];
    real<lower=0.0> sigma_t[n_regions,12];
    real<lower=0.0> sigma_p[n_regions,12];
    real<lower=-1.0, upper=+1.0> rho[n_regions,12];
    //real noise_sigma[n_regions];

}


model {
    real tmp;
    real d_y;
    for (n in 1:n_regions){
        for (m in 1:12){
            mu_t[n_regions,12] ~normal(20.0,100.0);
            mu_p[n_regions,12] ~normal(500.0,1000.0);
            y_norm[n_regions,12] ~normal(6.0,10.0);
            sigma_t[n_regions,12] ~normal(0.0,100.0);
            sigma_p[n_regions,12] ~normal(0.0,1000.0);
            rho[n_regions,12] ~normal(0.0,1000.0);
            //noise_sigma[n_regions] ~normal(0.0,100.0);
        }
    }
    for (n in 1:n_regions){
        
        for (y in 1:n_years){
            d_y = 0.0;
            tmp=0.0;
            for (m in 1:12){
                d_y = y_norm[n,m] * d_yield(d_temp[n,y,m], 
                                            d_precip[n,y,m], 
                                            mu_t[n,m], 
                                            mu_p[n,m], 
                                            sigma_t[n,m], 
                                            sigma_p[n,m], 
                                            rho[n,m]);
                tmp=tmp+d_y;
            }
            d_yields[n,y]~normal(tmp,1.0);
        }
    }
}


//generated quantities {
//    real d_yields_pred[n_regions,n_years];
//    real tmp;
//    for (n in 1:n_regions){
//        for (y in 1:n_years){
//            tmp=0.0;
//            for (m in 1:12){
//                tmp=tmp+s_temp[n,m]*d_temp[n,y,m] + s_precip[n,m]*d_precip[n,y,m];
//            }
//            d_yields_pred[n,y]=normal_rng(tmp,noise_sigma);
//        }
//    }
//}