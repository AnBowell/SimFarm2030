data {
int<lower=0> n_regions;
int<lower=0> n_years;
real d_temp[n_regions,n_years,12];
real d_yields[n_regions,n_years];
real var_temp[n_regions,12];
}

parameters {
real s_temp[n_regions,12];
real<lower=0> width[n_regions];
}



transformed parameters {

real v_val[n_regions];
for (v in 1:n_regions) {
        v_val[v] = 0;
}


for (n in 1:n_regions){
for (j in 1:12){
v_val[n]=v_val[n]+square(s_temp[n,j])*var_temp[n,j];
}
}
}


model {
real tmp;
for (n in 1:n_regions){
width[n] ~uniform(0, 10);
for (m in 1:12){
s_temp[n,m] ~normal(0.0,100.0);
}
}
for (n in 1:n_regions){
for (y in 1:n_years){
tmp=0.0;
for (m in 1:12){
tmp=tmp+s_temp[n,m]*d_temp[n,y,m];
}
d_yields[n,y]~normal(tmp,width[n]);
}
}
}


generated quantities {
real d_yields_pred[n_regions,n_years];
real tmp;
for (n in 1:n_regions){
for (y in 1:n_years){
tmp=0.0;
for (m in 1:12){
tmp=tmp+s_temp[n,m]*d_temp[n,y,m];
}
d_yields_pred[n,y]=normal_rng(tmp,width[n]);
}
}
}