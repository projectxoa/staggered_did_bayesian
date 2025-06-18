 -- FILE hospital_allocation.stan --

data {
  intlower=1 T;              Number of time periods (e.g., weeks)
  intlower=1 P;              Number of professions
  intlower=1 D;              Number of disease categories

  matrix[T, P] Availability;     Availability[t, p] - e.g., hours or FTEs
  intlower=0 DiseaseCounts[T, D];  DiseaseCounts[t, d] - observed counts
}

parameters {
   Allocation_row_p[p] is a simplex of length D for profession p
   Each element Allocation_row_p[p][d] is the proportion of time
   profession p allocates to disease category d.
  simplex[D] Allocation_row_p[P];

   Rate_p Overall activity rate or throughput for profession p.
   This helps scale the allocated availability to the magnitude of disease counts.
  vectorlower=0[P] Rate_p;

   Intercept_d Baseline expected count for disease category d,
   independent of profession availability.
  vectorlower=0[D] Intercept_d;

   If you later want to use Negative Binomial for overdispersion
   reallower=0 phi;  Overdispersion parameter
}

model {
   --- Priors ---
   Priors for Allocation Proportions
   A Dirichlet(1,1,...,1) prior is flat, meaning all combinations of
   allocations summing to 1 are equally likely a priori for each profession.
  for (p_idx in 1P) {
    Allocation_row_p[p_idx] ~ dirichlet(rep_vector(1.0, D));
  }

   Priors for Rates (must be positive)
   HalfNormal(0, 5) means most mass is near 0 but allows for larger values.
   Adjust the scale '5' based on your expected magnitude of rates.
  Rate_p ~ half_normal(0, 5);

   Priors for Intercepts (must be non-negative for Poisson mean)
   Adjust the scale '10' based on typical baseline disease counts.
  Intercept_d ~ half_normal(0, 10);

   Prior for overdispersion (if using Negative Binomial)
   phi ~ exponential(0.1);  Example prior

   --- Likelihood ---
   Loop through each time period and each disease category
  for (t_idx in 1T) {
    for (d_idx in 1D) {
      real mu_td = Intercept_d[d_idx];  Start with the baseline intercept

       Add contributions from each profession
      for (p_idx in 1P) {
        mu_td += Availability[t_idx, p_idx]   How much of profession p is available
                   Allocation_row_p[p_idx][d_idx]   What proportion of p's time goes to disease d
                   Rate_p[p_idx];                   How effectively p's allocated time translates to counts
      }

       Ensure mu_td is positive (Poisson mean must be non-negative).
       With positive inputs and parameters, this should generally hold,
       but a small floor can help numerical stability if mu_td gets very close to zero.
      if (mu_td  1e-9) {
        mu_td = 1e-9;
      }

       Model observed disease counts
      DiseaseCounts[t_idx, d_idx] ~ poisson(mu_td);
       If you suspect overdispersion (variance of counts  mean of counts),
       you would use the Negative Binomial distribution instead
       DiseaseCounts[t_idx, d_idx] ~ neg_binomial_2(mu_td, phi);
    }
  }
}

generated quantities {
   This block is for calculating quantities of interest from the posterior samples.
   Here, we reformat the Allocation_row_p array into a more conventional PxD matrix
   for easier extraction and interpretation in RPython.
  matrix[P, D] Allocation_matrix;  Allocation_matrix[p,d]
  for (p_idx in 1P) {
    for (d_idx in 1D) {
      Allocation_matrix[p_idx, d_idx] = Allocation_row_p[p_idx][d_idx];
    }
  }
}