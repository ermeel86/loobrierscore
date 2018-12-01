functions {

    vector row_sums(matrix X) {
        vector[rows(X)] s ;  
        s = X * rep_vector(1, cols(X)) ;
        return s ;
    }

    // gammas has nbasis rows and nsamples columns
    // betas has ncovariates rows and nsamples columns!!!
    matrix surv_mspline_t(real[] i_spline_basis_evals, matrix X, matrix gammas, real[] intercepts, matrix betas) {
        int npatients = rows(X);
        return exp(rep_matrix(-to_row_vector(i_spline_basis_evals) * gammas, npatients).*exp(X*betas + rep_matrix(to_row_vector(intercepts), npatients)));
    }

    // X has npatient rows and ncovariates columns
    // betas has ncovariates rows and nsamples columns!!!
    matrix surv_const_t(real time, matrix X, real[] intercepts, matrix betas) {
        int npatients = rows(X);
        return exp(-time*exp(X*betas + rep_matrix(to_row_vector(intercepts), npatients)));
    }
    // survs has npatients rows and nsamples columns
    // psis_ws has npatients rows and nsamples columns
    vector surv_loo(matrix survs, matrix psis_ws) {
        return row_sums(survs .* psis_ws);
    }
    /***************************************************************************************/
    /***************************************************************************************/
    matrix surv_mspline(matrix i_spline_basis_evals, matrix X, real[] gammas, real intercept,real[] betas) {
        vector[rows(i_spline_basis_evals)] basehaz;
        vector[rows(X)] hrs;
        matrix[rows(X), num_elements(basehaz)] survs;
        hrs = exp(intercept + X * to_vector(betas));
        basehaz = -i_spline_basis_evals * to_vector(gammas);
        survs = rep_matrix(to_row_vector(basehaz), rows(X));
        for(i in 1:rows(survs)) {
            survs[i,] = exp(survs[i,]*hrs[i]);
        }
        return survs;
        
  }

    matrix surv_const(real[] times, matrix X,  real intercept,real[] betas) {
        vector[size(times)] basehaz;
        vector[rows(X)] hrs;
        matrix[rows(X), num_elements(basehaz)] survs;
        hrs = exp(X * to_vector(betas));
        basehaz = -to_vector(times) * exp(intercept);
        survs = rep_matrix(to_row_vector(basehaz), rows(X));
        for(i in 1:rows(survs)) {
            survs[i,] = exp(survs[i,]*hrs[i]);
        }
        return survs;
        
  }

    // psis_ws has npatients rows and nsamps columns!!! 
    // X has npatients rows and ncovs columns
    // betas has nsamps rows and npatients columns
    matrix surv_const_loo(real [] times, matrix X, real[] intercepts, matrix betas, matrix psis_ws) {
        int ntimes = size(times);
        int npatients = rows(X);
        int nsamps = size(intercepts);
        matrix[npatients, nsamps] survs_tmp;
        matrix[npatients, ntimes] survs_loo;
        matrix[npatients,1] sm;
        real ts[1];
        for(t in 1:ntimes) {
            ts[1] = times[t];
            // todo: maybe transform survs_tmp for memory locality
            // i think this is fine, since matrix are column major in Stan!
            // see manual (version 2.17) chapter 26.3
            for(i in 1:nsamps) {
                sm = surv_const(ts, X, intercepts[i], to_array_1d(betas[i,]));
                survs_tmp[,i] = sm[,1];
            }
            survs_loo[,t] = row_sums(survs_tmp .* psis_ws);
        }
        return survs_loo; 
    }




    // i_spline_basis_evals has ntimes rows and nbasis columns
    // gammas has nsamps rows and nbasis columns
    // psis_ws has npatients rows and nsamps columns!!! 
    // X has npatients rows and ncovs columns
    // betas has nsamps rows and npatients columns
    matrix surv_mspline_loo(matrix i_spline_basis_evals, matrix X, matrix gammas, 
    real[] intercepts, matrix betas, matrix psis_ws) {
        int ntimes = rows(i_spline_basis_evals);
        int npatients = rows(X);
        int nbasis = cols(i_spline_basis_evals);
        int nsamps = size(intercepts);
        matrix[npatients, nsamps] survs_tmp;
        matrix[npatients, ntimes] survs_loo;
        matrix[npatients,1] sm;
        matrix[1,nbasis] isp;
        for(t in 1:ntimes) {
            isp[1,] = i_spline_basis_evals[t,];
            // todo: maybe transform survs_tmp for memory locality
            // i think this is fine, since matrix are column major in Stan!
            // see manual (version 2.17) chapter 26.3
            for(i in 1:nsamps) {
                sm = surv_mspline(isp, X, to_array_1d(gammas[i,]),intercepts[i],  to_array_1d(betas[i,]));
                survs_tmp[,i] = sm[,1];
            }
            survs_loo[,t] = row_sums(survs_tmp .* psis_ws);
        }
        return survs_loo; 
    }
  
}
