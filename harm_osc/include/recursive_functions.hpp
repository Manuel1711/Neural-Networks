#ifndef __REC_FUNC
#define __REC_FUNC


#ifdef ACT_ERF

inline complex mathT(const complex a, const complex b, const complex c, const complex d){
    return 2./datum::pi * asin( 2.* b / sqrt((1.+2.*a)*(1.+2.*d)));
}

inline complex mathTp(const complex a, const complex b, const complex c, const complex d){
    /*
    auxMat = np.empty((2,2))
    auxMat[0,0] = a
    auxMat[0,1] = b
    auxMat[1,0] = c
    auxMat[1,1] = d
    return 4./ np.pi * det(fractional_matrix_power(eye(2) + 2*auxMat, -1./2.))
    */
    cx_dmat auxMat(2, 2);

    // Fill auxMat with the given values
    auxMat(0, 0) = a;
    auxMat(0, 1) = b;
    auxMat(1, 0) = c;
    auxMat(1, 1) = d;

    // Calculate the result using Armadillo functions
    cx_dmat identityMat = eye<cx_dmat>(2, 2);
    cx_dmat sumMat = identityMat + 2 * auxMat;

    complex detResult = 4.0 / datum::pi * det(sqrtmat(inv_sympd(sumMat)));

    return detResult;
}

#endif

#ifdef ACT_RELU

inline complex mathT(const complex a, const complex b, const complex c, const complex d) {
    complex a1 = abs(a);
    complex d1 = abs(d);
    complex fact_or = b / sqrt(a1 * d1);
    complex output = sin(acos(fact_or)) + (datum::pi - acos(fact_or)) * fact_or;
    return 1.0 / 2.0 / datum::pi * sqrt(a1 * d1) * output;
}

inline complex mathTp(const complex a, const complex b, const complex c, const complex d) {
    complex a1 = abs(a);
    complex d1 = abs(d);
    return 1.0 / 2.0 / datum::pi * (datum::pi - acos(b / sqrt(a1 * d1)));
}

#endif

// Function to calculate Kappa
complex Kappa(const cx_dvec& arr1, const cx_dvec& arr2, const int l_index, const double sigma_w, const double sigma_b, const double n0) {
    if (l_index == 1) {
        return dot(arr1, arr2) / n0 * sigma_w * sigma_w + sigma_b * sigma_b;
    } else {
        complex T11 = Kappa(arr1, arr1, l_index - 1, sigma_w, sigma_b, n0);
        complex T12 = Kappa(arr1, arr2, l_index - 1, sigma_w, sigma_b, n0);
        complex T21 = T12;
        complex T22 = Kappa(arr2, arr2, l_index - 1, sigma_w, sigma_b, n0);
        return sigma_w * sigma_w * mathT(T11, T12, T21, T22) + sigma_b * sigma_b;
    }
}

// Function to calculate Theta
complex Theta(const cx_dvec& arr1, const cx_dvec& arr2, const int l_index, const double sigma_w, const double sigma_b, const double n0) {
    if (l_index == 1) {
        return dot(arr1, arr2) / n0 * sigma_w * sigma_w + sigma_b * sigma_b;
    } else {
        complex T11 = Kappa(arr1, arr1, l_index - 1, sigma_w, sigma_b, n0);
        complex T12 = Kappa(arr1, arr2, l_index - 1, sigma_w, sigma_b, n0);
        complex T21 = T12;
        complex T22 = Kappa(arr2, arr2, l_index - 1, sigma_w, sigma_b, n0);
        return sigma_w * sigma_w * mathTp(T11, T12, T21, T22) * Theta(arr1, arr2, l_index - 1, sigma_w, sigma_b, n0);
    }
}

#endif
