#ifndef __REC_FUNC
#define __REC_FUNC

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
