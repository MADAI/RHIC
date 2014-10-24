/// Efficient templates for calculating integer powers.
#pragma once

/// Calculate integer powers using squaring.
template<class T>
inline constexpr T pow(const T base, unsigned const exponent) {
    return (exponent == 0)     ? 1 :
           (exponent % 2 == 0) ? pow(base, exponent/2) * pow(base, exponent/2) :
                                 base * pow(base, exponent - 1);
}

template<class T>
inline constexpr T square(const T base) {
    return pow(base, 2);
}
