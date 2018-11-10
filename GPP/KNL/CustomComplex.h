/*
Templated CustomComplex class that represents a complex class comprised of  any type of real and imaginary types.
*/
#ifndef __CustomComplex
#define __CustomComplex

#include <iostream>
#include <cstdlib>
#include <memory>
#include <iomanip>
#include <cmath>
#include <omp.h>
#include <ctime>
#include <stdio.h>
#include <sys/time.h>

template<class re, class im>

class CustomComplex {

    private : 
    re x;
    im y;

    public:
    explicit CustomComplex () {
        x = 0.00;
        y = 0.00;
    }


    explicit CustomComplex(const double& a, const double& b) {
        x = a;
        y = b;
    }

    CustomComplex(const CustomComplex& src) {
        x = src.x;
        y = src.y;
    }

    CustomComplex& operator =(const CustomComplex& src) {
        x = src.x;
        y = src.y;

        return *this;
    }

    CustomComplex& operator +=(const CustomComplex& src) {
        x = src.x + this->x;
        y = src.y + this->y;

        return *this;
    }

    CustomComplex& operator -=(const CustomComplex& src) {
        x = src.x - this->x;
        y = src.y - this->y;

        return *this;
    }

    CustomComplex& operator -() {
        x = -this->x;
        y = -this->y;

        return *this;
    }

    CustomComplex& operator ~() {
        return *this;
    }

    void print() const {
        printf("( %f, %f) ", this->x, this->y);
        printf("\n");
    }

    double get_real() const
    {
        return this->x;
    }

    double get_imag() const
    {
        return this->y;
    }

    void set_real(double val)
    {
        this->x = val;
    }

    void set_imag(double val) 
    {
        this->y = val;
    }

// 6 flops
    template<class real, class imag>
    friend inline CustomComplex<real,imag> operator *(const CustomComplex<real,imag> &a, const CustomComplex<real,imag> &b) {
        real x_this = a.x * b.x - a.y*b.y ;
        imag y_this = a.x * b.y + a.y*b.x ;
        CustomComplex<real,imag> result(x_this, y_this);
        return (result);
    }

//2 flops
    template<class real, class imag>
    friend inline CustomComplex<real,imag> operator *(const CustomComplex<real,imag> &a, const double &b) {
       CustomComplex<real,imag> result(a.x*b, a.y*b);
       return result;
    }

//2 flops
    template<class real, class imag>
    friend inline CustomComplex<real,imag> operator -(const double &a, CustomComplex<real,imag>& src) {
        CustomComplex<real,imag> result(a - src.x, 0 - src.y);
        return result;
    }

    template<class real, class imag>
    friend inline CustomComplex<real,imag> operator +(const double &a, CustomComplex<real,imag>& src) {
        CustomComplex<real,imag> result(a + src.x, src.y);
        return result;
    }

    template<class real, class imag>
    friend inline CustomComplex<real,imag> CustomComplex_conj(const CustomComplex<real,imag>& src) ;

    template<class real, class imag>
    friend inline double CustomComplex_abs(const CustomComplex<real,imag>& src) ;

    template<class real, class imag>
    friend inline double CustomComplex_real( const CustomComplex<real,imag>& src) ;

    template<class real, class imag>
    friend inline double CustomComplex_imag( const CustomComplex<real,imag>& src) ;
};

/*
 * Return the conjugate of a complex number 
 1flop
 */
template<class re, class im>
inline CustomComplex<re, im> CustomComplex_conj(const CustomComplex<re,im>& src) {

    re re_this = src.x;
    im im_this = -1 * src.y;

    CustomComplex<re,im> result(re_this, im_this);
    return result;

}

/*
 * Return the absolute of a complex number 
 */
template<class re, class im>
inline double CustomComplex_abs(const CustomComplex<re,im>& src) {
    re re_this = src.x * src.x;
    im im_this = src.y * src.y;

    re result = sqrt(re_this+im_this);
    return result;
}

/*
 * Return the real part of a complex number 
 */
template<class re, class im>
inline double CustomComplex_real( const CustomComplex<re,im>& src) {
    return src.x;
}

/*
 * Return the imaginary part of a complex number 
 */
template<class re, class im>
inline double CustomComplex_imag( const CustomComplex<re,im>& src) {
    return src.y;
}

#endif
