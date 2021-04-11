#pragma once
#include <ostream>
#include <cmath>

namespace Linear {
    template<typename T>
    class Complex {
    public:
        Complex(T re = 0, T im = 0) {
            this->Re = re;
            this->Im = im;
        }
        template <typename U>
        Complex(const Complex<U>& other) {
            this->Re = (T)other.Re;
            this->Im = (T)other.Im;
        }

        bool IsReal() const {
            return this->Im == T(0);
        }
        bool IsPure() const {
            return this->Re == T(0);
        }
        T Abs() const {
            return std::sqrt(this->Re*this->Re+this->Im*this->Im);
        }
        T Modulus() const {
            return Abs();
        }
        T Arg() const {
            return std::atan2(this->Im, this->Re);
        }
        Complex<T> Conjugate() const {
            return Complex(this->Re, -this->Im);
        }

        static Complex<T> i() {
            return Complex(0,1);
        }
        static Complex<T> Sqrt(const Complex<T>& z) {
            T abs = z.Abs();
            T gamma = std::sqrt((z.Re+abs)/2);
            T delta = (z.Im > T(0) ? T(1) : z.Im == T(0) ? T(0) : T(-1))*std::sqrt((-z.Re+abs)/2);
            return Complex(gamma, delta);
        }
        static Complex<T> Exp(const Complex<T>& z) {
            T r = std::exp(z.Re);
            return Complex(r*std::cos(z.Im), r*std::sin(z.Im));
        }
        static Complex<T> Log(const Complex<T>& z) {
            return Complex(std::log(z.Abs()), z.Arg());
        }
        static Complex<T> Pow(const Complex<T>& z, const Complex<T>& w) {
            if (z.IsReal() && w.IsReal())
                return Complex(std::pow(z.Re, w.Re), 0);
            else
                return Complex::Exp(w*Complex::Log(z));
        }

        T Re, Im;
        // Type conversion.
        template<typename U, typename std::enable_if<std::is_convertible<T,U>::value>::type* = nullptr>
        operator Complex<U>() {
            return Complex<U>((U)this->Re, (U)this->Im);
        }
        // Assignment operators.
        Complex & operator=(const Complex& other) {
            this->Re = other.Re;
            this->Im = other.Im;
            return *this;
        }
        Complex & operator=(const T& other) {
            this->Re = other;
            this->Im = T(0);
            return *this;
        }
        Complex & operator+=(const Complex& other) {
            this->Re += other.Re;
            this->Im += other.Im;
            return *this;
        }
        Complex & operator+=(const T& other) {
            this->Re += other;
            return *this;
        }
        Complex & operator-=(const Complex& other) {
            this->Re -= other.Re;
            this->Im -= other.Im;
            return *this;
        }
        Complex & operator-=(const T& other) {
            this->Re -= other;
            return *this;
        }
        Complex & operator*=(const Complex& other) {
            T re = this->Re*other.Re - this->Im*other.Im;
            T im = this->Re*other.Im + this->Im*other.Re;
            this->Re = re;
            this->Im = im;
            return *this;
        }
        Complex & operator*=(const T& other) {
            this->Re *= other;
            this->Im *= other;
            return *this;
        }
        Complex & operator/=(const Complex& other) {
            T abs = other.Re*other.Re+other.Im*other.Im;
            T re = (this->Re*other.Re+this->Im*other.Im)/abs;
            T im = (this->Im*other.Re-this->Re*other.Im)/abs;
            this->Re = re;
            this->Im = im;
            return *this;
        }
        Complex & operator/=(const T& other) {
            this->Re /= other;
            this->Im /= other;
            return *this;
        }
        // Binary operators.
        friend Complex<T> operator+(Complex<T> z, const Complex<T>& w) { return z += w; }
        friend Complex<T> operator+(Complex<T> z, const T& x) { return z += x; }
        friend Complex<T> operator+(const T& x, Complex<T> z) { return z += x; }
        friend Complex<T> operator-(Complex<T> z, const Complex<T>& w) { return z -= w; }
        friend Complex<T> operator-(Complex<T> z, const T& x) { return z -= x; }
        friend Complex<T> operator-(const T& x, Complex<T> z) { return z -= x; }
        friend Complex<T> operator*(Complex<T> z, const Complex<T>& w) { return z *= w; }
        friend Complex<T> operator*(Complex<T> z, const T& x) { return z *= x; }
        friend Complex<T> operator*(const T& x, Complex<T> z) { return z *= x; }
        friend Complex<T> operator/(Complex<T> z, const Complex<T>& w) { return z /= w; }
        friend Complex<T> operator/(Complex<T> z, const T& x) { return z /= x; }
        friend Complex<T> operator/(const T& x, const Complex<T>& z) { return Complex<T>(x,0) /= z; }
        // Unary operators.
        Complex<T> operator-() {
            return Complex(-this->Re, -this->Im);
        }
        Complex<T> & operator++() { //prefix
            ++this->Re;
            return *this;
        }
        Complex<T> & operator++(int) { //postfix
            this->Re++;
            return *this;
        }
        Complex<T> & operator--() { //prefix
            --this->Re;
            return *this;
        }
        Complex<T> & operator--(int) { //postfix
            this->Re--;
            return *this;
        }
        // Output operators.
        friend std::ostream& operator<<(std::ostream &out, const Complex<T>& z) {
            if (z.Re == T(0) && z.Im == T(0)) {
                out << "0";
            }
            else if (z.IsPure()) {
                out << z.Im << "i";
            }
            else if (z.IsReal()) {
                out << z.Re;
            }
            else {
                out << z.Re << (z.Im < T(0) ? "-" : "+") << std::abs(z.Im) << "i";
            }
            return out;
        }

        // Comparison operators.
        friend bool operator==(const Complex<T>& z, const Complex<T>& w) {
            return (z.Re == w.Re && z.Im == w.Im);
        }
        friend bool operator==(const Complex<T>& z, const T& x) {
            return (z.IsReal() && z.Re == x);
        }
        friend bool operator==(const T& x, const Complex<T>& z) { return (z == x); }
        friend bool operator!=(const Complex<T>& z, const Complex<T>& w) { return !(z == w); }
        friend bool operator!=(const Complex<T>& z, const T& x) { return !(z == x); }
        friend bool operator!=(const T& x, const Complex<T>& z) { return !(z == x); }
    };

    using Complexi = Complex<int>;
    using Complexf = Complex<float>;
    using Complexd = Complex<double>;
}
