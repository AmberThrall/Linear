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

        static Complex<T> i() {
            return Complex<T>(T(0),T(1));
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
        friend Complex<T> operator+(const T& x, const Complex<T>& z) { return Complex<T>(x,0) += z; }
        friend Complex<T> operator-(Complex<T> z, const Complex<T>& w) { return z -= w; }
        friend Complex<T> operator-(Complex<T> z, const T& x) { return z -= x; }
        friend Complex<T> operator-(const T& x, const Complex<T>& z) { return Complex<T>(x,0) -= z; }
        friend Complex<T> operator*(Complex<T> z, const Complex<T>& w) { return z *= w; }
        friend Complex<T> operator*(Complex<T> z, const T& x) { return z *= x; }
        friend Complex<T> operator*(const T& x, Complex<T> z) { return z *= x; }
        friend Complex<T> operator/(Complex<T> z, const Complex<T>& w) { return z /= w; }
        friend Complex<T> operator/(Complex<T> z, const T& x) { return z /= x; }
        friend Complex<T> operator/(const T& x, const Complex<T>& z) { return Complex<T>(x,0) /= z; }
        // Unary operators.
        Complex<T> operator-() const {
            return Complex<T>(-this->Re, -this->Im);
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
            else if (z.Re == T(0)) {
                out << z.Im << "i";
            }
            else if (z.Im == T(0)) {
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
            return (z.Im == T(0) && z.Re == x);
        }
        friend bool operator==(const T& x, const Complex<T>& z) { return (z == x); }
        friend bool operator!=(const Complex<T>& z, const Complex<T>& w) { return !(z == w); }
        friend bool operator!=(const Complex<T>& z, const T& x) { return !(z == x); }
        friend bool operator!=(const T& x, const Complex<T>& z) { return !(z == x); }
    };

    using Complexi = Complex<int>;
    using Complexf = Complex<float>;
    using Complexd = Complex<double>;

    template <typename T>
    bool IsReal(T x) {
        return true;
    }
    template <typename T>
    bool IsReal(const Complex<T>& z) {
        return z.Im == T(0);
    }

    template <typename T>
    T Abs(T x) {
        if (x > T(0)) return x;
        return -x;
    }
    template <typename T>
    T Abs(const Complex<T>& z) {
        if (IsReal(z)) {
            return Abs(z.Re);
        }
        return std::sqrt(z.Re*z.Re+z.Im*z.Im);
    }

    template <typename T>
    T Sign(T x) {
        if (x > T(0)) return T(1);
        else if (x == T(0)) return T(0);
        else return T(-1);
    }
    template <typename T>
    Complex<T> Sign(const Complex<T>& z) {
        if (IsReal(z)) {
            return Sign(z.Re);
        }
        return z / Abs(z);
    }

    template <typename T>
    T Arg(T x) {
        return T(0);
    }
    template <typename T>
    T Arg(const Complex<T>& z) {
        return std::atan2(z.Im,z.Re);
    }

    template <typename T>
    T Conjugate(T x) {
        return x;
    }
    template <typename T>
    Complex<T> Conjugate(const Complex<T>& z) {
        return Complex<T>(z.Re, -z.Im);
    }

    template <typename T>
    T Sqrt(T x) {
        return std::sqrt(x);
    }
    template <typename T>
    Complex<T> Sqrt(const Complex<T>& z) {
        if (IsReal(z))
            return Sqrt(z.Re);
        T abs = Abs(z);
        T gamma = std::sqrt((z.Re+abs)/2);
        T delta = Sign(z.Im)*std::sqrt((-z.Re+abs)/2);
        return Complex<T>(gamma, delta);
    }

    template <typename T>
    T Exp(T x) {
        return std::exp(x);
    }
    template <typename T>
    Complex<T> Exp(const Complex<T>& z) {
        if (IsReal(z))
            return Exp(z.Re);
        T r = std::exp(z.Re);
        return Complex<T>(r*std::cos(z.Im), r*std::sin(z.Im));
    }

    template <typename T>
    T Log(T x) {
        return std::log(x);
    }
    template <typename T>
    Complex<T> Log(const Complex<T>& z) {
        if (IsReal(z))
            return Log(z.Re);
        return Complex<T>(std::log(Abs(z)), Arg(z));
    }
    template <typename T>
    T Log(T x, T base) {
        return std::log(x)/std::log(base);
    }
    template <typename T>
    Complex<T> Log(const Complex<T>& z, const Complex<T>& base) {
        if (IsReal(z) && IsReal(base))
            return Log(z.Re, base.Re);
        return Log(z)/Log(base);
    }
    template <typename T>
    Complex<T> Log(const Complex<T>& z, T base) {
        return Log(z, Complex<T>(base,0));
    }
    template <typename T>
    Complex<T> Log(T x, const Complex<T>& base) {
        return Log(Complex<T>(x,0), base);
    }

    template <typename T>
    T Pow(T x, T y) {
        return std::pow(x, y);
    }
    template <typename T>
    Complex<T> Pow(const Complex<T>& z, const Complex<T>& w) {
        if (IsReal(z) && IsReal(w))
            return Pow(z.Re, w.Re);
        else
            return Exp(w*Log(z));
    }

    template <typename T>
    T Mod(T x, T y) {
        return std::fmod(x, y);
    }
    template <typename T>
    Complex<T> Mod(const Complex<T>& z, T y) {
        return Complex<T>(Mod(z.Re, y), Mod(z.Im, y));
    }
    template <typename T>
    Complex<T> Mod(const Complex<T>& z, const Complex<T>& w) {
        return Complex<T>(Mod(z.Re, w.Re), Mod(z.Im, w.Im));
    }

    template <typename T>
    T Sin(T x) {
        return std::sin(x);
    }
    template <typename T>
    Complex<T> Sin(const Complex<T>& z) {
        if (IsReal(z))
            return Sin(z.Re);
        Complex<T> i = Complex<T>::i();
        return (Exp(i*z)-Exp(-i*z))/(2*i);
    }
    template <typename T>
    T Cos(T x) {
        return std::cos(x);
    }
    template <typename T>
    Complex<T> Cos(const Complex<T>& z) {
        if (IsReal(z))
            return Cos(z.Re);
        Complex<T> i = Complex<T>::i();
        return (Exp(i*z)+Exp(-i*z))/2;
    }
    template <typename T>
    T Tan(T x) {
        return std::tan(x);
    }
    template <typename T>
    Complex<T> Tan(const Complex<T>& z) {
        if (IsReal(z))
            return Tan(z.Re);
        Complex<T> i = Complex<T>::i();
        return ((Exp(2*i*z)-1)/(Exp(2*i*z)+1))/i;
    }
    template <typename T>
    T Csc(T x) {
        return T(1)/std::sin(x);
    }
    template <typename T>
    Complex<T> Csc(const Complex<T>& z) {
        if (IsReal(z))
            return Csc(z.Re);
        Complex<T> i = Complex<T>::i();
        return (2*i)/(Exp(i*z)-Exp(-i*z));
    }
    template <typename T>
    T Sec(T x) {
        return T(1)/std::cos(x);
    }
    template <typename T>
    Complex<T> Sec(const Complex<T>& z) {
        if (IsReal(z))
            return Sec(z.Re);
        Complex<T> i = Complex<T>::i();
        return 2/(Exp(i*z)+Exp(-i*z));
    }
    template <typename T>
    T Cot(T x) {
        return 1/std::tan(x);
    }
    template <typename T>
    Complex<T> Cot(const Complex<T>& z) {
        if (IsReal(z))
            return Cot(z.Re);
        Complex<T> i = Complex<T>::i();
        return i*((Exp(2*i*z)+1)/(Exp(2*i*z)-1));
    }

    template <typename T>
    T ASin(T x) {
        return std::asin(x);
    }
    template <typename T>
    Complex<T> ASin(const Complex<T>& z) {
        if (IsReal(z))
            return ASin(z.Re);
        Complex<T> i = Complex<T>::i();
        return Log(i*z+Sqrt(Abs(1-z*z))*Exp(Arg(1-z*z)*(i/2)))/i;
    }
    template <typename T>
    T ACos(T x) {
        return std::acos(x);
    }
    template <typename T>
    Complex<T> ACos(const Complex<T>& z) {
        if (IsReal(z))
            return ACos(z.Re);
        Complex<T> i = Complex<T>::i();
        return Log(z+i*Sqrt(Abs(1-z*z))*Exp(Arg(1-z*z)*i/2))/i;
    }
    template <typename T>
    T ATan(T x) {
        return std::atan(x);
    }
    template <typename T>
    Complex<T> ATan(const Complex<T>& z) {
        if (IsReal(z))
            return ATan(z.Re);
        Complex<T> i = Complex<T>::i();
        return Log((i-z)/(i+z))/(2*i);
    }
    template <typename T>
    T ACsc(T x) {
        return std::asin(T(1)/x);
    }
    template <typename T>
    Complex<T> ACsc(const Complex<T>& z) {
        if (IsReal(z))
            return ACsc(z.Re);
        return ASin(T(1)/z);
    }
    template <typename T>
    T ASec(T x) {
        return std::acos(T(1)/x);
    }
    template <typename T>
    Complex<T> ASec(const Complex<T>& z) {
        if (IsReal(z))
            return ASec(z.Re);
        return ACos(T(1)/z);
    }
    template <typename T>
    T ACot(T x) {
        return std::atan(T(1)/x);
    }
    template <typename T>
    Complex<T> ACot(const Complex<T>& z) {
        if (IsReal(z))
            return ACot(z.Re);
        Complex<T> i = Complex<T>::i();
        return Log((z+i)/(z-i))/(2*i);
    }

    template <typename T>
    T Sinh(T x) {
        return std::sinh(x);
    }
    template <typename T>
    Complex<T> Sinh(const Complex<T>& z) {
        if (IsReal(z))
            return Sinh(z.Re);
        return (Exp(z) - Exp(-z))/2;
    }
    template <typename T>
    T Cosh(T x) {
        return std::cosh(x);
    }
    template <typename T>
    Complex<T> Cosh(const Complex<T>& z) {
        if (IsReal(z))
            return Cosh(z.Re);
        return (Exp(z) + Exp(-z))/2;
    }
    template <typename T>
    T Tanh(T x) {
        return std::tanh(x);
    }
    template <typename T>
    Complex<T> Tanh(const Complex<T>& z) {
        if (IsReal(z))
            return Tanh(z.Re);
        return (Exp(z) - Exp(-z))/(Exp(z) + Exp(-z));
    }
    template <typename T>
    T Csch(T x) {
        return T(1)/std::sinh(x);
    }
    template <typename T>
    Complex<T> Csch(const Complex<T>& z) {
        if (IsReal(z))
            return Csch(z.Re);
        return 2/(Exp(z) - Exp(-z));
    }
    template <typename T>
    T Sech(T x) {
        return T(1)/std::cosh(x);
    }
    template <typename T>
    Complex<T> Sech(const Complex<T>& z) {
        if (IsReal(z))
            return Sech(z.Re);
        return 2/(Exp(z) + Exp(-z));
    }
    template <typename T>
    T Coth(T x) {
        return T(1)/std::tanh(x);
    }
    template <typename T>
    Complex<T> Coth(const Complex<T>& z) {
        if (IsReal(z))
            return Coth(z.Re);
        return (Exp(z) + Exp(-z))/(Exp(z) - Exp(-z));
    }

    template <typename T>
    T ASinh(T x) {
        return std::asinh(x);
    }
    template <typename T>
    Complex<T> ASinh(const Complex<T>& z) {
        if (IsReal(z))
            return ASinh(z.Re);
        Complex<T> i = Complex<T>::i();
        return Log(z+Sqrt(Abs(1+z*z))*Exp(Arg(1+z*z)*i/2));
    }
    template <typename T>
    T ACosh(T x) {
        return std::acosh(x);
    }
    template <typename T>
    Complex<T> ACosh(const Complex<T>& z) {
        if (IsReal(z))
            return ACosh(z.Re);
        Complex<T> i = Complex<T>::i();
        return Log(z+Sqrt(Abs(z*z-1))*Exp(Arg(z*z-1)*i/2));
    }
    template <typename T>
    T ATanh(T x) {
        return std::atanh(x);
    }
    template <typename T>
    Complex<T> ATanh(const Complex<T>& z) {
        if (IsReal(z))
            return ATanh(z.Re);
        return Log((1+z)/(1-z))/2;
    }
    template <typename T>
    T ACsch(T x) {
        return std::asinh(T(1)/x);
    }
    template <typename T>
    Complex<T> ACsch(const Complex<T>& z) {
        if (IsReal(z))
            return ACsch(z.Re);
        return ASinh(T(1)/z);
    }
    template <typename T>
    T ASech(T x) {
        return std::acos(T(1)/x);
    }
    template <typename T>
    Complex<T> ASech(const Complex<T>& z) {
        if (IsReal(z))
            return ASech(z.Re);
        return ACosh(T(1)/z);
    }
    template <typename T>
    T ACoth(T x) {
        return std::atanh(T(1)/x);
    }
    template <typename T>
    Complex<T> ACoth(const Complex<T>& z) {
        if (IsReal(z))
            return ACoth(z.Re);
        return Log((z+1)/(z-1))/2;
    }

    template <typename T>
    T Ceil(T x) {
        return std::ceil(x);
    }
    template <typename T>
    Complex<T> Ceil(const Complex<T>& z) {
        return Complex<T>(Ceil(z.Re), Ceil(z.Im));
    }
    template <typename T>
    T Floor(T x) {
        return std::floor(x);
    }
    template <typename T>
    Complex<T> Floor(const Complex<T>& z) {
        return Complex<T>(Floor(z.Re), Floor(z.Im));
    }
    template <typename T>
    T Round(T x) {
        return std::round(x);
    }
    template <typename T>
    Complex<T> Round(const Complex<T>& z) {
        return Complex<T>(Round(z.Re), Round(z.Im));
    }
}
