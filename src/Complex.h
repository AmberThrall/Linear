#pragma once
#include <ostream>
#include <cmath>

namespace Linear {
    /**
     * Class for complex numbers.
     * @param T Type to store the real and imaginary part as.
     */
    template<typename T>
    class Complex {
    public:
        /**
         * Constructor.
         * Creates the complex number re+im*i.
         * @param re The real part of the complex number. (default = 0)
         * @param im The imaginary part of the complex number. (default = 0)
         */
        Complex(T re = 0, T im = 0) {
            this->Re = re;
            this->Im = im;
        }
        /**
         * Constructor.
         * Copies a complex number.
         * @param other Complex number to copy.
         */
        template <typename U>
        Complex(const Complex<U>& other) {
            this->Re = (T)other.Re;
            this->Im = (T)other.Im;
        }

        /**
         * Static function.
         * Creates and returns the complex number i.
         * @return Complex
         */
        static Complex<T> i() {
            return Complex<T>(T(0),T(1));
        }

        T Re; /*!< Real part of the complex number */
        T Im; /*!< Imaginary part of the complex number */
        // Type conversion.
        template<typename U, typename std::enable_if<std::is_convertible<T,U>::value>::type* = nullptr>
        operator Complex<U>() {
            return Complex<U>((U)this->Re, (U)this->Im);
        }
        // Assignment operators.
        /**
         * Assigns the complex number to another complex number.
         * @param other Complex number.
         * @return Reference to complex number.
         */
        Complex & operator=(const Complex& other) {
            this->Re = other.Re;
            this->Im = other.Im;
            return *this;
        }
        /**
         * Assigns the complex number to the real number other.
         * @param other Real number
         * @return Reference to complex number.
         */
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
        /**
         * Takes two complex numbers \f$z=a+bi\f$ and \f$w=c+di\f$ and returns the complex number \f$z+w=(a+c)+(b+d)i\f$.
         * @param z Complex number.
         * @param w Complex number.
         * @return Complex number.
         */
        friend Complex<T> operator+(Complex<T> z, const Complex<T>& w) { return z += w; }
        /**
         * Takes a complex number \f$z=a+bi\f$ and real number \f$x\f$ and returns the complex number \f$z+x=(a+x)+bi\f$.
         * @param z Complex number.
         * @param x Real number.
         * @return Complex number.
         */
        friend Complex<T> operator+(Complex<T> z, const T& x) { return z += x; }
        /**
         * Takes a real number \f$x\f$ and complex number \f$z=a+bi\f$ returns the complex number \f$x+z=(x+a)+bi\f$.
         * @param x Real number.
         * @param z Complex number.
         * @return Complex number.
         */
        friend Complex<T> operator+(const T& x, const Complex<T>& z) { return Complex<T>(x,0) += z; }
        /**
         * Takes two complex numbers \f$z=a+bi\f$ and \f$w=c+di\f$ and returns the complex number \f$z-w=(a-c)+(b-d)i\f$.
         * @param z Complex number.
         * @param w Complex number.
         * @return Complex number.
         */
        friend Complex<T> operator-(Complex<T> z, const Complex<T>& w) { return z -= w; }
        /**
         * Takes a complex number \f$z=a+bi\f$ and real number \f$x\f$ and returns the complex number \f$z-x=(a-x)+bi\f$.
         * @param z Complex number.
         * @param x Real number.
         * @return Complex number.
         */
        friend Complex<T> operator-(Complex<T> z, const T& x) { return z -= x; }
        /**
         * Takes a real number \f$x\f$ and complex number \f$z=a+bi\f$ returns the complex number \f$x-z=(x-a)-bi\f$.
         * @param x Real number.
         * @param z Complex number.
         * @return Complex number.
         */
        friend Complex<T> operator-(const T& x, const Complex<T>& z) { return Complex<T>(x,0) -= z; }
        /**
         * Takes two complex numbers \f$z=a+bi\f$ and \f$w=c+di\f$ and returns the complex number \f$zw=(ac-bd)+(ad+bc)i\f$.
         * @param z Complex number.
         * @param w Complex number.
         * @return Complex number.
         */
        friend Complex<T> operator*(Complex<T> z, const Complex<T>& w) { return z *= w; }
        /**
         * Takes a complex number \f$z=a+bi\f$ and real number \f$x\f$ and returns the complex number \f$zx=ax+bxi\f$.
         * @param z Complex number.
         * @param x Real number.
         * @return Complex number.
         */
        friend Complex<T> operator*(Complex<T> z, const T& x) { return z *= x; }
        /**
         * Takes a real number \f$x\f$ and complex number \f$z=a+bi\f$ returns the complex number \f$xz=xa+xbi\f$.
         * @param x Real number.
         * @param z Complex number.
         * @return Complex number.
         */
        friend Complex<T> operator*(const T& x, Complex<T> z) { return z *= x; }
        /**
         * Takes two complex numbers \f$z=a+bi\f$ and \f$w=c+di\f$ and returns the complex number \f$\frac{z}{w}=\frac{1}{c^2+d^2}((ac+bd)+(bc-ad)i)\f$.
         * @param z Complex number.
         * @param w Complex number.
         * @return Complex number.
         */
        friend Complex<T> operator/(Complex<T> z, const Complex<T>& w) { return z /= w; }
        /**
         * Takes a complex number \f$z=a+bi\f$ and real number \f$x\f$ and returns the complex number \f$\frac{z}{x}=\frac{a}{x}+\frac{b}{x}i\f$.
         * @param z Complex number.
         * @param x Real number.
         * @return Complex number.
         */
        friend Complex<T> operator/(Complex<T> z, const T& x) { return z /= x; }
        /**
         * Takes a real number \f$x\f$ and complex number \f$z=a+bi\f$ returns the complex number \f$\frac{x}{z}=\frac{1}{a^2+b^2}(xa-xbi)\f$.
         * @param x Real number.
         * @param z Complex number.
         * @return Complex number.
         */
        friend Complex<T> operator/(const T& x, const Complex<T>& z) { return Complex<T>(x,0) /= z; }
        // Unary operators.
        Complex<T> operator-() const {
            return Complex<T>(-this->Re, -this->Im);
        }
        // Output operators.
        /**
         * Takes a complex number \f$z=a+bi\f$ and writes to out. If \f$z=0\f$, then simply \f$0\f$ will be written.
         * If \f$z=bi\f$, then it will write just \f$bi\f$. If \f$z=a\f$, then it will write just \f$a\f$. Otherwise the full
         * \f$a+bi\f$ will be written, with "+" replaced by "-" if \f$b<0\f$.
         * @param out Reference to std::ostream.
         * @param z Complex number.
         */
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
        /**
         * Compares two complex numbers \f$z=a+bi\f$ and \f$w=c+di\f$.
         * @param z Complex number.
         * @param w Complex number.
         * @return If \f$a=c\f$ and \f$b=d\f$ it returns true. Otherwise it returns false.
         */
        friend bool operator==(const Complex<T>& z, const Complex<T>& w) {
            return (z.Re == w.Re && z.Im == w.Im);
        }
        /**
         * Compares a complex number \f$z=a+bi\f$ and a real number \f$x\f$.
         * @param z Complex number.
         * @param x Real number.
         * @return If \f$a=x\f$ and \f$b=0\f$ it returns true. Otherwise it returns false.
         */
        friend bool operator==(const Complex<T>& z, const T& x) {
            return (z.Im == T(0) && z.Re == x);
        }
        /**
         * Compares a real number \f$x\f$ and a complex numbers \f$z=a+bi\f$.
         * @param x Real number.
         * @param z Complex number.
         * @return If \f$a=x\f$ and \f$b=0\f$ it returns true. Otherwise it returns false.
         */
        friend bool operator==(const T& x, const Complex<T>& z) { return (z == x); }
        /**
         * Compares two complex numbers \f$z=a+bi\f$ and \f$w=c+di\f$.
         * @param z Complex number.
         * @param w Complex number.
         * @return If \f$a\ne c\f$ or \f$b\ne d\f$ it returns true. Otherwise it returns false.
         */
        friend bool operator!=(const Complex<T>& z, const Complex<T>& w) { return !(z == w); }
        /**
         * Compares a complex number \f$z=a+bi\f$ and a real number \f$x\f$.
         * @param z Complex number.
         * @param x Real number.
         * @return If \f$a\ne x\f$ or \f$b\ne 0\f$ it returns true. Otherwise it returns false.
         */
        friend bool operator!=(const Complex<T>& z, const T& x) { return !(z == x); }
        /**
         * Compares a real number \f$x\f$ and a complex numbers \f$z=a+bi\f$.
         * @param x Real number.
         * @param z Complex number.
         * @return If \f$a\ne x\f$ or \f$b\ne 0\f$ it returns true. Otherwise it returns false.
         */
        friend bool operator!=(const T& x, const Complex<T>& z) { return !(z == x); }
    };

    using Complexi = Complex<int>;
    using Complexf = Complex<float>;
    using Complexd = Complex<double>;

    /**
     * Checks if a complex number \f$z=a+bi\f$ is real.
     * @param z Complex number.
     * @return Returns true if \f$b=0\f$.
     */
    template <typename T>
    bool IsReal(const Complex<T>& z) {
        return z.Im == T(0);
    }

    /**
     * Takes a real number \f$x\f$ and computes \f$|x|\f$.
     * @param x Real number.
     * @return Real number.
     */
    template <typename T>
    T Abs(T x) {
        if (x > T(0)) return x;
        return -x;
    }
    /**
     * Takes a complex number \f$z=a+bi\f$ and computes \f$|z|=\sqrt{a^2+b^2}\f$.
     * @param z Complex number.
     * @return Real number.
     */
    template <typename T>
    T Abs(const Complex<T>& z) {
        if (IsReal(z)) {
            return Abs(z.Re);
        }
        return std::sqrt(z.Re*z.Re+z.Im*z.Im);
    }

    /**
     * Takes a real number \f$x\f$ and computes \f$sgn(x)\f$.
     * If \f$x>0\f$ it returns \f$1\f$. If \f$x=0\f$ it returns \f$0\f$. Otherwise it returns \f$-1\f$.
     * @param x Real number.
     * @return -1, 0 or 1.
     */
    template <typename T>
    T Sign(T x) {
        if (x > T(0)) return T(1);
        else if (x == T(0)) return T(0);
        else return T(-1);
    }
    /**
     * Takes a complex number \f$z\f$ and computes \f$sgn(z)\f$.
     * If \f$z\f$ is real and if \f$z>0\f$ it returns \f$1\f$. If \f$z\f$ is real and \f$z=0\f$ it returns \f$0\f$.
     * Otherwise if \f$z\f$ is real it returns \f$-1\f$. If \f$z\f$ is not real, it returns \f$\frac{z}{|z|}\f$.
     * @param z Complex number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> Sign(const Complex<T>& z) {
        if (IsReal(z)) {
            return Sign(z.Re);
        }
        return z / Abs(z);
    }

    /**
     * Takes a complex number \f$z=a+bi\f$ and computes \f$arg(z)=atan2(b,a)\f$.
     * @param z Complex number.
     * @return Real number.
     */
    template <typename T>
    T Arg(const Complex<T>& z) {
        return std::atan2(z.Im,z.Re);
    }

    /**
     * Takes a complex number \f$z=a+bi\f$ and computes it's conjugate \f$\bar{z}=a-bi\f$.
     * @param z Complex number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> Conjugate(const Complex<T>& z) {
        return Complex<T>(z.Re, -z.Im);
    }

    /**
     * Computes the square root of a real number \f$x\f$. If \f$x<0\f$ it returns the complex number \f$\sqrt{|x|}i\f$.
     * @param x Real number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> Sqrt(T x) {
        if (x < T(0))
            return Complex<T>(T(0), std::sqrt(Abs(x)));
        return Complex<T>(std::sqrt(x), T(0));
    }
    /**
     * Takes a complex number \f$z=a+bi\f$ and computes it's principal square root \f$\sqrt{z}=\gamma+\delta i\f$,
     * where \f$\gamma=\sqrt{\frac{a+|z|}{2}}\f$ and \f$\delta=(sgn b)\sqrt{\frac{-a+|z|}{2}}\f$.
     * @param z Complex number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> Sqrt(const Complex<T>& z) {
        if (IsReal(z))
            return Sqrt(z.Re);
        T abs = Abs(z);
        T gamma = std::sqrt((z.Re+abs)/2);
        T delta = Sign(z.Im)*std::sqrt((-z.Re+abs)/2);
        return Complex<T>(gamma, delta);
    }

    /**
     * Takes a real number \f$x\f$ and computes \f$e^x\f$.
     * @param x Real number.
     * @return Real number.
     */
    template <typename T>
    T Exp(T x) {
        return std::exp(x);
    }
    /**
     * Takes a complex number \f$z\f$ and computes \f$e^z\f$.
     * This is done using Euler's formula. Thus \f$e^z=e^{a+bi}=e^ae^{bi}=e^a(\cos b+i\sin b)\f$.
     * @param z Complex number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> Exp(const Complex<T>& z) {
        if (IsReal(z))
            return Exp(z.Re);
        T r = std::exp(z.Re);
        return Complex<T>(r*std::cos(z.Im), r*std::sin(z.Im));
    }

    /**
     * Computes the natural log of a real number.
     * @param x Real number.
     * @return Real number.
     */
    template <typename T>
    T Log(T x) {
        return std::log(x);
    }
    /**
     * Computes the principal natural log of a complex number \f$z\f$.
     * If we use Euler's formula and rewrite \f$z=re^{\varphi i}\f$, then \f$\log z=\log r+\varphi i\f$. Generally this would be considered
     * a multivalued function as \f$\log z=\{\log r+(\varphi+2\pi k)i\mid k\in\mathbb{Z}\}\f$. But we just take the principal natural log, i.e.,
     * set \f$k=0\f$.
     *
     * (Read more: https://en.wikipedia.org/wiki/Complex_number#Complex_logarithm)
     * @param z Complex number
     * @return Complex number.
     */
    template <typename T>
    Complex<T> Log(const Complex<T>& z) {
        if (IsReal(z))
            return Log(z.Re);
        return Complex<T>(std::log(Abs(z)), Arg(z));
    }
    /**
     * Computes the \f$\log_{base}x\f$.
     * This is done via change-of-base formula, \f$\log_{base}x=\log x/\log base\f$.
     * @param x Real number.
     * @param base Real number.
     * @return Real number.
     */
    template <typename T>
    T Log(T x, T base) {
        return std::log(x)/std::log(base);
    }
    /**
     * Computes the \f$\log_{base}z\f$.
     * This is done via change-of-base formula, \f$\log_{base}z=\log z/\log base\f$.
     * @param z Complex number.
     * @param base Complex number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> Log(const Complex<T>& z, const Complex<T>& base) {
        if (IsReal(z) && IsReal(base))
            return Log(z.Re, base.Re);
        return Log(z)/Log(base);
    }
    /**
     * Computes the \f$\log_{base}z\f$.
     * This is done via change-of-base formula, \f$\log_{base}z=\log z/\log base\f$.
     * @param z Complex number.
     * @param base Real number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> Log(const Complex<T>& z, T base) {
        return Log(z, Complex<T>(base,0));
    }
    /**
     * Computes the \f$\log_{base}x\f$.
     * This is done via change-of-base formula, \f$\log_{base}x=\log x/\log base\f$.
     * @param x Real number.
     * @param base Complex number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> Log(T x, const Complex<T>& base) {
        return Log(Complex<T>(x,0), base);
    }

    /**
     * Computes the real number \f$x^y\f$.
     * Note this will not safeguard expressions such as \f$(-1)^{1/2}\f$.
     * @param x Real number.
     * @param y Real number.
     * @return Real number.
     */
    template <typename T>
    T Pow(T x, T y) {
        return std::pow(x, y);
    }
    /**
     * Computes the complex number \f$z^w\f$.
     * This is done using the exponential, \f$z^w:=e^{w\log z}\f$. Again this is a multivalued function but we take the principal value.
     * @param z Complex number.
     * @param w Complex number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> Pow(const Complex<T>& z, const Complex<T>& w) {
        if (IsReal(z) && IsReal(w))
            return Pow(z.Re, w.Re);
        else
            return Exp(w*Log(z));
    }

    /**
     * Computes the floating-point modulus \f$x\bmod y\f$.
     * @param x Real number.
     * @param y Real number.
     * @return Real number.
     */
    template <typename T>
    T Mod(T x, T y) {
        return std::fmod(x, y);
    }
    /**
     * Returns the complex number \f$(a\bmod y)+(b\bmod y)i\f$ where \f$z=a+bi\f$.
     * @param z Complex number.
     * @param y Real number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> Mod(const Complex<T>& z, T y) {
        return Complex<T>(Mod(z.Re, y), Mod(z.Im, y));
    }
    /**
     * Returns the complex number \f$(a\bmod c)+(b\bmod d)i\f$ where \f$z=a+bi\f$ and \f$w=b+di\f$.
     * @param z Complex number.
     * @param w Complex number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> Mod(const Complex<T>& z, const Complex<T>& w) {
        return Complex<T>(Mod(z.Re, w.Re), Mod(z.Im, w.Im));
    }

    /**
     * Returns the real number \f$\sin x\f$.
     * @param x Real number.
     * @return Real number.
     */
    template <typename T>
    T Sin(T x) {
        return std::sin(x);
    }
    /**
     * Returns the complex number \f$\sin z=\frac{e^{iz}-e^{-iz}}{2i}\f$.
     * @param z Complex number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> Sin(const Complex<T>& z) {
        if (IsReal(z))
            return Sin(z.Re);
        Complex<T> i = Complex<T>::i();
        return (Exp(i*z)-Exp(-i*z))/(2*i);
    }
    /**
     * Returns the real number \f$\cos x\f$.
     * @param x Real number.
     * @return Real number.
     */
    template <typename T>
    T Cos(T x) {
        return std::cos(x);
    }
    /**
     * Returns the complex number \f$\cos z=\frac{e^{iz}+e^{-iz}}{2}\f$.
     * @param z Complex number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> Cos(const Complex<T>& z) {
        if (IsReal(z))
            return Cos(z.Re);
        Complex<T> i = Complex<T>::i();
        return (Exp(i*z)+Exp(-i*z))/2;
    }
    /**
     * Returns the real number \f$\tan x\f$.
     * @param x Real number.
     * @return Real number.
     */
    template <typename T>
    T Tan(T x) {
        return std::tan(x);
    }
    /**
     * Returns the complex number \f$\tan z=\frac{1}{i}\left(\frac{e^{2iz}-1}{e^{2iz}+1}\right)\f$.
     * @param z Complex number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> Tan(const Complex<T>& z) {
        if (IsReal(z))
            return Tan(z.Re);
        Complex<T> i = Complex<T>::i();
        return ((Exp(2*i*z)-1)/(Exp(2*i*z)+1))/i;
    }
    /**
     * Returns the real number \f$\csc x=1/\sin x\f$.
     * @param x Real number.
     * @return Real number.
     */
    template <typename T>
    T Csc(T x) {
        return T(1)/std::sin(x);
    }
    /**
     * Returns the complex number \f$\csc z=\frac{2i}{e^{iz}-e^{-iz}}\f$.
     * @param z Complex number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> Csc(const Complex<T>& z) {
        if (IsReal(z))
            return Csc(z.Re);
        Complex<T> i = Complex<T>::i();
        return (2*i)/(Exp(i*z)-Exp(-i*z));
    }
    /**
     * Returns the real number \f$\sec x=1/\cos x\f$.
     * @param x Real number.
     * @return Real number.
     */
    template <typename T>
    T Sec(T x) {
        return T(1)/std::cos(x);
    }
    /**
     * Returns the complex number \f$\sec z=\frac{2}{e^{iz}+e^{-iz}}\f$.
     * @param z Complex number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> Sec(const Complex<T>& z) {
        if (IsReal(z))
            return Sec(z.Re);
        Complex<T> i = Complex<T>::i();
        return 2/(Exp(i*z)+Exp(-i*z));
    }
    /**
     * Returns the real number \f$\cot x=1/\tan x\f$.
     * @param x Real number.
     * @return Real number.
     */
    template <typename T>
    T Cot(T x) {
        return 1/std::tan(x);
    }
    /**
     * Returns the complex number \f$\cot z=i\left(\frac{e^{2iz}+1}{e^{2iz}-1}\right)\f$.
     * @param z Complex number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> Cot(const Complex<T>& z) {
        if (IsReal(z))
            return Cot(z.Re);
        Complex<T> i = Complex<T>::i();
        return i*((Exp(2*i*z)+1)/(Exp(2*i*z)-1));
    }

    /**
     * Returns the real number \f$\sin^{-1}x\f$.
     * @param x Real number.
     * @return Real number.
     */
    template <typename T>
    T ASin(T x) {
        return std::asin(x);
    }
    /**
     * Returns the complex number \f$\sin^{-1}z=\frac{1}{i}\log\left(iz+|1-z^2|^{1/2}e^{\frac{i}{2}arg(1-z^2)}\right)\f$.
     * @param z Complex number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> ASin(const Complex<T>& z) {
        if (IsReal(z))
            return ASin(z.Re);
        Complex<T> i = Complex<T>::i();
        return Log(i*z+Sqrt(Abs(1-z*z))*Exp(Arg(1-z*z)*(i/2)))/i;
    }
    /**
     * Returns the real number \f$\cos^{-1}x\f$.
     * @param x Real number.
     * @return Real number.
     */
    template <typename T>
    T ACos(T x) {
        return std::acos(x);
    }
    /**
     * Returns the complex number \f$\cos^{-1}z=\frac{1}{i}\log\left(z+i|1-z^2|^{1/2}e^{\frac{i}{2}arg(1-z^2)}\right)\f$.
     * @param z Complex number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> ACos(const Complex<T>& z) {
        if (IsReal(z))
            return ACos(z.Re);
        Complex<T> i = Complex<T>::i();
        return Log(z+i*Sqrt(Abs(1-z*z))*Exp(Arg(1-z*z)*i/2))/i;
    }
    /**
     * Returns the real number \f$\tan^{-1}x\f$.
     * @param x Real number.
     * @return Real number.
     */
    template <typename T>
    T ATan(T x) {
        return std::atan(x);
    }
    /**
     * Returns the complex number \f$\tan^{-1}z=\frac{1}{2i}\log\left(\frac{i-z}{i+z}\right)\f$.
     * @param z Complex number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> ATan(const Complex<T>& z) {
        if (IsReal(z))
            return ATan(z.Re);
        Complex<T> i = Complex<T>::i();
        return Log((i-z)/(i+z))/(2*i);
    }
    /**
     * Returns the real number \f$\csc^{-1}x=\sin^{-1}(1/x)\f$.
     * @param x Real number.
     * @return Real number.
     */
    template <typename T>
    T ACsc(T x) {
        return std::asin(T(1)/x);
    }
    /**
     * Returns the complex number \f$\csc^{-1}z=\sin^{-1}(1/z)\f$.
     * @param z Complex number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> ACsc(const Complex<T>& z) {
        if (IsReal(z))
            return ACsc(z.Re);
        return ASin(T(1)/z);
    }
    /**
     * Returns the real number \f$\sec^{-1}x=\cos^{-1}(1/x)\f$.
     * @param x Real number.
     * @return Real number.
     */
    template <typename T>
    T ASec(T x) {
        return std::acos(T(1)/x);
    }
    /**
     * Returns the complex number \f$\sec^{-1}z=\cos^{-1}(1/z)\f$.
     * @param z Complex number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> ASec(const Complex<T>& z) {
        if (IsReal(z))
            return ASec(z.Re);
        return ACos(T(1)/z);
    }
    /**
     * Returns the real number \f$\cot^{-1}x=\tan^{-1}(1/x)\f$.
     * @param x Real number.
     * @return Real number.
     */
    template <typename T>
    T ACot(T x) {
        return std::atan(T(1)/x);
    }
    /**
     * Returns the complex number \f$\cot^{-1}z=\frac{1}{2i}\log\left(\frac{z+i}{z-i}\right)\f$.
     * @param z Complex number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> ACot(const Complex<T>& z) {
        if (IsReal(z))
            return ACot(z.Re);
        Complex<T> i = Complex<T>::i();
        return Log((z+i)/(z-i))/(2*i);
    }

    /**
     * Returns the real number \f$\sinh x\f$.
     * @param x Real number.
     * @return Real number.
     */
    template <typename T>
    T Sinh(T x) {
        return std::sinh(x);
    }
    /**
     * Returns the complex number \f$\sinh z=\frac{e^z-e^{-z}}{2}\f$.
     * @param z Complex number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> Sinh(const Complex<T>& z) {
        if (IsReal(z))
            return Sinh(z.Re);
        return (Exp(z) - Exp(-z))/2;
    }
    /**
     * Returns the real number \f$\cosh x\f$.
     * @param x Real number.
     * @return Real number.
     */
    template <typename T>
    T Cosh(T x) {
        return std::cosh(x);
    }
    /**
     * Returns the complex number \f$\cosh z=\frac{e^z+e^{-z}}{2}\f$.
     * @param z Complex number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> Cosh(const Complex<T>& z) {
        if (IsReal(z))
            return Cosh(z.Re);
        return (Exp(z) + Exp(-z))/2;
    }
    /**
     * Returns the real number \f$\tanh x\f$.
     * @param x Real number.
     * @return Real number.
     */
    template <typename T>
    T Tanh(T x) {
        return std::tanh(x);
    }
    /**
     * Returns the complex number \f$\tanh z=\frac{e^z-e^{-z}}{e^z+e^{-z}}\f$.
     * @param z Complex number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> Tanh(const Complex<T>& z) {
        if (IsReal(z))
            return Tanh(z.Re);
        return (Exp(z) - Exp(-z))/(Exp(z) + Exp(-z));
    }
    /**
     * Returns the real number \f$csch x = 1/\sinh x\f$.
     * @param x Real number.
     * @return Real number.
     */
    template <typename T>
    T Csch(T x) {
        return T(1)/std::sinh(x);
    }
    /**
     * Returns the complex number \f$csch z=\frac{2}{e^z-e^{-z}}\f$.
     * @param z Complex number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> Csch(const Complex<T>& z) {
        if (IsReal(z))
            return Csch(z.Re);
        return 2/(Exp(z) - Exp(-z));
    }
    /**
     * Returns the real number \f$sech x = 1/\cosh x\f$.
     * @param x Real number.
     * @return Real number.
     */
    template <typename T>
    T Sech(T x) {
        return T(1)/std::cosh(x);
    }
    /**
     * Returns the complex number \f$sech z=\frac{2}{e^z+e^{-z}}\f$.
     * @param z Complex number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> Sech(const Complex<T>& z) {
        if (IsReal(z))
            return Sech(z.Re);
        return 2/(Exp(z) + Exp(-z));
    }
    /**
     * Returns the real number \f$\coth x = 1/\tanh x\f$.
     * @param x Real number.
     * @return Real number.
     */
    template <typename T>
    T Coth(T x) {
        return T(1)/std::tanh(x);
    }
    /**
     * Returns the complex number \f$\coth z=\frac{e^z+e^{-z}}{e^z-e^{-z}}\f$.
     * @param z Complex number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> Coth(const Complex<T>& z) {
        if (IsReal(z))
            return Coth(z.Re);
        return (Exp(z) + Exp(-z))/(Exp(z) - Exp(-z));
    }

    /**
     * Returns the real number \f$\sinh^{-1} x\f$.
     * @param x Real number.
     * @return Real number.
     */
    template <typename T>
    T ASinh(T x) {
        return std::asinh(x);
    }
    /**
     * Returns the complex number \f$\sinh^{-1} z=\log\left(z+|1+z^2|^{1/2}e^{\frac{i}{2}arg(1+z^2)}\right)\f$.
     * @param z Complex number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> ASinh(const Complex<T>& z) {
        if (IsReal(z))
            return ASinh(z.Re);
        Complex<T> i = Complex<T>::i();
        return Log(z+Sqrt(Abs(1+z*z))*Exp(Arg(1+z*z)*i/2));
    }
    /**
     * Returns the real number \f$\cosh^{-1} x\f$.
     * @param x Real number.
     * @return Real number.
     */
    template <typename T>
    T ACosh(T x) {
        return std::acosh(x);
    }
    /**
     * Returns the complex number \f$\cosh^{-1} z=\log\left(z+|z^2-1|^{1/2}e^{\frac{i}{2}arg(z^2-1)}\right)\f$.
     * @param z Complex number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> ACosh(const Complex<T>& z) {
        if (IsReal(z))
            return ACosh(z.Re);
        Complex<T> i = Complex<T>::i();
        return Log(z+Sqrt(Abs(z*z-1))*Exp(Arg(z*z-1)*i/2));
    }
    /**
     * Returns the real number \f$\tanh^{-1} x\f$.
     * @param x Real number.
     * @return Real number.
     */
    template <typename T>
    T ATanh(T x) {
        return std::atanh(x);
    }
    /**
     * Returns the complex number \f$\tanh^{-1} z=\frac{1}{2}\log\left(\frac{1+z}{1-z}\right)\f$.
     * @param z Complex number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> ATanh(const Complex<T>& z) {
        if (IsReal(z))
            return ATanh(z.Re);
        return Log((1+z)/(1-z))/2;
    }
    /**
     * Returns the real number \f$csch^{-1} x=\sinh^{-1}(1/x)\f$.
     * @param x Real number.
     * @return Real number.
     */
    template <typename T>
    T ACsch(T x) {
        return std::asinh(T(1)/x);
    }
    /**
     * Returns the complex number \f$csch^{-1} z=\sinh^{-1}(1/z)\f$.
     * @param z Complex number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> ACsch(const Complex<T>& z) {
        if (IsReal(z))
            return ACsch(z.Re);
        return ASinh(T(1)/z);
    }
    /**
     * Returns the real number \f$sech^{-1} x=\cosh^{-1}(1/x)\f$.
     * @param x Real number.
     * @return Real number.
     */
    template <typename T>
    T ASech(T x) {
        return std::acos(T(1)/x);
    }
    /**
     * Returns the complex number \f$sech^{-1} z=\cosh^{-1}(1/z)\f$.
     * @param z Complex number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> ASech(const Complex<T>& z) {
        if (IsReal(z))
            return ASech(z.Re);
        return ACosh(T(1)/z);
    }
    /**
     * Returns the real number \f$\coth^{-1} x=\tanh^{-1}(1/x)\f$.
     * @param x Real number.
     * @return Real number.
     */
    template <typename T>
    T ACoth(T x) {
        return std::atanh(T(1)/x);
    }
    /**
     * Returns the complex number \f$\coth^{-1} z=\frac{1}{2}\log\left(\frac{z+1}{z-1}\right)\f$.
     * @param z Complex number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> ACoth(const Complex<T>& z) {
        if (IsReal(z))
            return ACoth(z.Re);
        return Log((z+1)/(z-1))/2;
    }

    /**
     * Returns the real number \f$\lceil x\rceil\f$.
     * @param x Real number.
     * @return Real number.
     */
    template <typename T>
    T Ceil(T x) {
        return std::ceil(x);
    }
    /**
     * Returns the complex number \f$\lceil a\rceil+\lceil b\rceil i\f$ where \f$z=a+bi\f$.
     * @param z Complex number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> Ceil(const Complex<T>& z) {
        return Complex<T>(Ceil(z.Re), Ceil(z.Im));
    }
    /**
     * Returns the real number \f$\lfloor x\rfloor\f$.
     * @param x Real number.
     * @return Real number.
     */
    template <typename T>
    T Floor(T x) {
        return std::floor(x);
    }
    /**
     * Returns the complex number \f$\lfloor a\rfloor+\lfloor b\rfloor i\f$ where \f$z=a+bi\f$.
     * @param z Complex number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> Floor(const Complex<T>& z) {
        return Complex<T>(Floor(z.Re), Floor(z.Im));
    }
    /**
     * Rounds \f$x\f$ to the nearest integer.
     * @param x Real number.
     * @return Real number.
     */
    template <typename T>
    T Round(T x) {
        return std::round(x);
    }
    /**
     * Rounds the real and imaginary parts of \f$z\f$ to the nearest integer.
     * @param z Complex number.
     * @return Complex number.
     */
    template <typename T>
    Complex<T> Round(const Complex<T>& z) {
        return Complex<T>(Round(z.Re), Round(z.Im));
    }
}
