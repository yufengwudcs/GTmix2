//
//  PolynomialUtils.hpp
//  
//
//  Created by Yufeng Wu on 7/14/24.
//

#ifndef PolynomialUtils_hpp
#define PolynomialUtils_hpp

#define _USE_MATH_DEFINES
#include <complex>
#include <vector>
#include <cmath>
#include <iostream>

// The following is taken from: https://algoteka.com/samples/45/polynomial-multiplication-c-plus-plus-o%2528n-log-n%2529-fast-fourier-transform-based-solution

// Assumes the size of x is a power of 2
std::vector<std::complex<double> > fast_fourier_transform(std::vector<std::complex<double> > x, bool inverse = false) {
    
    std::vector<std::complex<double> > w(x.size(), 0.0);  // w[k] = e^{-i2\pi k/N}
    // Precalculate w[k] for faster FFT computation, do it in a 'binary-search' way to provide decent numeric accuracy
    w[0] = 1.0;
    for(int pow_2 = 1; pow_2 < (int)x.size(); pow_2 *= 2) {
        w[pow_2] = std::polar(1.0, 2*M_PI * pow_2/x.size() * (inverse ? 1 : -1) );
    }
    for(int i=3, last=2; i < (int)x.size(); i++) {
        // This way of computing w[k] guarantees that each w[k] is computed in at most log_2 N multiplications
        if(w[i] == 0.0) {
            w[i] = w[last] * w[i-last];
        } else {
            last = i;
        }
    }
    
    for(int block_size = x.size(); block_size > 1; block_size /= 2) {
        // Do the rearrangement for 'recursive' call block by block for each level'
        std::vector<std::complex<double> > new_x(x.size());
        
        for(int start = 0; start < (int)x.size(); start += block_size) {
            for(int i=0; i<block_size; i++) {
                new_x[start + block_size/2 * (i%2) + i/2] = x[start + i];
            }
        }
        x = new_x;
    }
    
    for(int block_size = 2; block_size <= (int)x.size(); block_size *= 2) {
        // Now compute the FFT 'recursively' level by level bottom-up
        std::vector<std::complex<double> > new_x(x.size());
        int w_base_i = x.size() / block_size;  // w[w_base_i] is the e^{-i2\pi / N} value for the N of this level
        
        for(int start = 0; start < (int)x.size(); start += block_size) {
            for(int i=0; i < block_size/2; i++) {
                new_x[start+i]              = x[start+i] + w[w_base_i*i] * x[start + block_size/2 + i];
                new_x[start+block_size/2+i] = x[start+i] - w[w_base_i*i] * x[start + block_size/2 + i];
            }
        }
        x = new_x;
    }
    return x;
}

// ********************************************************

struct Polynomial {
    std::vector<double> a;
    Polynomial(std::vector<double> new_a) : a(new_a) {}
    
    int GetMaxDegree() const { return (int)a.size()-1; }
    
    Polynomial operator*(const Polynomial &r) {  // An FFT based polynomial multiplication algorithm
        int power_2 = 1;
        while(power_2 < (int)(a.size() + r.a.size() - 1)) {
            power_2 *= 2;  // Find the lowest power of 2 that can hold the size of the resulting polynomial
        }
        
        // Perform the FFT
        std::vector<std::complex<double> > x_l(power_2, 0.0);
        std::vector<std::complex<double> > x_r(power_2, 0.0);
        std::vector<std::complex<double> > product(power_2, 0.0);
        
        for(int i=0; i<(int)a.size(); i++) {
            x_l[i] = a[i];
        }
        for(int i=0; i<(int)r.a.size(); i++) {
            x_r[i] = r.a[i];
        }
        x_l = fast_fourier_transform(x_l);  // Do FFT, a.k.a find polynomial values at designated positions
        x_r = fast_fourier_transform(x_r);
        for(int i=0; i<power_2; i++) {  // Multiply the polynomial values at these positions
            product[i] = x_l[i] * x_r[i];
        }
        product = fast_fourier_transform(product, true);  // Do the inverse FFT, a.k.a polynomial interpolation
        
        std::vector<double> result_a(a.size() + r.a.size() - 1);
        for(int i=0; i<(int)result_a.size(); i++) {
            result_a[i] = product[i].real() / power_2;
        }
        // fix the lowest position by direct multiply
        result_a[0] = this->a[0] * r.a[0];
        
        return result_a;
    }
    // for comparison, the slower version
    Polynomial multiplySlow(const Polynomial &rhs) const
    {
        // do direct multiply
        std::vector<double> res( GetMaxDegree()+rhs.GetMaxDegree() + 1 );
        for(int d=0; d<(int)res.size(); ++d)
        {
            res[d] = 0.0;
            for(int k=0; k<=GetMaxDegree() && k <= d; ++k)
            {
                if( d - k <= rhs.GetMaxDegree() )
                {
                    res[d] += a[k]*rhs.a[d-k];
                }
            }
        }
        Polynomial pres(res);
        return pres;
    }
    
    void Dump() const
    {
        std::cout << "Number of coefficients: " << GetNumCoeffs() << " :  ";
        for(int i=0; i<GetNumCoeffs(); ++i)
        {
            std::cout << GetCoeffAt(i) << "  ";
        }
        std::cout << std::endl;
    }
    
    // access polynomial
    int GetNumCoeffs() const { return a.size(); }
    double GetCoeffAt(int i) const { return a[i]; }
};
// ********************************************************
// FFT in log-space

static double LogSumOf2Imp(double x, double y)
{
    return x + log( 1.0 + exp( y-x ) );
}
static double LogSubtractOf2Imp(double x, double y)
{
    return x + log( 1.0 - exp( y-x ) );
}

// Assumes the size of x is a power of 2, x in log-space
std::vector<std::complex<double> > fast_fourier_transform_log(std::vector<std::complex<double> > x, bool inverse = false) {
    
    std::vector<std::complex<double> > w(x.size(), 0.0);  // w[k] = e^{-i2\pi k/N}
    // Precalculate w[k] for faster FFT computation, do it in a 'binary-search' way to provide decent numeric accuracy
    w[0] = 1.0;
    for(int pow_2 = 1; pow_2 < (int)x.size(); pow_2 *= 2) {
        w[pow_2] = std::polar(1.0, 2*M_PI * pow_2/x.size() * (inverse ? 1 : -1) );
    }
    for(int i=3, last=2; i < (int)x.size(); i++) {
        // This way of computing w[k] guarantees that each w[k] is computed in at most log_2 N multiplications
        if(w[i] == 0.0) {
            w[i] = w[last] * w[i-last];
        } else {
            last = i;
        }
    }
    
    for(int block_size = x.size(); block_size > 1; block_size /= 2) {
        // Do the rearrangement for 'recursive' call block by block for each level'
        std::vector<std::complex<double> > new_x(x.size());
        
        for(int start = 0; start < (int)x.size(); start += block_size) {
            for(int i=0; i<block_size; i++) {
                new_x[start + block_size/2 * (i%2) + i/2] = x[start + i];
            }
        }
        x = new_x;
    }
    
    for(int block_size = 2; block_size <= (int)x.size(); block_size *= 2) {
        // Now compute the FFT 'recursively' level by level bottom-up
        std::vector<std::complex<double> > new_x(x.size());
        int w_base_i = x.size() / block_size;  // w[w_base_i] is the e^{-i2\pi / N} value for the N of this level
        
        for(int start = 0; start < (int)x.size(); start += block_size) {
            for(int i=0; i < block_size/2; i++) {
                //new_x[start+i]              = x[start+i] + w[w_base_i*i] * x[start + block_size/2 + i];
                
                //new_x[start+i]              = x[start+i] + log(std::complex(1.0)+w[w_base_i*i] * exp( x[start + block_size/2 + i] - x[start+i] ) );
                
                // tempoary: it seems to give correct results but it doesn't reduce underflow
                new_x[start+i]              = log(exp(x[start+i]) + w[w_base_i*i] * exp( x[start + block_size/2 + i] ) );
              
                // for some reason this doesn't work, need to figure out why. YW: TBD
                //new_x[start+i]              = std::complex(x[start+i].real(), 0.0) +  log(exp(x[start+i]-x[start+i].real() )  + w[w_base_i*i] * exp(  x[start + block_size/2 + i]-x[start+i].real() ) );
                
                //new_x[start+i]              = std::complex(x[start+i].real(), 0.0) +  log(exp( std::complex(0.0, x[start+i].imag() ) ) + w[w_base_i*i] * exp( std::complex( x[start + block_size/2 + i].real()-x[start+i].real(), x[start + block_size/2 + i].imag() ) ) );
                //new_x[start+i]              = x[start+i].real() +  log(exp( std::complex(0.0, x[start+i].imag() ) ) + w[w_base_i*i] * exp( std::complex( x[start + block_size/2 + i] - std::complex( x[start+i].real(), 0.0)  ) ) );
                
                //new_x[start+i]              = x[start+i] + log(std::complex(1.0, 0.0) + w[w_base_i*i] * exp( x[start + block_size/2 + i]-x[start+i] ) );
                
                                                            
                //new_x[start+block_size/2+i] = x[start+i] - w[w_base_i*i] * x[start + block_size/2 + i];
                //new_x[start+block_size/2+i][0] = LogSubtractOf2Imp( x[start+i][0], log(w[w_base_i*i][0]*exp(x[start + block_size/2 + i][0] ) -w[w_base_i*i][1]*exp(x[start + block_size/2 + i][1] ) ) );
                //new_x[start+block_size/2+i][1] = LogSubtractOf2Imp( x[start+i][1], log(w[w_base_i*i][0]*exp(x[start + block_size/2 + i][1] ) + w[w_base_i*i][1]*exp(x[start + block_size/2 + i][0] ) ) );
                
                //new_x[start+block_size/2+i]              = x[start+i] + log(1.0-w[w_base_i*i] * exp( x[start + block_size/2 + i] - x[start+i] ) );
                
                // tempoary: it seems to give correct results but it doesn't reduce underflow
                new_x[start+block_size/2+i]              = log(exp(x[start+i]) - w[w_base_i*i] * exp( x[start + block_size/2 + i] ) );
                
                
                // for some reason this doesn't work, need to figure out why. YW: TBD
                //new_x[start+block_size/2+i]              = std::complex(x[start+i].real(), 0.0) +  log(exp(x[start+i]-x[start+i].real() )  - w[w_base_i*i] * exp(  x[start + block_size/2 + i]-x[start+i].real()));
                
                //new_x[start+block_size/2+i]              = x[start+i].real() +  log(exp( std::complex(0.0, x[start+i].imag() ) ) - w[w_base_i*i] * exp( std::complex( x[start + block_size/2 + i] - std::complex( x[start+i].real(), 0.0)  ) ) );
                //new_x[start+block_size/2+i]              = std::complex(x[start+i].real(), 0.0) +  log(exp( std::complex(0.0, x[start+i].imag() ) ) - w[w_base_i*i] * exp( std::complex( x[start + block_size/2 + i].real()-x[start+i].real(), x[start + block_size/2 + i].imag() ) ) );
                                                                                          
                //new_x[start+block_size/2+i]              = x[start+i] + log( std::complex(1.0, 0.0) - w[w_base_i*i] * exp( x[start + block_size/2 + i] - x[start+i] ) );
                                                        
            }
        }
        x = new_x;
    }
    return x;
}


// coefficient represented as logs
struct LogPolynomial {
    std::vector<double> a;
    LogPolynomial(std::vector<double> new_a) : a(new_a) {}
    
    int GetMaxDegree() const { return (int)a.size()-1; }
    
    LogPolynomial operator*(const LogPolynomial &r) {  // An FFT based polynomial multiplication algorithm
        int power_2 = 1;
        while(power_2 < (int)(a.size() + r.a.size() - 1)) {
            power_2 *= 2;  // Find the lowest power of 2 that can hold the size of the resulting polynomial
        }
        std::cout << "power_2: " << power_2 << std::endl;
        
        // Perform the FFT
        std::vector<std::complex<double> > x_l(power_2, -1e100);
        std::vector<std::complex<double> > x_r(power_2, -1e100);
        std::vector<std::complex<double> > product(power_2, -1e100);
        
        for(int i=0; i<(int)a.size(); i++) {
            x_l[i] = a[i];
        }
        for(int i=0; i<(int)r.a.size(); i++) {
            x_r[i] = r.a[i];
        }
        x_l = fast_fourier_transform_log(x_l);  // Do FFT, a.k.a find polynomial values at designated positions
        x_r = fast_fourier_transform_log(x_r);
        for(int i=0; i<power_2; i++) {  // Multiply the polynomial values at these positions
            product[i] = x_l[i] + x_r[i];
        }
        product = fast_fourier_transform_log(product, true);  // Do the inverse FFT, a.k.a polynomial interpolation
        
        std::vector<double> result_a(a.size() + r.a.size() - 1);
        for(int i=0; i<(int)result_a.size(); i++) {
            result_a[i] = product[i].real() - log( power_2);
        }
        // fix the lowest position by direct multiply
        //result_a[0] = this->a[0] + r.a[0];
        
        return result_a;
    }
    
    void Dump() const
    {
        std::cout << "Number of coefficients: " << GetNumCoeffs() << " :  ";
        for(int i=0; i<GetNumCoeffs(); ++i)
        {
            std::cout << GetCoeffAt(i) << "  ";
        }
        std::cout << std::endl;
    }
    
    // access polynomial
    int GetNumCoeffs() const { return a.size(); }
    double GetCoeffAt(int i) const { return a[i]; }
};


// ********************************************************
// another implementation
using namespace std;

typedef long double DOUBLE;
typedef complex<DOUBLE> COMPLEX;
typedef vector<DOUBLE> VD;
typedef vector<COMPLEX> VC;

struct FFT {
  VC A;
  int n, L;

  int ReverseBits(int k) {
    int ret = 0;
    for (int i = 0; i < L; i++) {
      ret = (ret << 1) | (k & 1);
      k >>= 1;
    }
    return ret;
  }

  void BitReverseCopy(VC a) {
    for (n = 1, L = 0; n < a.size(); n <<= 1, L++) ;
    A.resize(n);
    for (int k = 0; k < n; k++)
      A[ReverseBits(k)] = a[k];
  }
  
  VC DFT(VC a, bool inverse) {
    BitReverseCopy(a);
    for (int s = 1; s <= L; s++) {
      int m = 1 << s;
      COMPLEX wm = exp(COMPLEX(0, 2.0 * M_PI / m));
      if (inverse) wm = COMPLEX(1, 0) / wm;
      for (int k = 0; k < n; k += m) {
    COMPLEX w = 1;
    for (int j = 0; j < m/2; j++) {
      COMPLEX t = w * A[k + j + m/2];
      COMPLEX u = A[k + j];
      A[k + j] = u + t;
      A[k + j + m/2] = u - t;
      w = w * wm;
    }
      }
    }
    if (inverse) for (int i = 0; i < n; i++) A[i] /= n;
    return A;
  }

  // c[k] = sum_{i=0}^k a[i] b[k-i]
  VD Convolution(VD a, VD b) {
    int L = 1;
    while ((1 << L) < a.size()) L++;
    while ((1 << L) < b.size()) L++;
    int n = 1 << (L+1);

    VC aa, bb;
    for (size_t i = 0; i < n; i++) aa.push_back(i < a.size() ? COMPLEX(a[i], 0) : 0);
    for (size_t i = 0; i < n; i++) bb.push_back(i < b.size() ? COMPLEX(b[i], 0) : 0);
    
    VC AA = DFT(aa, false);
    VC BB = DFT(bb, false);
    VC CC;
    for (size_t i = 0; i < AA.size(); i++) CC.push_back(AA[i] * BB[i]);
    VC cc = DFT(CC, true);

    VD c;
    for (int i = 0; i < a.size() + b.size() - 1; i++) c.push_back(cc[i].real());
    return c;
  }

};

// ********************************************************
// another implementation (log-space)

struct FFTLog {
  VC A;
  int n, L;

  int ReverseBits(int k) {
    int ret = 0;
    for (int i = 0; i < L; i++) {
      ret = (ret << 1) | (k & 1);
      k >>= 1;
    }
    return ret;
  }

  void BitReverseCopy(VC a) {
    for (n = 1, L = 0; n < a.size(); n <<= 1, L++) ;
    A.resize(n);
    for (int k = 0; k < n; k++)
      A[ReverseBits(k)] = a[k];
  }
  
  VC DFT(VC a, bool inverse) {
    BitReverseCopy(a);
    for (int s = 1; s <= L; s++) {
      int m = 1 << s;
      COMPLEX wm = exp(COMPLEX(0, 2.0 * M_PI / m));
      if (inverse) wm = COMPLEX(1, 0) / wm;
      for (int k = 0; k < n; k += m) {
    COMPLEX w = 1;
    for (int j = 0; j < m/2; j++) {
      //COMPLEX t = w * exp(A[k + j + m/2]);
      //COMPLEX u = exp(A[k + j]);
      //A[k + j] = log(u + t);
      //A[k + j + m/2] = log(u - t);
        COMPLEX u = A[k + j];
        A[k + j] = u + log(COMPLEX(1.0) + w*exp( A[k + j + m/2] - u ));
        A[k + j + m/2] = u + log(COMPLEX(1.0) - w*exp( A[k + j + m/2] - u ));
      w = w * wm;
    }
      }
    }
    if (inverse) for (int i = 0; i < n; i++) A[i] -= log(n);
    return A;
  }

  // c[k] = sum_{i=0}^k a[i] b[k-i]
  VD Convolution(VD a, VD b) {
    int L = 1;
    while ((1 << L) < a.size()) L++;
    while ((1 << L) < b.size()) L++;
    int n = 1 << (L+1);

    VC aa, bb;
    for (size_t i = 0; i < n; i++) aa.push_back(i < a.size() ? COMPLEX(a[i], 0) : 0);
    for (size_t i = 0; i < n; i++) bb.push_back(i < b.size() ? COMPLEX(b[i], 0) : 0);
    
    VC AA = DFT(aa, false);
    VC BB = DFT(bb, false);
    VC CC;
    for (size_t i = 0; i < AA.size(); i++) CC.push_back(AA[i] + BB[i]);
    VC cc = DFT(CC, true);

    VD c;
    for (int i = 0; i < a.size() + b.size() - 1; i++) c.push_back(cc[i].real());
    return c;
  }

};


/*
int main() {
  double a[] = {1, 3, 4, 5, 7};
  double b[] = {2, 4, 6};

  FFT fft;
  VD c = fft.Convolution(VD(a, a + 5), VD(b, b + 3));

  // expected output: 2 10 26 44 58 58 42
  for (int i = 0; i < c.size(); i++) cerr << c[i] << " ";
  cerr << endl;
  
  return 0;
}
*/


// utilities for convolution
void UtilConvoluteTwoVecs(const std::vector<double> &vec1, const std::vector<double> &vec2, std::vector<double> &vecProduct)
{
    Polynomial p1(vec1), p2(vec2);
    Polynomial pprod = p1 * p2;
    vecProduct = pprod.a;
}

// log-based
void UtilConvoluteTwoLogVecs(const std::vector<double> &vec1, const std::vector<double> &vec2, std::vector<double> &vecProduct)
{
    // the numbers are in logs
    // first get the min of each vec of logs, and then use the exponent of the differences
    double vMin1 = vec1[0];
    double vMin2 = vec2[0];
#if 0
    // use median
    std::vector<double> vec1Sort = vec1, vec2Sort = vec2;
    std::sort(vec1Sort.begin(), vec1Sort.end());
    std::sort(vec2Sort.begin(), vec2Sort.end());
    vMin1 = vec1Sort[ vec1Sort.size()/2 ];
    vMin2 = vec1Sort[ vec2Sort.size()/2 ];
#endif
//#if 0
    for(unsigned int i=0; i<vec1.size(); ++i)
    {
        if( vec1[i] < vMin1)
        {
            vMin1 = vec1[i];
        }
    }
    for(unsigned int i=0; i<vec2.size(); ++i)
    {
        if( vec2[i] < vMin2)
        {
            vMin2 = vec2[i];
        }
    }
//#endif
    std::vector<long double> vecExpDiff1(vec1.size()), vecExpDiff2(vec2.size());
    for(unsigned int i=0; i<vec1.size(); ++i)
    {
        vecExpDiff1[i] = std::exp(vec1[i] - vMin1);
    }
    for(unsigned int i=0; i<vec2.size(); ++i)
    {
        vecExpDiff2[i] = std::exp(vec2[i] - vMin2);
    }
    // do a multiplication
    std::vector<double> vecExpProd;
    //UtilConvoluteTwoVecs(vecExpDiff1, vecExpDiff2, vecExpProd);
    FFT fft;
    VD c = fft.Convolution(vecExpDiff1, vecExpDiff2);
    vecProduct.resize( c.size() );
    for(unsigned int i=0; i<c.size(); ++i)
    {
        double valUse = c[i];
        if( c[i] < 0.0)
        {
            // assign a very small number
            valUse = 1.0e-300;
        }
        vecProduct[i] = log(valUse) + vMin1 + vMin2;
    }
    //
    //vecProduct.resize( vecExpProd.size() );
    //for(unsigned int i=0; i<vecExpProd.size(); ++i)
    //{
    //    vecProduct[i] = log(vecExpProd[i]) + vMin1 + vMin2;
    //}
    // for lowest term, just do a direct multiplication
    if( vec1.size() >0 && vec2.size() > 0)
    {
        vecProduct[0] = vec1[0] + vec2[0];
    }
}

// another attempt for log-based polynomial
// coefficient represented as logs
class LogPolynomial2 {
public:
    LogPolynomial2()
        : degree{0}, coeffs{nullptr} {
            this->coeffs = new double[1];
        //this->coeffs[0] = 0;
            this->coeffs[0] = MIN_VAL;
    }

    LogPolynomial2(long degree)
        : degree{degree}, coeffs{nullptr} {
        this->coeffs = new double[degree + 1];
        for(auto i = 0; i <= this->degree; ++i)
            //this->coeffs[i] = 0;
            this->coeffs[i] = MIN_VAL;
    }

    LogPolynomial2(const LogPolynomial2 &src)
        : degree{src.degree}, coeffs{nullptr} {
        this->coeffs = new double[src.degree + 1];
        for(auto i = 0; i <= this->degree; ++i)
            this->coeffs[i] = src.coeffs[i];
    }

    LogPolynomial2(LogPolynomial2 &&src)
        : degree(src.degree), coeffs(src.coeffs) {
        src.degree = -1;
        src.coeffs = nullptr;
    }

    ~LogPolynomial2() {
        delete this->coeffs;
    }

    LogPolynomial2 &operator=(const LogPolynomial2 &rhs) {
        this->degree = rhs.degree;
        delete [] this->coeffs;
        this->coeffs = new double[rhs.degree + 1];
        for(auto i = 0; i <= this->degree; ++i)
            this->coeffs[i] = rhs.coeffs[i];
        return *this;
    }

    LogPolynomial2 &operator=(LogPolynomial2 &&rhs) {
        this->degree = rhs.degree;
        rhs.degree = -1;
        delete [] this->coeffs;
        this->coeffs = rhs.coeffs;
        rhs.coeffs = nullptr;
        return *this;
    }

    bool operator==(const LogPolynomial2 &rhs) const {
        bool res = (this->degree == rhs.degree);

        for(auto i = 0; res && i <= this->degree; ++i)
            res = (this->coeffs[i] == rhs.coeffs[i]);

        return res;
    }

    LogPolynomial2 operator+(const LogPolynomial2 &rhs) const {
        long size = (this->degree > rhs.degree) ? this->degree : rhs.degree;
        LogPolynomial2 temp(size);
        for(int i = 0; i <= size; ++i) {
            if(i <= this->degree)
                //temp.coeffs[i] += this->coeffs[i];
                temp.coeffs[i] = GetLogSumOfTwo( temp.coeffs[i], this->coeffs[i] );
            if(i <= rhs.degree)
                //temp.coeffs[i] += rhs.coeffs[i];
                temp.coeffs[i] = GetLogSumOfTwo( temp.coeffs[i], rhs.coeffs[i] );
        }
        
        long new_degree = temp.degree;
        //while(new_degree >= 0 && temp.coeffs[new_degree] == 0)
        while(new_degree >= 0 && temp.coeffs[new_degree] <= MIN_VAL)
            new_degree--;

        LogPolynomial2 res(new_degree);
        for(auto i = 0; i <= new_degree; ++i)
            res.coeffs[i] = temp.coeffs[i];
        return res;
    }

    LogPolynomial2 operator-(const LogPolynomial2 &rhs) const {
        long size = (this->degree > rhs.degree) ? this->degree : rhs.degree;
        LogPolynomial2 temp(size);
        for(int i = 0; i <= size; ++i) {
            if(i <= this->degree)
                //temp.coeffs[i] += this->coeffs[i];
                temp.coeffs[i] = GetLogSumOfTwo(temp.coeffs[i], this->coeffs[i]);
            if(i <= rhs.degree)
                //temp.coeffs[i] -= rhs.coeffs[i];
                temp.coeffs[i] = GetLogDiffOfTwo( temp.coeffs[i], rhs.coeffs[i] );
        }
        
        long new_degree = temp.degree;
        //while(new_degree >= 0 && temp.coeffs[new_degree] == 0)
        while(new_degree >= 0 && temp.coeffs[new_degree] <= MIN_VAL )
            new_degree--;

        LogPolynomial2 res(new_degree);
        for(auto i = 0; i <= new_degree; ++i)
            res.coeffs[i] = temp.coeffs[i];
        return res;
    }

    LogPolynomial2 operator*(const LogPolynomial2 &rhs) const {
        LogPolynomial2 res(this->degree + rhs.degree);

        for(auto i = 0; i <= this->degree; ++i)
        {
            for(auto j = 0; j <= rhs.degree; ++j)
            {
//cout << "Before Update: res.coeffs[" << i+j << "]: " << res.coeffs[i+j] << ", coeffs[" << i << "]:" << this->coeffs[i] << ", rhs.coeffs[" << j << "]:" << rhs.coeffs[j] << endl;
                //res.coeffs[i + j] += this->coeffs[i] * rhs.coeffs[j];
                res.coeffs[i+j] = GetLogSumOfTwo(res.coeffs[i+j], this->coeffs[i] + rhs.coeffs[j] );
                
//cout << "Update res.coeffs[" << i+j << "]: " << res.coeffs[i+j] << endl;
            }
        }
        return res;
    }

    LogPolynomial2 increaseExponent(long n) {
        this->setDegree(this->degree + n);

        for(auto i = this->degree; i >= 0; --i) {
            if(i < n) //this->coeffs[i] = 0;
                this->coeffs[i] = MIN_VAL;
            else this->coeffs[i] = this->coeffs[i - n];
        }
        return *this;
    }

    LogPolynomial2 karatsuba(const LogPolynomial2 &rhs) const {
        if(this->degree <= 1 || rhs.degree >= 1)
            return ((*this) * rhs);

        LogPolynomial2 res(this->degree + rhs.degree);

        long half = floor(this->degree / 2);
        LogPolynomial2 A0(half - 1), A1(this->degree - half), B0(half - 1), B1(this->degree - half);
        for(auto i = 0; i < half; ++i) {
            A0.setCoeff(this->coeffs[i], i);
            B0.setCoeff(rhs.coeffs[i], i);
            A1.setCoeff(this->coeffs[i + half], i);
            B1.setCoeff(rhs.coeffs[i + half], i);
            if(this->degree % 2 == 1 && i + 1 == half) {
                A1.setCoeff(this->coeffs[i + half], i);
                B1.setCoeff(rhs.coeffs[i + half], i);
            }
        }

        LogPolynomial2 U = A0.karatsuba(B0);
        LogPolynomial2 Y = (A0 + A1).karatsuba(B0 + B1);
        LogPolynomial2 Z = A1.karatsuba(B1);

        LogPolynomial2 W = Y - U - Z;
        res = U + (W.increaseExponent(half)) + (Z.increaseExponent(this->degree));
        return res;
    }

    void displayCoeffs() const {
        for(auto i = 0; i <= this->degree; ++i)
            std::cout << this->coeffs[i] << " ";
        std::cout << std::endl;
    }

    void inputPolynomial() {
        for(auto i = 0; i <= this->degree; ++i)
            std::cin >> this->coeffs[i];
    }

    void randomPolynomial() {
        for(auto i = 0; i <= this->degree; ++i)
            this->coeffs[i] = rand() % 100;
    }

    double getCoeff(long n) {
        return this->coeffs[n];
    }

    void setCoeff(long n, int coeff) {
        this->coeffs[n] = coeff;
    }

    long getDegree() {
        return this->degree;
    }

    void setDegree(long n) {
        long prev_degree = this->degree;
        double prev_coeffs[prev_degree + 1];
        for(auto i = 0; i <= prev_degree; ++i)
            prev_coeffs[i] = this->coeffs[i];

        this->degree = n;
        this->coeffs = new double[this->degree + 1];

        for(auto i = 0; i <= this->degree; ++i) {
            if(i < prev_degree)
                //this->coeffs[i] = 0;
                this->coeffs[i] = MIN_VAL;
            else this->coeffs[i] = prev_coeffs[i];
        }
    }

private:
    // utility
    double GetLogSumOfTwo(double logv1, double logv2) const
    {
        // calculate it directly
        if( logv1 < logv2 )
        {
            std::swap(logv1, logv2);
        }
        return logv1 +  log( 1.0 + exp( logv2-logv1 ) );
    }
    double GetLogDiffOfTwo(double logv1, double logv2) const
    {
        // calculate it directly
        if( logv1 < logv2 )
        {
            std::swap(logv1, logv2);
        }
        return logv1 +  log( 1.0 - exp( logv2-logv1 ) );
    }
    const double MIN_VAL = -1.0e100;
    
    long degree;
    double *coeffs;

};

#endif /* PolynomialUtils_hpp */
