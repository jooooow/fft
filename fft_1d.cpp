#include <iostream>
#include <math.h>
#include <complex>
#include <chrono>

#define TIME_BEGIN(a) auto time_begin_##a = std::chrono::high_resolution_clock::now()

#define TIME_END(a)   auto time_end_##a = std::chrono::high_resolution_clock::now();\
					  auto elapse_##a = std::chrono::duration_cast<std::chrono::nanoseconds>(time_end_##a - time_begin_##a);\
                      printf("[%s time measured : %.5f seconds.]\n", #a, elapse_##a.count() * 1e-9)


typedef std::complex<double> complexd;

complexd W(double N, double k)
{
    return std::pow(M_E, -complexd(0, 2 * M_PI * k / N));
}

complexd M(double N, double n)
{
    return std::pow(M_E, complexd(0, 2 * M_PI * n / N));
}

void FFT(double* x, complexd* X, int len)
{
    if(len == 1)
    {
        X[0] = x[0] * W(1, 0);
    }
    else
    {
        double* x_even = new double[len / 2];
        double* x_odd = new double[len / 2];
        complexd* E = new complexd[len / 2];
        complexd* O = new complexd[len / 2];

        for(int k = 0; k < len / 2; k++)
        {
            x_even[k] = x[2 * k];
            x_odd[k] = x[2 * k + 1];
        }
        
        FFT(x_even, E, len / 2);
        FFT(x_odd, O, len / 2);

        for(int k = 0; k < len / 2; k++)
        {
            complexd w = W(len, k);
            X[k] = E[k] + w * O[k];
            X[len / 2 + k] = E[k] - w * O[k];
        }
    }
}

void IFFT_core(complexd* X, complexd* x, int len)
{
    if(len == 1)
    {
        x[0] = X[0] * M(1, 0);
    }
    else
    {
        complexd* X_even = new complexd[len / 2];
        complexd* X_odd = new complexd[len / 2];
        complexd* E = new complexd[len / 2];
        complexd* O = new complexd[len / 2];

        for(int n = 0; n < len / 2; n++)
        {
            X_even[n] = X[2 * n];
            X_odd[n] = X[2 * n + 1];
        }
        
        IFFT_core(X_even, E, len / 2);
        IFFT_core(X_odd, O, len / 2);

        for(int n = 0; n < len / 2; n++)
        {
            complexd m = M(len, n);
            x[n] = E[n] + m * O[n];
            x[len / 2 + n] = E[n] - m * O[n];
        }
    }
}

void IFFT(complexd* X, double* x, int len)
{
    complexd* res = new complexd[len];
    IFFT_core(X,res,len);
    for(int i = 0; i < len; i++)
    {
        x[i] = std::real(res[i]) / len;
    }
}

void DFT(double* x, complexd* X, int len)
{
    for(int k = 0; k < len; k++)
    {
        complexd res(0,0);
        for(int n = 0; n < len; n++)
        {
            res += x[n] * W(len, k * n);
        }
        X[k] = res;
    }
}

void IDFT(complexd* X, double* x, int len)
{
    for(int n = 0; n < len; n++)
    {
        complexd res(0,0);
        for(int k = 0; k < len; k++)
        {
            res += X[k] * M(len, k * n);
        }
        x[n] = std::real(res) / (double)len;
    }
}

int main(int argc, char* argv[])
{
    int N = 32;
    if(argc > 1)
    {
        N = std::atoi(argv[1]);
    }
    printf("fft_1d(N=%d)\n",N);
    double* x = new double[N];
    for(int i = 0; i < N; i++)
    {
        //x[i] = std::sin((double)i / N * 2 * M_PI) + std::cos((double)i / 8 * 2 * M_PI) +  + std::cos((double)i / 6 * 2 * M_PI);;
        x[i] = drand48();
        //std::cout<<x[i]<<",";
    }
    //std::cout<<std::endl;
    complexd* X1 = new complexd[N];
    TIME_BEGIN(dft);
    DFT(x, X1, N);
    TIME_END(dft);

    complexd* X2 = new complexd[N];
    TIME_BEGIN(fft);
    FFT(x, X2, N);
    TIME_END(fft);
    double err1 = 0.0f;
    for(int i = 0; i < N; i++)
    {
        err1 += std::norm(X1[i] - X2[i]);
    }
    printf("fft error = %.8f\n", err1);

    double* x1 = new double[N];
    TIME_BEGIN(idft);
    IDFT(X1, x1, N);
    TIME_END(idft);

    double* x2= new double[N];
    TIME_BEGIN(ifft);
    IFFT(X2, x2, N);
    TIME_END(ifft);

    double err2 = 0.0f;
    for(int i = 0; i < N; i++)
    {
        err2 += std::abs(x1[i] - x2[i]);
    }
    printf("ifft error = %.8f\n", err2);

    return 0;
}