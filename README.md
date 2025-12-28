


created: 2017

author copyright marcin matysek (r)ewertyn.PL

marcin.rewertyn@gmail.com

open-source

http://mixedradixfastfouriertransformifft.blogspot.com/
#
http://www.mediafire.com/file/hyz4dbski4w00pb/inverse+fourier+transform+iDFT+ifft+4+methods+in+open+office+++something+extra.ods
#
http://www.mediafire.com/file/59bpnci966ulec9/DFT+FFT+RADIX-2+DIT+algorytm+Transformacja+Fouriera+analitycznie+v3.4.xlsx 
#
    //assumption: if N==signal period in table tab[] then resolution = 1 Hz but N=2^b;
    //when we have signal period 22000 Hz and N=2^15=32768 then MAYBE
    //the fundamental frequency is 0,671 Hz
    //that means that in F(1) is 1*0,671 Hz =0,671 Hz
    //in F(2) is 2*0,671 Hz =1,342 Hz
    //in F(9) is 9*0,671 Hz =6,039 Hz
    //in F(30) is 30*0,671 Hz =20,13 Hz that means in F(30) we will see what value has our signal in 20,13 Hz
//-------------------------------<br />
//when you want to have equal results that are in that false modificator in normal FFT then change this:<br /><br /><br />
/*<br /><br />
 void fun_fourier_transform_FFT_mixed_radix(int N,std::complex<double> tab[])<br />
{<br />
	//new:<br />
    for(int j=0;j<N;j++)<br />
    {<br />
     tab[j].real(tab[j].real()*2/N);<br />
     tab[j].imag(tab[j].imag()*2/N);<br />
    }<br />
}<br />
//and:<br />
void fun_inverse_fourier_transform_FFT_mixed_radix(int N,std::complex<double> tab[])<br />
{<br />
	//new:<br />
    for(int j=0;j<N;j++)<br />
    {<br />
     tab[j].real(tab[j].real()*0.5);<br />
     tab[j].imag(tab[j].imag()*0.5);<br />
    }<br />
}<br />
///<br />
///<br />
//for official modificator that is only in inverse FFT:<br />

 void fun_fourier_transform_FFT_mixed_radix(int N,std::complex<double> tab[])<br />
{ <br />
    for(int j=0;j<N;j++)<br />
    {<br />
      //nothing<br />
    }<br />
}<br />
void fun_inverse_fourier_transform_FFT_mixed_radix(int N,std::complex<double> tab[])<br />
{<br />
    for(int j=0;j<N;j++)<br />
    {<br />
     tab[j].real(tab[j].real()*1/(float)N);//??<br />
     tab[j].imag(tab[j].imag()*1/(float)N);//??<br />
    }<br />
<br />
}<br />
*/<br />

function: void fun_inverse_bits_radix_4(int N,std::complex<double> tab[]) 
is not needed for working. 
In this source code i have my own implementation of inverse_bits that is universal for all radix:

# FFT Inverse Table Permutation

This module provides a C++ implementation for permuting an array of complex numbers. It is a key step in Fast Fourier Transform (FFT) algorithms, specifically performing the index reordering (digit-reversal) necessary for mixed-radix FFT implementations.

## Functionality

The `fun_inverse_table_FFT` function:
1. Accepts the array size `M` and an array of complex numbers `tab`.
2. Calculates the new index order based on the prime factor decomposition (radix).
3. Permutes the elements of `tab` in place using temporary buffers.

## Code

```cpp
#include <complex>
#include <iostream>

// NOTE: This function requires the definition of a helper function 'radix_base'.
// int radix_base(int M, int* stg);

void fun_inverse_table_FFT(int M, std::complex<double> tab[])
{
    int rx5=5,rx4=4,rx3=3,rx2=2,rx7=7,rx11=11;
    int stg[100]={};
    int *tab8 = new int[M];
    int *tab9 = new int[M];
    std::complex<double> *tab11 = new std::complex<double>[M];
    int nb_stages=0;
    int nb1=0;
    int nb2=0;
    int nb3=0;

    nb_stages=6;

    // External function to determine radices
    nb_stages=radix_base(M,stg);

    for(int i=0;i<M;i++)
    {
        //tab9[i]=tab2[i];
        //tab8[i]=tab2[i];
        tab9[i]=i;
        tab8[i]=i;
    }

    nb3=1;
    for(int op=nb_stages;op>=2;op--)
    {
        nb1=stg[op];
        nb3=nb3*stg[op];

        if(op==nb_stages)
        {
            nb2=stg[0];
        }
        else
        {
            nb2=nb2*stg[op+1];
        }

           for(int i=0,n=0,p=0;i<M;i=i+M/nb3,n++)
        {
            if(n>=nb1)
            {
                n=0,p=p+M/nb2;
            }
            for(int j=0,k=0;j<M/nb3;j++,k=k+nb1)
            {
                if(op%2==0)
                {
                    tab8[i+j]=tab9[k+n+p];
                }
                else
                {
                    tab9[i+j]=tab8[k+n+p];
                }
            }
        }
    }

    for(int i=0;i<M;i++)
    {
      tab11[i]=tab[tab8[i]];

    }
    for(int i=0;i<M;i++)
    {
      tab[i]=tab11[i];
    }

    delete [] tab8;
    delete [] tab9;
    delete [] tab11;
}
