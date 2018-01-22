/*
 * Apollonius.cpp
 *
 *  Created on: 2017/11/26
 *      Author: stomo
 */

#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <vector>
#include <complex>

#define EPS 0.000001  // 停止判定
#define MAXIT 30      // 最大反復回数
#define N 5         // 多項式の次数

// Set Zeros
template <typename T>
void Set_Zeros( std::vector< std::complex<T> > &Zrs )
{
	Zrs[0] = std::complex<T> (  0.0,  1.0 );
	Zrs[1] = std::complex<T> (  1.0,  2.0 );
	Zrs[2] = std::complex<T> ( -1.0,  2.0 );
	Zrs[3] = std::complex<T> (  3.0, -3.0 );
	Zrs[4] = std::complex<T> ( -3.0, -3.0 );
}

// Set Coefficients
template <typename T>
void Set_Coef( std::vector< std::complex<T> > &Cef )
{
	Cef[0] = std::complex<T> (  1.0,   0.0 );  // z^5
	Cef[1] = std::complex<T> (  0.0,   1.0 );  // Z^4
	Cef[2] = std::complex<T> (  3.0,   0.0 );  // Z^3
	Cef[3] = std::complex<T> (  0.0,  41.0 );  // z^2
	Cef[4] = std::complex<T> (132.0,   0.0 );  // z^1
	Cef[5] = std::complex<T> (  0.0, -90.0 );  // z^0
}

// Horner method for polynomial
template<typename T>
void Horner( const std::vector< std::complex<T> > Cef, const std::complex<T> z,
		std::complex<T> &vf, std::complex<T> &df )
{
	vf = Cef[0];
	df = std::complex<T> (0.0,0.0);
	std::complex<T> tmp;

    for(auto itr = Cef.begin()+1; itr < Cef.end(); ++itr)
    {
    	tmp = vf;
    	vf = vf*z + *itr;
    	df = df*z + tmp;
    }
}

// Nourein subfunction
template <typename T>
std::complex<T> vc( const int P, const std::vector< std::complex<T> > Zrs, const std::vector< std::complex<T> > Cef,
		const std::complex<T> z )
{
	std::complex<T> tmp = std::complex<T> (0.0,0.0);;

	for (auto itr = Zrs.begin(); itr < Zrs.end(); ++itr )
	{
		std::complex<T> vf, df;
		Horner( Cef, *itr, vf, df );

		// tmp *= (1/f'(z_i) (-1 / (z_i -z)^{P+1})
		tmp += ( (T)(1.0) / df )*( (T)(-1.0) / pow( (*itr - z), (T)(P+1) ));
	}
	return tmp;
}

template <typename T>
std::complex<T> Nourein( const int P, const std::vector< std::complex<T> > Zrs, const std::vector< std::complex<T> > Cef,
		const std::complex<T> z, int &count, T &er )
{
	assert(P>=2);

	std::complex<T> vf, df;
	Horner( Cef, z, vf, df );
	count = 0;

	while ((count < MAXIT) && (abs(vf) > EPS))
	{
		z += vc(P-2,z) / vc(P-1,z);
		Horner( Cef, z, vf, df );
		count++;
	}
	er = abs(vf);

	return z;
}

template <typename T>
void SetGamma( const std::vector< std::complex<T> > Zrs, const std::vector< std::complex<T> > Cef, std::vector<T> &Gam )
{
	for (int i=0; i<Zrs.size(); i++)
	{
		T max = (T)(0.0);
		for (int j=0; j<Zrs.size(); j++)
		{
			if (i != j)
			{
				std::complex<T> vf, dfi, dfj;

				Horner( Cef, Zrs[i], vf, dfi );
				Horner( Cef, Zrs[j], vf, dfj );

				if (max < abs(dfi / dfj))
				{
					max = abs(dfi / dfj);
				}
			}
		}
		Gam[i] = max;
	}
}

template <typename T>
T fAlp( const int P,  const std::vector< std::complex<T> > Zrs, const T Gamma, const T alp )
{
	assert(P>=2);

	return ((T)(Zrs.size()) -1.0) * Gamma * pow(alp,(T)(P-1)) - (1.0-alp)/(1.0+3.0*alp);
}

template <typename T>
T dAlp( const int P,  const std::vector< std::complex<T> > Zrs, const T Gamma, const T alp )
{
	assert(P>=2);

	T tmp = ((T)(Zrs.size()) -1.0) * Gamma * ((T)(P) -1.0) * pow(alp,(T)(P-2));
	return  tmp + 4.0/(1.0 + 6.0*alp + 9.0*alp*alp);
}

template <typename T>
void GetAlpha( const int P, const std::vector< std::complex<T> > Zrs, const std::vector<T> Gam, std::vector<T> &Alp )
{
	for (int i=0; i<Zrs.size(); i++)
	{
		T alp = (T)(1.0);
		int count = 0;

		while ((count < MAXIT) && (abs(fAlp(P,Zrs,Gam[i],alp)) > EPS))
		{
			alp -= fAlp(P,Zrs,Gam[i],alp) / dAlp(P,Zrs,Gam[i],alp);
			count++;
		}
		if (count == MAXIT)
		{
			std::cerr << "No convergence in GetAlpha\n";
			std::exit (EXIT_FAILURE);
		}
		Alp[i] = alp;
	}
}

int main(int argc, char *argv[])
{
	if (argc <2)
	{
		std::cerr << "Usage: a.out [order of Nourein's method]\n";
		exit (EXIT_FAILURE);
	}
	const int P = atoi(argv[1]);    // Order of Nourein method

	std::vector< std::complex<double> > Zrs( N );       // Zeros
	std::vector< std::complex<double> > Cef( N+1);   // Coefficients

	Set_Zeros( Zrs );
	Set_Coef( Cef );

	std::vector<double> Gam( N );  // Γ
	std::vector<double> Alp( N );    // α

	SetGamma( Zrs, Cef, Gam );
	GetAlpha( P, Zrs, Gam, Alp );

	std::cout << P << ", ";
	for(auto itr = Alp.begin(); itr < Alp.end(); itr++ )
	{
//		std::cout << *itr << ", ";
		std::cout << (1.0 - (*itr)) / (log((double)(P)) / (double)(P) ) << ", ";
	}
	std::cout << std::endl;

	for (int p=2; p<=1024; p=p*2)
	{
		GetAlpha( p, Zrs, Gam, Alp );
		std::cout << p << ", ";
		for(auto itr = Alp.begin(); itr < Alp.end(); itr++ )
		{
			std::cout << (1.0 - (*itr)) / (log((double)(p)) / (double)(p) ) << ", ";
		}
		std::cout << std::endl;
	}
	return EXIT_SUCCESS;
}
