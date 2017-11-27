/*
 * Apollonius.cpp
 *
 *  Created on: 2017/11/26
 *      Author: stomo
 */

#include <iostream>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <complex>

#define EPS 0.000001  // 停止判定
#define MAXIT 30      // 最大反復回数

// Zeros
std::vector< std::complex<double> > Zrs {
	std::complex<double> (  0.0,  1.0 ),
	std::complex<double> (  1.0,  2.0 ),
	std::complex<double> ( -1.0,  2.0 ),
	std::complex<double> (  3.0, -3.0 ),
	std::complex<double> ( -3.0, -3.0 )
};

// Coefficients
std::vector< std::complex<double> > Cef {
	std::complex<double> (  1.0,   0.0 ),  // z^5
	std::complex<double> (  0.0,   1.0 ),  // Z^4
	std::complex<double> (  3.0,   0.0 ),  // Z^3
	std::complex<double> (  0.0,  41.0 ),  // z^2
	std::complex<double> (132.0,   0.0 ),  // z^1
	std::complex<double> (  0.0, -90.0 )   // z^0
};

std::vector<double> Gam( Zrs.size() );  // Γ
std::vector<double> Alp( Zrs.size() );  // α

// Horner method for polynomial
template<typename T> void Horner( std::vector< std::complex<T> > cf, std::complex<T> z,
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
template <typename T> std::complex<T> vc( const int K, std::complex<T> z )
{
	std::complex<T> tmp = std::complex<T> (0.0,0.0);;

	for (auto itr = Zrs.begin(); itr < Zrs.end(); ++itr )
	{
		std::complex<T> vf, df;
		Horner( Cef, *itr, vf, df );

		// tmp *= (1/f'(z_i) (-1 / (z_i -z)^{K+1})
		tmp += ( (T)(1.0) / df )*( (T)(-1.0) / pow( (*itr - z), (T)(K+1) ));
	}
	return tmp;
}

template <typename T> std::complex<T> Nourein( const int p, std::complex<T> z, int &count, T &er )
{
	assert(p>=2);

	std::complex<T> vf, df;
	Horner( Cef, z, vf, df );
	count = 0;

	while ((count < MAXIT) && (abs(vf) > EPS))
	{
		z += vc(p-2,z) / vc(p-1,z);
		Horner( Cef, z, vf, df );
		count++;
	}
	er = abs(vf);

	return z;
}

template <typename T> void SetGamma( std::vector<T> &Gam )
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

template <typename T> T fAlp( const int p, const T Gamma, T alp )
{
	assert(p>=2);

	return ((T)(Zrs.size()) -1.0) * Gamma * pow(alp,(T)(p-1)) - (1.0-alp)/(1.0+3.0*alp);
}

template <typename T> T dAlp( const int p, const T Gamma, T alp )
{
	assert(p>=2);

	T tmp = ((T)(Zrs.size()) -1.0) * Gamma * ((T)(p) -1.0) * pow(alp,(T)(p-2));
	return  tmp + 4.0/(1.0+6.0*alp+9.0*alp*alp);
}

template <typename T> void GetAlpha( const int P, const std::vector<T> Gam, std::vector<T> &Alp )
{
	for (int i=0; i<Zrs.size(); i++)
	{
		T alp = (T)(1.0);
		int count = 0;

		while ((count < MAXIT) && (abs(fAlp(P,Gam[i],alp)) > EPS))
		{
			alp -= fAlp(P,Gam[i],alp) / dAlp(P,Gam[i],alp);
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
	int P;  // Order of Nourein method

	SetGamma( Gam );       // Γ

    for(int P=2; P<=128; P++)
    {
    	GetAlpha( P, Gam, Alp );

    	std::cout << P << ", ";
    	for(auto itr = Alp.begin(); itr < Alp.end(); itr++ )
    	{
    		std::cout << *itr << ", ";
    	}
    	std::cout << std::endl;
    }

	return EXIT_SUCCESS;
}
