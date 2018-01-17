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
#define N 5           // 多項式の次数

// Zrs: ゼロ点ベクトル
template <typename T>
void Set_Zeros( std::vector< std::complex<T> > &Zrs )
{
	Zrs[0] = std::complex<T> (  0.0,  1.0 );
	Zrs[1] = std::complex<T> (  1.0,  2.0 );
	Zrs[2] = std::complex<T> ( -1.0,  2.0 );
	Zrs[3] = std::complex<T> (  3.0, -3.0 );
	Zrs[4] = std::complex<T> ( -3.0, -3.0 );
}

// Cef: 多項式の係数ベクトル
// z^5 + i z^4 + 3 z^3 + 41i z^2 + 132 z -90i <- Phoenix
template <typename T>
void Set_Coef( std::vector< std::complex<T> > &Cef )
{
	Cef[0] = std::complex<T> (  1.0,   0.0 );    // z^5
	Cef[1] = std::complex<T> (  0.0,   1.0 );    // Z^4
	Cef[2] = std::complex<T> (  3.0,   0.0 );    // Z^3
	Cef[3] = std::complex<T> (  0.0,  41.0 );    // z^2
	Cef[4] = std::complex<T> (132.0,   0.0 );    // z^1
	Cef[5] = std::complex<T> (  0.0, -90.0 );    // z^0
}

// Horner 法
//  vf = f(z)
//  df = f'(z)
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

// Gam: Γ_i = max_{j \neq i} | f'(z_i) / f'(z_j) |
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

// Ap比の方程式
//   次数P、ゼロ点i毎にAp比が求められる
template <typename T>
T fAlp( const int P,  const std::vector< std::complex<T> > Zrs, const T Gamma, const T alp )
{
	assert(P>=2);

	return ((T)(Zrs.size()) -1.0) * Gamma * pow(alp,(T)(P-1)) - (1.0-alp)/(1.0+3.0*alp);
}

// Ap比の方程式の導関数
template <typename T>
T dAlp( const int P,  const std::vector< std::complex<T> > Zrs, const T Gamma, const T alp )
{
	assert(P>=2);

	T tmp = ((T)(Zrs.size()) -1.0) * Gamma * ((T)(P) -1.0) * pow(alp,(T)(P-2));
	return  tmp + 4.0/(1.0 + 6.0*alp + 9.0*alp*alp);
}

// Newton法でAp比を求める
// Alp: Ap比ベクトル
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

// Apollonius 円クラス
template <typename T>
class ApCircle
{
protected:
	std::complex<T> c_; // center
	T r_;  // radius

public:
	ApCircle( const T alp, const std::complex<T> base, const std::complex<T> oppo )
	{
		assert( alp > 0.0 );

		c_ = (base - alp*alp*oppo) / (1 - alp*alp);
		r_ = alp*abs(base - oppo) / (1 - alp*alp);
	}

	// Getters
	std::complex<T> center() const { return c_; }
	T radius() const { return r_; }
};

template <typename T>
class ApArc : public ApCircle<T>
{
	protected:
		T stA_;             // 始点の偏角
		T edA_;             // 終点の偏角

	public:
		// Constructor
		ApArc( const T alp, const std::complex<T> base, const std::complex<T> oppo, const T stA, const T edA ) : ApCircle<T>(alp,base,oppo)
		{
			assert( stA > -M_PI/2.0 && stA < M_PI/2.0 );
			assert( edA > -M_PI/2.0 && edA < M_PI/2.0 );

			stA_ = stA;
			edA_ = edA;
		}

		// Getters
		T startArg() const { return stA_; }
		T endArg() const { return edA_; }

		std::complex<T> startPt() const
		{
			std::complex<T> tmp;

			tmp.real() = (this->c_).real() + cos( stA_ );
			tmp.imag() = (this->c_).imag() + sin( stA_ );

			return tmp;
		}

		std::complex<T> endPt() const
		{
			std::complex<T> tmp;

			tmp.real() = (this->c_).real() + cos( edA_ );
			tmp.imag() = (this->c_).imag() + sin( edA_ );

			return tmp;
		}
};

int main(int argc, char *argv[])
{
	if (argc <2)
	{
		std::cerr << "Usage: a.out [order of Nourein's method]\n";
		exit (EXIT_FAILURE);
	}
	const int P = atoi(argv[1]);    // Order of Nourein method

	std::vector< std::complex<double> > Zrs( N );    // Zeros
	std::vector< std::complex<double> > Cef( N+1);   // Coefficients

	Set_Zeros( Zrs );
	Set_Coef( Cef );

	std::vector<double> Gam( N );  // Γ
	std::vector<double> Alp( N );  // Apollonius比

	SetGamma( Zrs, Cef, Gam );
	GetAlpha( P, Zrs, Gam, Alp );

	std::cout << P << ", ";
	for(auto itr = Alp.begin(); itr < Alp.end(); itr++ )
	{
		std::cout << *itr << ", ";
	}
	std::cout << std::endl;

	int i = 0;
	std::complex<double> base = Zrs[i];
	double alp = Alp[i];
	for (int j=0; j<N; j++)
	{
		if (i != j)
		{
			ApCircle<double> tmp( alp, base, Zrs[j] );
			std::cout << j << ": " << tmp.center() << ", " << tmp.radius() << std::endl;
		}
	}

	return EXIT_SUCCESS;
}
