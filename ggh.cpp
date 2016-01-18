#include <iostream>
#include <sstream>
#include <string>
#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/vec_RR.h>
#include <NTL/mat_RR.h>
#include <NTL/LLL.h>

using namespace std;
using namespace NTL;

long genParam_l;
RR genParam_k;
long sigma;

//L*U : Lower/Upper triangular, {-1,1} on diag, {-1,0,1} in triangle
void getUnimodular(Mat<RR>& uni, const long dim_n)
{
	Mat<RR> upper = ident_mat_RR(dim_n);
	Mat<RR> lower = ident_mat_RR(dim_n);
	for(long i = 0; i < dim_n; ++i)
	{
		for(long j = 0; j < dim_n; ++j)
		{
			if(i > j)	// upper
			{
				upper[j][i] = RandomBnd(3) - 1;	// {-1,0,1}
				lower[j][i] = 0;
			}
			else if(i < j)	// lower
			{
				upper[j][i] = 0;
				lower[j][i] = RandomBnd(3) - 1;	// {-1,0,1}
			}
			else if(i == j)	// diag
			{
				upper[j][i] = RandomBnd(2) * 2 - 1;	// {-1,1}
				lower[j][i] = RandomBnd(2) * 2 - 1;	// {-1,1}
			}
		}
	}
	mul(uni, lower, upper);
}

//fill with random integers (-l,l)
void generatePrivateRandom(Mat<RR>& priv, const long dim_n)
{
	priv = ident_mat_RR(dim_n);
	for(long i = 0; i < dim_n; ++i)
	{
		for(long j = 0; j < dim_n; ++j)
		{
			priv[j][i] = RandomBnd(2 * genParam_l + 1) - genParam_l;
		}
	}
}

//fill with random integers (-l,l)
//add k * I
void generatePrivateRectangular(Mat<RR>& priv, const long dim_n)
{
	priv = diag(dim_n, genParam_k);
	for(long i = 0; i < dim_n; ++i)
	{
		for(long j = 0; j < dim_n; ++j)
		{
			priv[j][i] += RandomBnd(2 * genParam_l + 1) - genParam_l;
		}
	}
}

// Add/subtract random vectors to each vector.  Do >=2 passes.
void generatePublicMixing(Mat<RR>& pub, const Mat<RR>& priv, const long dim_n)
{
	long rand;
	pub = priv;
	for(long pass = 0; pass < 2; ++pass)
	{
		for(long i = 0; i < dim_n; ++i)
		{
			for(long mix = 0; mix < dim_n; ++mix)
			{
				if(i == mix)
				{
					continue;
				}
				
				rand = RandomBnd(7);
				if(rand == 0)
				{
					for(long j = 0; j < dim_n; ++j)
					{
						pub[j][i] += pub[j][mix];
					}
				}
				else if(rand == 1)
				{
					for(long j = 0; j < dim_n; ++j)
					{
						pub[j][i] -= pub[j][mix];
					}
				}
			}
		}
	}
}

// Multiply private key by at least 4 random unimodular matrices.
void generatePublicTransforms(Mat<RR>& pub, const Mat<RR>& priv, const long dim_n)
{
	Mat<RR> transform;
	pub = priv;
	for(long i = 0; i < 4; ++i)
	{
		getUnimodular(transform, dim_n);
		mul(pub, pub, transform);
	}
}

void errVec(Vec<RR>& v, const long dim_n)
{
	RR r;
	v.SetLength(dim_n);

	for(long i = 0; i < dim_n; ++i)
	{
		v[i] = sigma * (RandomBnd(2) * 2 - 1);	// {-sigma,sigma}
	}
}

void roundVec(Vec<RR>& v, const long dim_n)
{
	for(long i = 0; i < dim_n; ++i)
	{
		round(v[i], v[i]);
	}
}

int runTest()
{
	Mat<RR> priv_R;	// Private key (low dual-orthogonality-defect)
	Mat<RR> pub_B;	// Public key (high dual-orthogonality-defect)
	
	const long dim_n = 17;
	genParam_l = 4;
	genParam_k = genParam_l * ceil(1 + sqrt(dim_n));
	sigma = 3;
	
	cout << "Generating private key...\n";
//	generatePrivateRandom(priv_R, dim_n);	// BUG 0% success rate
	generatePrivateRectangular(priv_R, dim_n);	// Preferred by GGH
//	cout << "priv_R=\n" << priv_R << "\n";

	cout << "Generating public key...\n";
	generatePublicMixing(pub_B, priv_R, dim_n);	// Preferred by GGH
//	generatePublicTransforms(pub_B, priv_R, dim_n);
//	cout << "pub_B=\n" << pub_B << "\n";

	cout << "Generating private inverse...\n";
	Mat<RR> priv_R_Inv;
	inv(priv_R_Inv, priv_R);

	cout << "Generating public inverse...\n";
	Mat<RR> pub_B_Inv;
	inv(pub_B_Inv, pub_B);

//	cout << "priv_R_Inv=\n" << priv_R_Inv << "\n";
//	cout << "pub_B_Inv=\n" << pub_B_Inv << "\n";

	// Error vector
	Vec<RR> err_e;
	errVec(err_e, dim_n);
//	cout << "err_e=" << err_e << "\n";

	// Get message
	Vec<RR> msg_v;
//	cin >> msg_v;
	stringstream ss;
	long rand;
	ss << "[";
	for(int i = 0; i < dim_n; ++i)
	{
		rand = RandomBnd(2 * dim_n + 1) - dim_n;
		ss << rand << " ";
	}
	ss << "0]";
	ss >> msg_v;
	msg_v.SetLength(dim_n);	// Pad or truncate as necessary.
//	cout << "msg_v=" << msg_v << "\n";

	// Encrypt
	// c = B * v + e
	// c = R * U * v + e
	cout << "Encrypting...\n";
	Vec<RR> crypt_c;
	mul(crypt_c, pub_B, msg_v);
	add(crypt_c, crypt_c, err_e);
//	cout << "crypt_c=" << crypt_c << "\n";

	// Decrypt
	// T = B^-1 * R
	// v = T * round(R^-1 * c)

	// v = B^-1 * R * round(R^-1 * (R * U * v + e))
	// v = B^-1 * R * round(U * v + R^-1 * e)
	// v = B^-1 * R * U * v
	// v = B^-1 * B * v
	// v = v
	cout << "Decrypting...\n";
	Mat<RR> trans_T;
	Vec<RR> dec_v;
	mul(trans_T, pub_B_Inv, priv_R);
	mul(dec_v, priv_R_Inv, crypt_c);
	roundVec(dec_v, dim_n);
	mul(dec_v, trans_T, dec_v);
	roundVec(dec_v, dim_n);
//	cout << "dec_v=" << dec_v << "\n";

	Vec<RR> diff;
	diff = dec_v - msg_v;
	roundVec(diff, dim_n);
//	cout << "diff=" << diff << "\n";
	if(IsZero(diff))
	{
		cout << "\n\nSuccess\n";
		return 1;
	}
	else
	{
		cout << "\n\nFail\n";
		return 0;
	}

	return 0;
}

int main()
{
	long pass = 0;
	long total = 0;

	RR::SetPrecision(300);
	RR::SetOutputPrecision(300);

	for(long i = 0; i < 100; ++i)
	{
		++total;
		cout << "\n\nTest: " << total << "\n";
		pass += runTest();
		cout << "Success rate = " << ((float) pass/total) << "\n";
	}

	return 0;
}

