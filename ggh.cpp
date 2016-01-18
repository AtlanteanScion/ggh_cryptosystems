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


class GGH
{
private:
	long dimension_n;	// Lattice dimension n
	long sigma;	// Error vector generation parameter sigma
	long parameter_l;	// Key generation parameter l
	long parameter_k;	// Key generation parameter k
	Mat<ZZ> private_R;	// Private basis R (low dual-orthogonality-defect)
	Mat<ZZ> public_B;	// Public basis B (high dual-orthogonality-defect)
	Mat<RR> private_R_Inv;	// R Inverse
	Mat<RR> public_B_Inv;	// B Inverse
	Mat<RR> transform_T;	// Private transform T

public:
	GGH(long dimension_n) : dimension_n(dimension_n)
	{
		sigma = 3;
		parameter_l = 4;
		parameter_k = parameter_l * ceil(1 + sqrt(dimension_n));

		cout << "Generating private key...\n";
		GeneratePrivate();

		cout << "Generating public key...\n";
		GeneratePublic();

		cout << "Generating private inverse...\n";
		GeneratePrivateInverse();

		cout << "Generating public inverse...\n";
		GeneratePublicInverse();

		cout << "Generating private transform...\n";
		GeneratePrivateTransform();
	}

	void Encrypt(Vec<ZZ>& cyphertext_c, Vec<ZZ>& plaintext_v)
	{
		// Error vector
		Vec<ZZ> error_e;
		GenerateErrorVector(error_e);

		// Encrypt
		// c = B * v + e
		// c = R * U * v + e
		cout << "Encrypting...\n";
		mul(cyphertext_c, public_B, plaintext_v);
		add(cyphertext_c, cyphertext_c, error_e);
	}

	void Decrypt(Vec<ZZ>& plaintext_v, Vec<ZZ>& cyphertext_c)
	{
		// Decrypt
		// T = B^-1 * R
		// v = T * round(R^-1 * c)

		// v = B^-1 * R * round(R^-1 * (R * U * v + e))
		// v = B^-1 * R * round(U * v + R^-1 * e)
		// v = B^-1 * R * U * v
		// v = B^-1 * B * v
		// v = v
		cout << "Decrypting...\n";
		Vec<RR> temp;
		mul(temp, private_R_Inv, conv<Vec<RR>>(cyphertext_c));
		RoundVector(temp, temp);
		mul(temp, transform_T, temp);
		RoundVector(plaintext_v, temp);
	}

	void GenerateRandomPlaintext(Vec<ZZ>& plaintext_v)
	{
		stringstream ss;
		long rand;
		ss << "[";
		for(int i = 0; i < dimension_n; ++i)
		{
			rand = RandomBnd(2 * dimension_n + 1) - dimension_n;
			ss << rand << " ";
		}
		ss << "0]";
		ss >> plaintext_v;
		plaintext_v.SetLength(dimension_n);	// Pad or truncate as necessary.
	}

	long GetDimension()
	{
		return dimension_n;
	}

private:
	// Fill with random integers (-l,l).
	// Add k * I.
	void GeneratePrivate()
	{
		diag(private_R, dimension_n, conv<ZZ>(parameter_k));
		for(long i = 0; i < dimension_n; ++i)
		{
			for(long j = 0; j < dimension_n; ++j)
			{
				private_R[j][i] += RandomBnd(2 * parameter_l + 1) - parameter_l;
			}
		}
	}

	// Add/subtract random vectors to each vector.  Do >=2 passes.
	void GeneratePublic()
	{
		long rand;
		transpose(public_B, private_R);
		for(long pass = 0; pass < 2; ++pass)
		{
			for(long i = 0; i < dimension_n; ++i)
			{
				for(long mix = 0; mix < dimension_n; ++mix)
				{
					if(i == mix)
					{
						continue;
					}

					rand = RandomBnd(7);
					if(rand == 0)
					{
						public_B[i] += public_B[mix];
					}
					else if(rand == 1)
					{
						public_B[i] -= public_B[mix];
					}
				}
			}
		}

		// LLL-reduction
//		ZZ det2;
//		LLL(det2, public_B);	// exact-aritmetic (slow but sure)
//		LLL_FP(public_B);	// floating-point (fast but only approximate)
		G_LLL_FP(public_B);	// floating-point Givens rotations (somewhat slower than LLL_FP, but more stable)

		transpose(public_B, public_B);
	}

	void GeneratePrivateInverse()
	{
		inv(private_R_Inv, conv<Mat<RR>>(private_R));
	}

	void GeneratePublicInverse()
	{
		inv(public_B_Inv, conv<Mat<RR>>(public_B));
	}

	// T = B^-1 * R
	void GeneratePrivateTransform()
	{
		mul(transform_T, public_B_Inv, conv<Mat<RR>>(private_R));
	}

	void GenerateErrorVector(Vec<ZZ>& v)
	{
		v.SetLength(dimension_n);

		for(long i = 0; i < dimension_n; ++i)
		{
			v[i] = sigma * (RandomBnd(2) * 2 - 1);	// {-sigma,sigma}
		}
	}

	void RoundVector(Vec<ZZ>& vec_rounded, Vec<RR>& vec)
	{
		vec_rounded.SetLength(vec.length());
		for(long i = 0; i < vec.length(); ++i)
		{
			RoundToZZ(vec_rounded[i], vec[i]);
		}
	}

	void RoundVector(Vec<RR>& vec_rounded, Vec<RR>& vec)
	{
		vec_rounded.SetLength(vec.length());
		for(long i = 0; i < vec.length(); ++i)
		{
			round(vec_rounded[i], vec[i]);
		}
	}
};

int RunTest(GGH& ggh)
{
	Vec<ZZ> plaintext_v;
	ggh.GenerateRandomPlaintext(plaintext_v);

	Vec<ZZ> cyphertext_c;
	ggh.Encrypt(cyphertext_c, plaintext_v);

	Vec<ZZ> decrypted_v;
	ggh.Decrypt(decrypted_v, cyphertext_c);

	if(plaintext_v == decrypted_v)
	{
		cout << "\nSuccess\n";
		return 1;
	}
	else
	{
		cout << "\nplaintext_v = \n" << plaintext_v << "\n";
		cout << "\ndecrypted_v = \n" << decrypted_v << "\n";
		cout << "\nFail\n";
		return 0;
	}

	return 0;
}

int main()
{
	long pass = 0;
	long total = 0;

	RR::SetPrecision(300);	// Default: 150 bits;  Minimum: 53 bits;  Maximum: 268435455 bits (NTL_OVFBND-1)

	for(long i = 0; i < 100; ++i)
	{
		GGH ggh = GGH(80);
		for(long j = 0; j < 100; ++j)
		{
			++total;
			cout << "\n\nTest: " << total << "\n";
			pass += RunTest(ggh);
			cout << "Success rate = " << ((float) pass/total) << "\n";
		}
	}

	return 0;
}

