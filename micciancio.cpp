#include <time.h>
#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/RR.h>
#include <NTL/vec_RR.h>
#include <NTL/mat_RR.h>
#include <NTL/HNF.h>

using namespace std;
using namespace NTL;

class GGH
{
private:
	long dimension_n;	// Lattice dimension n
	long parameter_sigma;	// Error vector generation parameter sigma
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
		long prec = conv<long>(CeilToZZ(3.33 * pow(conv<RR>(dimension_n), conv<RR>(1.2))));
		RR::SetPrecision(prec);

		parameter_sigma = 3;
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

	void Encrypt(Vec<ZZ>& cyphertext_c, Vec<ZZ>& plaintext_r)
	{
		// Encrypt
		// c = r mod B
		// c = r - Bx
		// xi = floor((ri - SUMj>i(bij * xj)) / (bii))
		cout << "Encrypting...\n";

		Vec<ZZ> x;
		ComputeEncryptionReductionVector(x, plaintext_r, public_B);

		mul(cyphertext_c, public_B, x);
		sub(cyphertext_c, plaintext_r, cyphertext_c);
	}

	void Decrypt(Vec<ZZ>& plaintext_r, Vec<ZZ>& cyphertext_c)
	{
		// Decrypt
		// T = B^-1 * R
		// -x = T * round(R^-1 * c)
		// r = c - B * -x
		//    Proof:
		//    -x = T * round(R^-1 * c)
		//    -x = B^-1 * R * round(R^-1 * (r - R * U * x))
		//    -x = B^-1 * R * round(R^-1 * r - U * x)
		//    -x = B^-1 * R * U * -x
		//    -x = B^-1 * B * -x
		//    -x = -x

		cout << "Decrypting...\n";
		Vec<RR> temp;
		Vec<ZZ> negX;
		mul(temp, private_R_Inv, conv<Vec<RR>>(cyphertext_c));
		RoundVector(temp, temp);
		mul(temp, transform_T, temp);
		RoundVector(negX, temp);
		mul(plaintext_r, public_B, negX);
		sub(plaintext_r, cyphertext_c, plaintext_r);
	}

	void GenerateRandomErrorVectorPlaintext(Vec<ZZ>& r)
	{
		r.SetLength(dimension_n);

		for(long i = 0; i < dimension_n; ++i)
		{
			r[i] = parameter_sigma * (RandomBnd(2) * 2 - 1);	// {-sigma,sigma}
		}
	}

	long GetDimension()
	{
		return dimension_n;
	}

private:
	// Fill with random integers in range [-l,l].
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

	// Compute Hermite Normal Form of private key.
	void GeneratePublic()
	{
		transpose(public_B, private_R);

		// Hermite Normal Form
		ZZ D;
		determinant(D, public_B);
		HNF(public_B, public_B, D);

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

	// xi = floor((ri - SUMj>i(bij * xj)) / (bii))
	void ComputeEncryptionReductionVector(Vec<ZZ>& x, Vec<ZZ>& r, Mat<ZZ>& B)
	{
		x.SetLength(dimension_n);
		for(long i = dimension_n - 1; i >= 0; --i)
		{
			ZZ sum = ZZ(0);
			for(long j = dimension_n - 1; j > i; --j)
			{
				sum += B[i][j] * x[j];
			}
			x[i] = (r[i] - sum) / B[i][i];
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
	ggh.GenerateRandomErrorVectorPlaintext(plaintext_v);

	Vec<ZZ> cyphertext_c;
	ggh.Encrypt(cyphertext_c, plaintext_v);

	Vec<ZZ> decrypted_v;
	ggh.Decrypt(decrypted_v, cyphertext_c);

	if(plaintext_v == decrypted_v)
	{
		cout << "Success\n";
		return 1;
	}
	else
	{
		cout << "Fail\n";
		cout << "\nplaintext_v = \n" << plaintext_v << "\n";
		cout << "\ndecrypted_v = \n" << decrypted_v << "\n";
		return 0;
	}

	return 0;
}

int main()
{
	long pass = 0;
	long total = 0;

	SetSeed(conv<ZZ>(time(0)));

	for(long i = 0; i < 10; ++i)
	{
		GGH ggh = GGH(80);
		for(long j = 0; j < 10; ++j)
		{
			++total;
			cout << "\nTest: " << total << "\n";
			pass += RunTest(ggh);
		}
		cout << "\nSuccess rate = " << ((float) pass/total) << "\n\n";
	}

	return 0;
}

