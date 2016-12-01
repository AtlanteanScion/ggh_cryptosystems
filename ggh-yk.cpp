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
	long parameter_s;	// Error vector generation parameter s
	long parameter_alpha;	// Key generation parameter alpha
	RR parameter_gamma;	// Key generation parameter gamma = alpha * n ^ beta (beta = 1/2)
	ZZ parameter_h;	// Message encoding parameter h
	long parameter_k;	// Message encoding parameter k
	Mat<ZZ> private_R;	// Private basis R (low dual-orthogonality-defect)
	Mat<ZZ> public_B;	// Public basis B (high dual-orthogonality-defect)
	Mat<RR> private_R_Inv;	// R Inverse
	Mat<RR> public_B_Inv;	// B Inverse

public:
	GGH(long dimension_n) : dimension_n(dimension_n)
	{
		// Selecting h and k can be problematic at low dimensions.
		if(dimension_n < 38)
		{
			dimension_n = 38;
		}

		long prec = conv<long>(CeilToZZ(3.33 * pow(conv<RR>(dimension_n), conv<RR>(1.2))));
		RR::SetPrecision(prec);

		parameter_s = 3;
		parameter_alpha = 4;
		parameter_gamma = conv<RR>(parameter_alpha) * sqrt(conv<RR>(dimension_n));

		// Select values h and k that satisfy conditions (6), (7), & (8).
		RR temp = conv<RR>(dimension_n) - 6 * sqrt(conv<RR>(dimension_n));
		parameter_k = conv<long>(floor(sqrt(temp)));
		parameter_h = FloorToZZ((temp - 1) / conv<RR>(parameter_k));
		parameter_h = max(parameter_h, conv<ZZ>(parameter_s + 1));
		parameter_k = conv<long>((temp - 1) / conv<RR>(parameter_h));

		cout << "Generating private key...\n";
		GeneratePrivate();

		cout << "Generating public key...\n";
		GeneratePublic();

		cout << "Generating public inverse...\n";
		GeneratePublicInverse();
	}

	void Encode(Vec<ZZ>& ptext, Vec<uint8_t>& msg)
	{
		long remain_h = parameter_k;
		long remain_s = msg.length() - parameter_k;
		long remain_zero = dimension_n - msg.length();

		ptext.SetLength(dimension_n);

		for(long ptext_index = 0, msg_index = 0; ptext_index < dimension_n; ++ptext_index)
		{
			ZZ val;
			long rnd = RandomBnd(remain_h + remain_s + remain_zero);

			if(rnd < remain_h)
			{
				val = parameter_h;
				--remain_h;
			}
			else if(rnd < remain_h + remain_s)
			{
				val = RandomBnd(parameter_s) + 1;
				--remain_s;
			}
			else
			{
				val = 0;
				--remain_zero;
			}

			if(val != 0)
			{
				if(msg[msg_index++] == 0)
				{
					val *= -1;
				}
			}

			ptext[ptext_index] = val;
		}
	}

	void Decode(Vec<uint8_t>& msg, Vec<ZZ>& ptext)
	{
		msg.SetLength(0);

		for(long ptext_index = 0; ptext_index < ptext.length(); ++ptext_index)
		{
			if(ptext[ptext_index] < 0)
			{
				msg.append(0);
			}
			else if(ptext[ptext_index] > 0)
			{
				msg.append(1);
			}
		}
	}

	void Encrypt(Vec<ZZ>& cyphertext_c, Vec<ZZ>& plaintext_e)
	{
		// Encrypt
		// c = e mod B
		// c = e - Bx
		// xi = floor((ei - SUMj>i(bij * xj)) / (bii))
		cout << "Encrypting...\n";

		Vec<ZZ> x;
		ComputeEncryptionReductionVector(x, plaintext_e, public_B);

		mul(cyphertext_c, public_B, x);
		sub(cyphertext_c, plaintext_e, cyphertext_c);
	}

	void Decrypt(Vec<ZZ>& plaintext_e, Vec<ZZ>& cyphertext_c)
	{
		// Decrypt

		// _s = R^-1 * c - round(R^-1 * c)
		cout << "Decrypting...\n";
		Vec<RR> temp;
		Vec<RR> rInv_c;
		Vec<RR> vec_s;
		mul(rInv_c, private_R_Inv, conv<Vec<RR>>(cyphertext_c));
		RoundVector(temp, rInv_c);
		sub(vec_s, rInv_c, temp);

		// e' = R * _s
		Vec<RR> e_prime;
		mul(e_prime, conv<Mat<RR>>(private_R), vec_s);

			// Calculate a:
			Vec<RR> vec_a;
			vec_a.SetLength(e_prime.length());
			// if e'i < -h - k, then ai = 1
			// if e'i > h + k, then ai = -1
			// else ai = 0
			RR h_plus_k = conv<RR>(parameter_h) + conv<RR>(parameter_k);
			for(long i = 0; i < e_prime.length(); ++i)
			{
				if(e_prime[i] < -h_plus_k)
				{
					vec_a[i] = 1;
				}
				else if(e_prime[i] > h_plus_k)
				{
					vec_a[i] = -1;
				}
				else
				{
					vec_a[i] = 0;
				}
			}


		// e = R * (_s + a)
		add(temp, vec_s, vec_a);
		mul(temp, conv<Mat<RR>>(private_R), temp);
		RoundVector(plaintext_e, temp);
	}

	void GenerateRandomMessage(Vec<uint8_t>& msg)
	{
		msg.SetLength(RandomBnd(dimension_n) + 1);

		for(long i = 0; i < msg.length(); ++i)
		{
			msg[i] = RandomBnd(2);
		}
	}

	long GetDimension()
	{
		return dimension_n;
	}

private:
	bool PrivateKeyIsValid()
	{
		cout << "  Generating private inverse...\n";
		GeneratePrivateInverse();

		// Diagonal entries, q, of inverse:  abs(q) <= 1/gamma
		// Non-diagonal entries, q, of inverse:  abs(q) < 2/gamma^2
		RR max_diag = conv<RR>(1) / parameter_gamma;
		RR max_non_diag = conv<RR>(2) / power(parameter_gamma, 2L);
		for(long i = 0; i < dimension_n; ++i)
		{
			for(long j = 0; j < dimension_n; ++j)
			{
				if(i == j && abs(private_R_Inv[j][i]) > max_diag)
				{
					// Constraint not satisfied; new key must be generated.
					return false;
				}
				if(i != j && abs(private_R_Inv[j][i]) >= max_non_diag)
				{
					// Constraint not satisfied; new key must be generated.
					return false;
				}
			}
		}

		// All entries in inverse satisfy constraints.
		return true;
	}

	// Fill with random integers {-1,0,1}.
	// Add gamma * I.
	void GeneratePrivate()
	{
		for(unsigned int attempt = 1; ; ++attempt)
		{
			cout << "  Attempt #" << attempt << "\n";
			diag(private_R, dimension_n, CeilToZZ(parameter_gamma));
			for(long i = 0; i < dimension_n; ++i)
			{
				for(long j = 0; j < dimension_n; ++j)
				{
					private_R[j][i] += RandomBnd(3) - 1;
				}
			}

			// If key satisfies constraints, then return, else, try again.
			if(PrivateKeyIsValid())
			{
				break;
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
	Vec<uint8_t> msg_in;
	ggh.GenerateRandomMessage(msg_in);

	Vec<ZZ> plaintext_e;
	ggh.Encode(plaintext_e, msg_in);

	Vec<ZZ> cyphertext_c;
	ggh.Encrypt(cyphertext_c, plaintext_e);

	Vec<ZZ> decrypted_e;
	ggh.Decrypt(decrypted_e, cyphertext_c);

	Vec<uint8_t> msg_out;
	ggh.Decode(msg_out, decrypted_e);

	if(msg_in == msg_out)
	{
		cout << "Success\n";
		return 1;
	}
	else
	{
		cout << "Fail\n";
		cout << "\nmsg_in  = \n" << msg_in << "\n";
		cout << "\nmsg_out = \n" << msg_out << "\n";
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

