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
	ZZ parameter_sigma;	// Error vector generation parameter sigma
	long parameter_alpha;	// Key generation parameter alpha
	ZZ parameter_gamma;	// Key generation parameter gamma = alpha * n
	ZZ parameter_h;	// Message encoding parameter h
	long parameter_k;	// Message encoding parameter k
	Mat<ZZ> private_B;	// Private basis B (low dual-orthogonality-defect)
	Mat<ZZ> public_W;	// Public basis W (high dual-orthogonality-defect)
	Mat<RR> private_B_Inv;	// B Inverse
	Mat<RR> public_W_Inv;	// W Inverse

public:
	GGH(long dimension_n) : dimension_n(dimension_n)
	{
		if(dimension_n < 4)
		{
			dimension_n = 4;
		}

		long prec = conv<long>(CeilToZZ(3.33 * pow(conv<RR>(dimension_n), conv<RR>(1.2))));
		RR::SetPrecision(prec);

		parameter_alpha = 3;
		parameter_gamma = conv<ZZ>(parameter_alpha * dimension_n);
		parameter_h = parameter_gamma / 2 + 1;
		parameter_k = dimension_n / 5;
		parameter_sigma = 2 * ((39 * dimension_n) / 152);

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
		long remain_sigma = dimension_n - parameter_k;

		msg.SetLength(remain_sigma);
		ptext.SetLength(dimension_n);

		for(long ptext_index = 0, msg_index = 0; ptext_index < dimension_n; ++ptext_index)
		{
			ZZ val;
			long rnd = RandomBnd(remain_h + remain_sigma);

			if(rnd < remain_h)
			{
				val = parameter_h;
				--remain_h;
			}
			else if(rnd < remain_h + remain_sigma)
			{
				val = RandomBnd(parameter_sigma / 2) + 1;
				if(msg[msg_index++] != 0)
				{
					val += parameter_sigma / 2;
				}
				--remain_sigma;
			}

			ptext[ptext_index] = val;
		}
	}

	void Decode(Vec<uint8_t>& msg, Vec<ZZ>& ptext)
	{
		msg.SetLength(0);

		for(long ptext_index = 0; ptext_index < ptext.length(); ++ptext_index)
		{
			if(1 <= ptext[ptext_index] && ptext[ptext_index] <= parameter_sigma / 2)
			{
				msg.append(0);
			}
			else if(parameter_sigma / 2 + 1 <= ptext[ptext_index] && ptext[ptext_index] <= parameter_sigma)
			{
				msg.append(1);
			}
		}
	}

	void Encrypt(Vec<ZZ>& cyphertext_c, Vec<ZZ>& plaintext_r)
	{
		// Encrypt
		// c = r mod W
		// c = Wx + r
		// x = -floor(W^-1 * r)
		cout << "Encrypting...\n";

		Vec<ZZ> x;
		Vec<RR> temp;
		mul(temp, public_W_Inv, conv<Vec<RR>>(plaintext_r));
		FloorVector(x, temp);
		mul(x, -1L, x);

		mul(cyphertext_c, public_W, x);
		add(cyphertext_c, cyphertext_c, plaintext_r);
	}

	void Decrypt(Vec<ZZ>& plaintext_r, Vec<ZZ>& cyphertext_c)
	{
		// Decrypt

		// u = B^-1 * c - round(B^-1 * c)
		cout << "Decrypting...\n";
		Vec<RR> temp;
		Vec<RR> bInv_c;
		Vec<RR> vec_u;
		mul(bInv_c, private_B_Inv, conv<Vec<RR>>(cyphertext_c));
		RoundVector(temp, bInv_c);
		sub(vec_u, bInv_c, temp);

		// r' = Bu = r - Be
		Vec<RR> r_prime;
		mul(r_prime, conv<Mat<RR>>(private_B), vec_u);

			// Calculate e:
			Vec<RR> vec_e;
			vec_e.SetLength(r_prime.length());
			// if r'i < 0, then ei = 1
			// else ei = 0
			for(long i = 0; i < r_prime.length(); ++i)
			{
				if(r_prime[i] < 0)
				{
					vec_e[i] = 1;
				}
				else
				{
					vec_e[i] = 0;
				}
			}

		// r = r' + Be
		mul(temp, conv<Mat<RR>>(private_B), vec_e);
		add(temp, r_prime, temp);
		RoundVector(plaintext_r, temp);
	}

	void GenerateRandomMessage(Vec<uint8_t>& msg)
	{
		msg.SetLength(dimension_n - parameter_k);

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

		// (4) Non-diagonal entries, q, of inverse:  abs(q) < 2/gamma^2
		// (8) Diagonal entries, q, of inverse:  abs(q) <= 2/gamma
		// (11) Diagonal entries, q, of inverse:  abs(q) > 1/gamma
		RR max_non_diag = conv<RR>(2) / conv<RR>(power(parameter_gamma, 2L));
		RR max_diag = conv<RR>(2) / conv<RR>(parameter_gamma);
		RR min_diag = conv<RR>(1) / conv<RR>(parameter_gamma);
		for(long i = 0; i < dimension_n; ++i)
		{
			for(long j = 0; j < dimension_n; ++j)
			{
				if(i == j)
				{
					if(abs(private_B_Inv[j][i]) > max_diag)
					{
						// Constraint not satisfied; new key must be generated.
						return false;
					}
					if(abs(private_B_Inv[j][i]) <= min_diag)
					{
						// Constraint not satisfied; new key must be generated.
						return false;
					}
				}
				else if(abs(private_B_Inv[j][i]) >= max_non_diag)
				{
					// Constraint not satisfied; new key must be generated.
					return false;
				}
			}
		}

		// All entries in inverse satisfy constraints.
		return true;
	}

	// Fill with random integers {-1,0}.
	// Add gamma * I.
	void GeneratePrivate()
	{
		for(unsigned int attempt = 1; ; ++attempt)
		{
			cout << "  Attempt #" << attempt << "\n";
			diag(private_B, dimension_n, parameter_gamma);
			for(long i = 0; i < dimension_n; ++i)
			{
				for(long j = 0; j < dimension_n; ++j)
				{
					private_B[j][i] += RandomBnd(2) - 1;
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
		transpose(public_W, private_B);

		// Hermite Normal Form
		ZZ D;
		determinant(D, public_W);
		HNF(public_W, public_W, D);

		transpose(public_W, public_W);
	}

	void GeneratePrivateInverse()
	{
		inv(private_B_Inv, conv<Mat<RR>>(private_B));
	}

	void GeneratePublicInverse()
	{
		inv(public_W_Inv, conv<Mat<RR>>(public_W));
	}

	void FloorVector(Vec<ZZ>& vec_floored, Vec<RR>& vec)
	{
		vec_floored.SetLength(vec.length());
		for(long i = 0; i < vec.length(); ++i)
		{
			FloorToZZ(vec_floored[i], vec[i]);
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

	Vec<ZZ> plaintext_r;
	ggh.Encode(plaintext_r, msg_in);

	Vec<ZZ> cyphertext_c;
	ggh.Encrypt(cyphertext_c, plaintext_r);

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

