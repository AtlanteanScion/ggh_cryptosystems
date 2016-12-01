#include <time.h>
#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
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
	Mat<ZZ> private_B_Inv;	// B Inverse
	Mat<ZZ> public_W_Inv;	// W Inverse
	ZZ private_B_Det;	// B Determinant
	ZZ public_W_Det;	// W Determinant

public:
	GGH(long dimension_n) : dimension_n(dimension_n)
	{
		if(dimension_n < 4)
		{
			dimension_n = 4;
		}

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
		mul(x, public_W_Inv, plaintext_r);
		DivideFloorVector(x, x, public_W_Det);
		mul(x, -1L, x);

		mul(cyphertext_c, public_W, x);
		add(cyphertext_c, cyphertext_c, plaintext_r);
	}

	void Decrypt(Vec<ZZ>& plaintext_r, Vec<ZZ>& cyphertext_c)
	{
		// Decrypt

		// u = B^-1 * c - round(B^-1 * c)
		cout << "Decrypting...\n";
		Vec<ZZ> temp;
		Vec<ZZ> bInv_c;
		Vec<ZZ> vec_u;
		mul(bInv_c, private_B_Inv, cyphertext_c);
		ModVector(bInv_c, bInv_c, private_B_Det);
		DivideRoundVector(temp, bInv_c, private_B_Det);
		mul(temp, temp, private_B_Det);
		sub(vec_u, bInv_c, temp);

		// r' = Bu = r - Be
		Vec<ZZ> r_prime;
		mul(r_prime, private_B, vec_u);
		DivideRoundVector(r_prime, r_prime, private_B_Det);

			// Calculate e:
			Vec<ZZ> vec_e;
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
		mul(temp, private_B, vec_e);
		add(plaintext_r, r_prime, temp);
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
		for(long i = 0; i < dimension_n; ++i)
		{
			for(long j = 0; j < dimension_n; ++j)
			{
				if(i == j)
				{
					if(abs(private_B_Inv[j][i]) * parameter_gamma > 2 * abs(private_B_Det))
					{
						// Constraint not satisfied; new key must be generated.
						return false;
					}
					if(abs(private_B_Inv[j][i]) * parameter_gamma <= 1 * abs(private_B_Det))
					{
						// Constraint not satisfied; new key must be generated.
						return false;
					}
				}
				else if(abs(private_B_Inv[j][i]) * parameter_gamma * parameter_gamma >= 2 * abs(private_B_Det))
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
		inv(private_B_Det, private_B_Inv, private_B);
	}

	void GeneratePublicInverse()
	{
		inv(public_W_Det, public_W_Inv, public_W);
	}

	void DivideFloorVector(Vec<ZZ>& vector_quotient, Vec<ZZ>& vector_dividend, ZZ& divisor)
	{
		vector_quotient.SetLength(vector_dividend.length());
		for(long i = 0; i < vector_dividend.length(); ++i)
		{
			div(vector_quotient[i], vector_dividend[i], divisor);
		}
	}

	void DivideRoundVector(Vec<ZZ>& quotient, Vec<ZZ>& dividend, ZZ& divisor)
	{
		quotient.SetLength(dividend.length());
		for(long i = 0; i < dividend.length(); ++i)
		{
			ZZ remainder;
			DivRem(quotient[i], remainder, dividend[i], divisor);
			if(remainder > divisor - remainder)
			{
				++quotient[i];
			}
		}
	}

	void ModVector(Vec<ZZ>& remainder, Vec<ZZ>& dividend, ZZ& divisor)
	{
		remainder.SetLength(dividend.length());
		for(long i = 0; i < dividend.length(); ++i)
		{
			rem(remainder[i], dividend[i], divisor);
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

