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
//	cout << "upper=\n" << upper << "\n\n";
//	cout << "lower=\n" << lower << "\n\n";
	}
	mul(uni, lower, upper);

//	cout << "determinant=" << determinant(upper) << "\n";
//	cout << "upper=\n" << upper << "\n\n";

//	cout << "determinant=" << determinant(lower) << "\n";
//	cout << "lower=\n" << lower << "\n\n";

//	cout << "determinant=" << determinant(uni) << "\n";
//	cout << "uni=\n" << uni << "\n\n";
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

	Mat<RR> fullTransform;
	Mat<RR> fullTransformInv;
// Multiply private key by at least 4 random unimodular matrices.
void generatePublicTransforms(Mat<RR>& pub, const Mat<RR>& priv, const long dim_n)
{
	Mat<RR> transform;
	fullTransform = ident_mat_RR(dim_n);
	Mat<RR> recPriv;
	pub = priv;
	for(long i = 0; i < 4; ++i)
	{
		getUnimodular(transform, dim_n);
//	cout << "pub=\n" << pub << "\n";
//	cout << "transform=\n" << transform << "\n\n";
		mul(pub, pub, transform);
		mul(fullTransform, fullTransform, transform);
	}
//	cout << "pub=\n" << pub << "\n";
//	cout << "fullTransform=\n" << fullTransform << "\n";
	inv(fullTransformInv, fullTransform);
//	cout << "fullTransformInv=\n" << fullTransformInv << "\n";
//	mul(recPriv, pub, fullTransformInv);
//	cout << "recPriv=\n" << recPriv << "\n";
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

RR maxVal(Mat<RR>& mtx, const long dim_n)
{
	RR maxElt;
	maxElt = 0;
	for(long i = 0; i < dim_n; ++i)
	{
		for(long j = 0; j < dim_n; ++j)
		{
			maxElt = max(maxElt, abs(mtx[j][i]));
		}
	}
	return maxElt;
}

void matRRtoZZ(Mat<ZZ>& matZZ, Mat<RR>& matRR, const long dim_n)
{
	matZZ = ident_mat_ZZ(dim_n);
	for(long i = 0; i < dim_n; ++i)
	{
		for(long j = 0; j < dim_n; ++j)
		{
			matZZ[j][i] = RoundToZZ(matRR[j][i]);
		}
	}
}

int runTest()
{
	Mat<RR> priv_R;	// Private key (low dual-orthogonality-defect)
	Mat<RR> pub_B;	// Public key (high dual-orthogonality-defect)
	
	const long dim_n = 4;
	genParam_l = 4;
//	genParam_k = round(genParam_l * sqrt(dim_n));
	genParam_k = genParam_l * ceil(1 + sqrt(dim_n));
	sigma = 3;
	
	cout << "Generating private key...\n";
//	generatePrivateRandom(priv_R, dim_n);	// BUG
	generatePrivateRectangular(priv_R, dim_n);
//	cout << "priv_R=\n" << priv_R << "\n";

	cout << "Generating public key...\n";
	if(1)//RandomBnd(2))
	{
		cout << "Mixing...\n";
		generatePublicMixing(pub_B, priv_R, dim_n);
	}
	else
	{
		cout << "Transform...\n";
		generatePublicTransforms(pub_B, priv_R, dim_n);	// BUG
	}
//cout << "max pub_B=\n" << maxVal(pub_B, dim_n) << "\n";
//	cout << "pub_B=\n" << pub_B << "\n";
ZZ det2;
mat_ZZ inB;
mat_ZZ outU;
mat_ZZ newB;
matRRtoZZ(inB, pub_B, dim_n);
det2 = determinant(inB);
cout << "determinant=" << det2 << "\n";
det2 *= det2;
cout << "inB=\n" << inB << "\n";
LLL(det2, inB, outU);
cout << "outU=\n" << outU << "\n";
mul(newB, outU, inB);
cout << "newB=\n" << newB << "\n";
return 0;
	cout << "Generating private inverse...\n";
	Mat<RR> priv_R_Inv;
	inv(priv_R_Inv, priv_R);
	cout << "max priv_R_Inv=\n" << maxVal(priv_R_Inv, dim_n) << "\n";

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

// Alt decrypt
// T == U^-1
// U^-1 * round(R^-1 * c)
//cout << "trans_T=" << trans_T << "\n";
//cout << "fullTransformInv=" << fullTransformInv << "\n";
//cout << "Alt Decrypting...\n";
//Vec<RR> dec_v2;
//mul(dec_v2, priv_R_Inv, crypt_c);
//roundVec(dec_v2, dim_n);
//mul(dec_v2, fullTransformInv, dec_v2);
//roundVec(dec_v2, dim_n);
//cout << "dec_v2=" << dec_v2 << "\n";

//cout << "fullTransform=\n" << fullTransform << "\n";
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
//cout << "priv_R=\n" << priv_R << "\n";
//cout << "pub_B=\n" << pub_B << "\n";
//cout << "priv_R_Inv=\n" << priv_R_Inv << "\n";
//cout << "pub_B_Inv=\n" << pub_B_Inv << "\n";
//cout << "fullTransform=\n" << fullTransform << "\n";
//cout << "fullTransformInv=" << fullTransformInv << "\n";
//cout << "trans_T=" << trans_T << "\n";
//cout << "err_e=" << err_e << "\n";
//cout << "msg_v=" << msg_v << "\n";
//cout << "crypt_c=" << crypt_c << "\n";
//cout << "dec_v=" << dec_v << "\n";
Vec<RR> sigmaCheck;
mul(sigmaCheck, priv_R_Inv, err_e);
cout << "sigmaCheck=" << sigmaCheck << "\n";
roundVec(sigmaCheck, dim_n);
cout << "sigmaCheck=" << sigmaCheck << "\n";
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
cout << "precision() = " << RR::precision() << "\n";

	for(long i = 0; i < 100; ++i)
	{
		++total;
		cout << "\n\nTest: " << total << "\n";
		pass += runTest();
		cout << "Success rate = " << ((float) pass/total) << "\n";
if(pass!=total)return 1;		
	}

	return 0;
}

