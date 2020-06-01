// Xiaoqi duan <duanxq17@mails.tsinghua.edu.cn>

#ifndef COMPARISION_H_
#define COMPARISION_H_

#include <libscapi/include/infra/Measurement.hpp>
#include "ProtocolParty.h"
#include "Interpolate.h"
#include <vector>
#include<ctime>
#include<thread>
#include <algorithm>    // std::min
//append B to A
#define _append(A,B) A.insert(A.end(),B.begin(),B.end())

//#define in_range(x) true
//#define in_range(x) (((x) < field->GetElement(1ll<<63)) || ((-x) < field->GetElement(1ll<<63)))
template<class FieldType>
class CompareGate 
{
	private:
		ProtocolParty<FieldType> *helper;
		Interpolate<FieldType> interp;
		int eleSize,blank_size;
		vector<vector<FieldType> > OrVector; //the share of coefficients of fan-in Or function
		vector<FieldType> _bitSharesValue;
		vector<vector<FieldType> > _bitSharesBits;
		//vector<FieldType> _zeroShares;
                vector<FieldType> chkCVec;
                vector<vector<FieldType>> chkAVec, chkBVec; // for inner products
		//int _zeroShareOffset = 0;
		int _bitShareOffset = 0;
		//for safety, we use 40 bit floats
		//fixed-point floats, 32 decimals
		//can it work?
		int _k = 62;
		int _m = 48;
		int _kappa = 0;
		int iteration;
		int m_partyID;
		int numThreads = 16;
		TemplateField<FieldType> *field; // used to call some field functions
		vector<vector<FieldType> > _Ai;
		vector<FieldType> _bi;
		int times; //times of rerunning
		int n_iter; //times of iteration
		string inputsFile;
		int TruncPR_t = 0;
		int Compare_t = 0; 
		int Mult_t = 0;
		int Veri_t = 0;
	public:
		CompareGate(ProtocolParty<FieldType> *ptr,int siz,int blk_siz,int m_partyID,TemplateField<FieldType> *field,int n,string inputf);
		~CompareGate();
		long transElement(FieldType tmp);
  void run();
  void inputPhase();
  void runOnline();
  void runOffline();
  void verificationPhase();
  void writeOutput(vector<FieldType> &res,string f);



  
		// Compgate functions
		void getRandomBitShare(int num,vector<FieldType> &res, vector<vector<FieldType> > &bits);
		//void genRandom01(int num,vector<vector<FieldType>> &bin,vector<FieldType> &val);
		void compRandom(vector<FieldType> &a,vector<FieldType> &b,vector<FieldType> &res);
		//void compRandom(vector<FieldType> &a, vector<FieldType> &b, vector<FieldType> &dest);
		// given known element c and binary share a, store (c<a) to dest
		void compGiven(vector<vector<FieldType> > &a, vector<FieldType> &c, vector<FieldType> &dest);
		void compGivenNoPrefixOr(vector<vector<FieldType> > &a, vector<FieldType> &c, vector<FieldType> &dest);
		void computeAnswer(vector<vector<FieldType> > &b, vector<vector<FieldType> > &c,vector<FieldType> &res);
		void compP(vector<vector<FieldType> > &a, vector<FieldType> &dest);
		// compare the function (a < p/2)
		void compHalfP(vector<FieldType> &a, vector<FieldType> &dest);
		// compute square inverse of a given element
		void compSqrInverse(vector<FieldType> &a, vector<FieldType> &divs);
		//notice: only supported for Mersenne128
		void compSqrInverse_SingleThread(vector<FieldType> &a, vector<FieldType> &divs, int threadId);
		//FieldType compSqrInverse(FieldType &a);
		//save the prefix or of array a to array res
		void compPrefixOr(vector<vector<FieldType> > &a, vector<vector<FieldType> > &res);
		void compFanInOr(int num, vector<vector<FieldType> > &a, vector<FieldType> &res);
		int generateBitShares(int numOfBitShares);
		void getRandomInvPairs(int num, vector<FieldType> &b, vector<FieldType> &invb, vector<FieldType> &tmpc);
		void getBitDecomp(FieldType &e, vector<FieldType> &dest);
		///////////////////////////////////////
		//only used for multiplication prod. Par: m and k
		void TruncPR(vector<FieldType> &a,vector<FieldType> &res);//vector<int> &k,vector<int> &m);
		void TruncPRSecure(vector<FieldType> &a,vector<FieldType> &res);
		void doubleVecMult(vector<FieldType> &a,vector<FieldType> &b,vector<FieldType> &res, int groupSize = 1);
		void SoftThres(vector<FieldType> &thres, vector<FieldType> &a, vector<FieldType> &res);

		//compute 1/a[i] under double
		void doubleInverse(FieldType a,FieldType &res);
		void getBits(FieldType a,vector<FieldType> &bits);
		//return a mod r. Require: r < 2^64
		long getMod(FieldType a,long r);
		FieldType getDiv(FieldType a,long r);
		// TODO:
		void readLassoInputs();
		//void readLassoInput(vector<vector<FieldType> > & Ai, vector<FieldType> &bi);
		void runLasso(int iter,FieldType lambda, FieldType rho, vector<vector<FieldType> > & Ai, vector<FieldType> &bi, vector<FieldType> &res);
		//this is for debug
		bool in_range(FieldType tmp);
};

/////////////////////////////////////////////
//START OF COMPAREGATE
template<class FieldType>
bool CompareGate<FieldType>::in_range(FieldType tmp)
{
	cout<<"Shouldn't be here!"<<endl;
	abort();
	return true;
}

template<>
inline bool CompareGate<ZpMersenne127Element>::in_range(ZpMersenne127Element tmp)
{
	if(tmp.elem < (1ull<<63)) return true;
	tmp = tmp + (1ull<<63);
	if(tmp.elem < (1ull<<63)) return true;
	return false;
}

template<class FieldType>
FieldType CompareGate<FieldType>::getDiv(FieldType a,long r)
{
	long tmp = transElement(a);
	return field->GetElement(tmp / r);
}

template<>
inline ZpMersenne127Element CompareGate<ZpMersenne127Element>::getDiv(ZpMersenne127Element a, long r)
{
	/*if(a.elem > (1ull<<63))
	{
		cout<<"Invalid getDiv!"<<endl;
		abort();
	}*/
	__int128_t elem = a.elem;
	elem = elem / r;
	return ZpMersenne127Element(elem);
}

template<class FieldType>
long CompareGate<FieldType>::getMod(FieldType a,long r)
{
	long tmp = transElement(a);
	return tmp % r;
}

template<>
inline long CompareGate<ZpMersenne127Element>::getMod(ZpMersenne127Element a,long r)
{
	/*if(a.elem > (1ull<<63))
	{
		cout<<"Invalid getMod!"<<endl;
		abort();
	}*/
	__int128_t tmp = a.elem;
	return (long)(tmp % r);
}

template<class FieldType>
void CompareGate<FieldType>::getBits(FieldType a,vector<FieldType> &bits)
{
	long tmp = transElement(a);
	bits.resize(eleSize);
	for(int i=eleSize-1; i>=0; i--)
	{
		bits[i] = field->GetElement(tmp % 2);
		tmp /= 2;
	}
}

template<>
inline void CompareGate<ZpMersenne127Element>::getBits(ZpMersenne127Element a, vector<ZpMersenne127Element> &bits)
{
	__uint128_t tmp = a.elem;
	bits.resize(eleSize);
	for(int i=eleSize-1; i>=0; i--)
	{
		bits[i] = field->GetElement(tmp % 2);
		tmp /=2;
	}
}

	template<>
inline long CompareGate<ZpMersenneLongElement>::transElement(ZpMersenneLongElement tmp)
{
	return tmp.elem;
}

	template<>
inline long CompareGate<ZpMersenneIntElement>::transElement(ZpMersenneIntElement tmp)
{
	return tmp.elem;
}

	template<>
inline long CompareGate<ZZ_p>::transElement(ZZ_p tmp)
{
	return conv<long>(tmp);
}

	template<>
inline long CompareGate<ZpMersenne127Element>::transElement(ZpMersenne127Element tmp)
{
	return (long)tmp.elem;
}

	template<class FieldType>
CompareGate<FieldType>::CompareGate(ProtocolParty<FieldType> *ptr,int siz,int bsiz,int id,TemplateField<FieldType> *f,int n,string inp)
{
	if(flag_print)
		cout<<"testing FieldType():"<<FieldType(1)<<endl;
	m_partyID = id;
	helper = ptr;
	_bitShareOffset=0;
	eleSize = siz;
	blank_size = bsiz;
	field = f;
	//N = n;
	//T = (N-1)/2;
	times = 1;
	n_iter = 10;
	inputsFile = inp;	
}

	template<class FieldType>
CompareGate<FieldType>::~CompareGate()
{
}
template<class FieldType>
void CompareGate<FieldType>::compSqrInverse_SingleThread(vector<FieldType> &a, vector<FieldType> &divs, int threadId)
{
}

template<>
inline void CompareGate<ZpMersenne127Element>::compSqrInverse_SingleThread(vector<ZpMersenne127Element> &a, vector<ZpMersenne127Element> & divs,int threadId)
{
	int tot = a.size();
	for(int i = threadId; i<tot; i+=numThreads)
	{
		ZpMersenne127Element tmp(a[i]),tmp2;
		for(int j=0; j<125; j++)
			tmp = tmp * tmp;
		tmp2 = tmp + tmp;
		if(tmp2.elem < tmp.elem) tmp = ZpMersenne127Element(0) - tmp;
		if(tmp!=ZpMersenne127Element(0)) divs[i] = ZpMersenne127Element(1) / tmp;
		else divs[i] = 0;
	}
}
//warning: this may not be pure random. Please Use Specified version of ZZ_p and ZpMersenne127
	template <class FieldType>
void CompareGate<FieldType>::compSqrInverse(vector<FieldType> &num, vector<FieldType> &divs)
{
	//if(res.size()<num.size()) res.resize(num.size());
	if(divs.size()<num.size()) divs.resize(num.size());
	for(int i=0; i<num.size(); i++)
	{
		//res[i] = num[i].sqrt();
		FieldType tmp = num[i].sqrt();
		if(tmp!=FieldType(0)) divs[i] = FieldType(1) / tmp; //res[i];
		else divs[i] = 0;
		//if(res[i] + res[i] < res[i]) res[i] = FieldType(0) - res[i];
	}
}

//since we don't use ZZ_p,
//this is badly implemented, with very low efficiency.
	template<>
inline void CompareGate<ZZ_p>::compSqrInverse(vector<ZZ_p> &num, vector<ZZ_p> &divs)
{
	//the mod is 2147483647
	//we don't need to init ZZ (?)
	//if(flag_print)
	//	cout<<"Computing Sqr Inv of"<<num<<":"<<endl;;ZZ md = _ZZPRIME; //rep(_PRIME);
	ZZ md((1ll<<61)-1); //rep(_PRIME);
	//	cout<<"Prime:"<<md<<endl;
	ZZ tmp = (md+1)/4;
	//if(flag_print)
	//	cout<<"Times:"<<tmp<<endl;
	int tot=num.size();
	//if(res.size()<tot)
	//	res.resize(tot);
	if(divs.size()<tot)
		divs.resize(tot);
	for(int i=0; i<num.size(); i++)
	{
		ZZ_p t0 = power(num[i],tmp);
		if(conv<ZZ>(t0+t0) < conv<ZZ>(t0)) t0 = ZZ_p(0) - t0;
		if(t0!=0) divs[i] = 1 / t0;
		else divs[i] = 0;
	}
	
	//ZZ_p x = power(num,tmp); //qpow(num,tmp); //power(num,tmp);
	//if(flag_print)
	//	cout<<"Res:"<<x<<endl<<"Check:"<<x * num<<endl;
	//return x;
}

template<>
inline void CompareGate<ZpMersenne127Element>::compSqrInverse(vector<ZpMersenne127Element> &num, vector<ZpMersenne127Element> &divs)
{
	int tot = num.size();
	//if(res.size()<tot) res.resize(tot);
	if(divs.size()<tot) divs.resize(tot);
	//multithreading here
	vector<thread> threads;
	if(numThreads > 1)
	{
		threads.resize(numThreads-1);
		for(int t=0; t<numThreads - 1; t++)
			threads[t] = thread(&CompareGate::compSqrInverse_SingleThread, this, ref(num),  ref(divs), t+1);
	}
	//(md+1)/4: **2 for 125 times
	
	for(int i=0; i<tot; i+=numThreads)
	{
		ZpMersenne127Element tmp(num[i]),tmp2;
		for(int j=0; j<125; j++)
			tmp = tmp * tmp;
		tmp2 = tmp + tmp;
		if(tmp2.elem < tmp.elem) tmp = ZpMersenne127Element(0) - tmp;
		//res[i]=tmp;
		if(tmp!=ZpMersenne127Element(0)) divs[i] = ZpMersenne127Element(1) / tmp;
		else divs[i]=0;
	}
	if(numThreads > 1)
	{
		for(int t=0; t<numThreads - 1; t++)
			threads[t].join();
	}
}

//return value: succesfully generated shares
	template <class FieldType>
int CompareGate<FieldType>::generateBitShares(int num)
{
	_bitShareOffset = 0;
	_bitSharesValue.clear();
	auto _t1 = high_resolution_clock::now();
	int suc_cnt=0;
	_bitShareOffset = 0;
	//we pick out 2*required bits
	int tot = eleSize * num;
	if(flag_print)
		cout<<"Required #bitshare:"<<tot<<endl;
	vector<FieldType> tempShares,resShares,secrets;
	_bitSharesValue.resize(tot);
	helper->getPadShares(tot, tempShares);
	if(flag_print)
		cout<<"Generating random shares, #"<<tot<<endl;
	/*vector<FieldType> PlainText;
	if(flag_print)
	{
		helper->openShare(tot,tempShares,PlainText);
		cout<<"plaintext opened"<<endl;
	}*/
	//now do multiplications and rule out 0 shares
	helper->DNMultVec(tempShares,tempShares,resShares,1);
        for (int i=0; i<tot; i++) {
          chkAVec.push_back(vector<FieldType>(1, tempShares[i]));
          chkBVec.push_back(vector<FieldType>(1, tempShares[i]));
        }
	_append(chkCVec,resShares);
	if(flag_print)
		cout<<"DNMult finished"<<endl;

        verificationPhase();    // fix: verify before open mult results
	helper->openTShares(tot, true, resShares,secrets);
	auto _t2 = high_resolution_clock::now();
	if(flag_print)
	{
		cout<<"Share Opened in GenB(),"<<tot<<endl;
		cout<<secrets[0]<<","<<secrets[1]<<endl;
	}
	//for each opened, check if is 0, if not we generate the share
	compSqrInverse(secrets,resShares);
	auto _t3 = high_resolution_clock::now();
	if(flag_print)
		cout<<"Inverse computed"<<endl;
	FieldType inv2 = FieldType(1) / FieldType(2);
	for(int i=0; i<tot; i++)
		if(resShares[i] != field->GetElement(0)) //valid
		{
			/*if(flag_print)
			{
				if(i%1000000==0)
					cout<<"Generating bit "<<i<<endl;
			}*/
			//if(flag_print)
			//	cout<<"Trying to calc inv for "<<secrets[i]<<endl;
			//FieldType cur(invs[i]); // = compSqrInverse(secrets[i]);
			//if(flag_print)
			//	cout<<"find invs:"<<cur<<","<<PlainText[i]<<endl;
			//cout<<"Inv2:"<<FieldType(1) / FieldType(2)<<endl;
			//cur = tempShares[i] / cur * inv2 + inv2;
			//cur = tempShares[i] / cur + inv2; //FieldType(1))/ FieldType(2);
			//if(flag_print)
			//	cout<<"found share:"<<cur<<endl;
			/*if(flag_print)
			  {
			  vector<FieldType> tp(1),res;
			  tp[0] = cur;
			  helper->openShare(1,tp,res);
			  cout<<"Generated a "<<res[0]<<endl;
			  }*/
			_bitSharesValue[_bitShareOffset] = tempShares[i] * resShares[i] * inv2 + inv2;
			suc_cnt++;
			_bitShareOffset++;
			/*if(flag_print)
			{
				vector<FieldType> d0(1),d1;
				d0[0] = cur;
				helper->openShare(1,d0,d1);
				//cout<<"generate a bit "<<d1[0]<<endl;
				if(d1[0]!=FieldType(0) && d1[0]!=FieldType(1))
				{
					cout<<"invalid generated bit!"<<endl;
					abort();
				}
			}*/
			//if(flag_print)
			//	cout<<"Suc count:"<<suc_cnt<<endl;
		}
	//else if(flag_print)
	//	cout<<i<<" failed"<<endl;
	if(flag_print)
		cout<<"Bit Generate done"<<endl;
	auto _t4 = high_resolution_clock::now();
	cout<<"Mult time:"<<duration_cast<microseconds>(_t2-_t1).count()<<endl;
	cout<<"Inv time:"<<duration_cast<microseconds>(_t3-_t2).count()<<endl;
	cout<<"Mult time:"<<duration_cast<microseconds>(_t4-_t3).count()<<endl;
	_bitShareOffset = 0;
	return suc_cnt;
}


	template<class FieldType>
void CompareGate<FieldType>::compRandom(vector<FieldType> &a, vector<FieldType> &b, vector<FieldType> &res)
{
	auto _t0 = high_resolution_clock::now();
	//auto _t0 = clock();
	//here a,b,c means the *share* [a],[b],[c] that are hold by this party
	int cnt = a.size();
	if(flag_print)
		cout<<"Entering compRandom with size "<<cnt<<endl;
	//vector<FieldType> c,w,x,y,tmp,tmp2,tmp3;
	vector<FieldType> c,y;
	for(int l=0; l<cnt; l++)
	{
		//auto gate = gates[l];
		//a.push_back(helper->gateShareArr[gate.input1]);
		//b.push_back(helper->gateShareArr[gate.input2]);
		c.push_back(a[l] - b[l]);
	}
	//FieldType c = a - b;
	//same notation as the paper, share of [a<p/2],[b<p/2],[a-b<p/2]
	//FieldType w,x,y,tmp,tmp2,tmp3;
	//if(flag_print)
	//	cout<<"Entering comphalf"<<endl;
	/*if(flag_print)
	  {
	  cout<<"checking input:"<<endl;
	  vector<FieldType> d0(3),d1;
	  d0[0] = a; d0[1] = b; d0[2]=c;
	  helper->openShare(3,d0,d1);
	  cout<<d1[0]<<","<<d1[1]<<","<<d1[2]<<endl;
	  }*/
	// compHalfP(a,w);
	//	cout<<"Entering comphalf 2"<<endl;
	// compHalfP(b,x);
	compHalfP(c,y);
        // NOTE: since we always trunc values, will never get > p/2
        FieldType _one = field->GetElement(1);
        res.resize(cnt);
        for (int i=0; i<cnt; i++) {
          res[i] = _one - y[i];
        }
	/*if(flag_print)
	  {
	  cout<<"checking compHalfP results:"<<endl;
	  vector<FieldType> d0(3),d1;
	  d0[0]=w; d0[1]=x; d0[2]=y;
	  helper->openShare(3,d0,d1);
	  cout<<d1[0]<<","<<d1[1]<<","<<d1[2]<<endl;
	  }*/
	//    cout<<"Exiting compHalf"<<endl;
	//two innovations since we need to compute xy and w*(x+y-2xy)
	//first round: compute xy
	// helper->DNMultVec(x,y,tmp,1);
	// _append(chkA,x);
	// _append(chkB,y);
	// _append(chkC,tmp);
	// for(int l=0; l<cnt; l++)
	// 	tmp2.push_back(x[l] + y[l] - FieldType(2) * tmp[l]);
	// helper->DNMultVec(w,tmp2,tmp3,1);
	// _append(chkA,w);
	// _append(chkB,tmp2);
	// _append(chkC,tmp3);
	// if(res.size()<cnt)
	// 	res.resize(cnt);
	// for(int l=0; l<cnt; l++)
	// {
	// 	//auto gate = gates[l];
	// 	//helper->gateShareArr[gate.output] = tmp3[l] + FieldType(1) - y[l] - x[l] + tmp[l];
	// 	res[l] = tmp3[l] + FieldType(1) - y[l] - x[l] + tmp[l];
	// }
	auto _t1 = high_resolution_clock::now();
	Compare_t += duration_cast<microseconds>(_t1 - _t0).count();
}
// given known element c and binary share a, store (c<a) to dest
/*	template<class FieldType>
	void CompareGate<FieldType>::compGiven(vector<vector<FieldType> > &a, vector<FieldType> &c, vector<FieldType> &dest)
	{
	}
	template<>
	inline void CompareGate<ZZ_p>::compGiven(vector<vector<ZZ_p> > &a, vector<ZZ_p> &c, vector<ZZ_p> &dest)
 */
	template<class FieldType>
void CompareGate<FieldType>::compGiven(vector<vector<FieldType> > &a, vector<FieldType> &c, vector<FieldType> &dest)
{
	if(a.size()!=c.size())
	{
		cout<<"INVALID PARM for compGiven"<<endl;
		abort();
	}
	vector<FieldType> d0,d1;
	//vector<ZZ_p> d0,d1;
	int tot = a.size();
	if(flag_print)
	{
		cout<<"Entering compGiven"<<endl;
	}
	vector<vector<FieldType> > cBits;
	//vector<vector<ZZ_p> > cBits;
	cBits.resize(tot);
	for(int j=0; j<tot; j++)
	{
		getBits(c[j],cBits[j]);
		/*long tmp = transElement(c[j]);
		//ZZ tmp = rep(c[j]);
		cBits[j].resize(eleSize);
		for(int i=eleSize-1; i>=0; i--)
		{
			cBits[j][i] = field->GetElement(tmp % 2);
			tmp /=2;
		}*/
		/*if(flag_print)
		  {
		  for(int i=0; i<eleSize; i++)
		  cout<<cBits[i]<<",";
		  cout<<endl;
		  for(int i=0; i<eleSize; i++)
		  cout<<d1[i]<<",";
		  cout<<endl;
		  }*/
	}
	if(flag_print)
		cout<<"compGiven::xorSum"<<endl;
	//step 1: compute bitwise xor
	vector<vector<FieldType> > xorSum;
	//vector<vector<ZZ_p> > xorSum;
	xorSum.resize(tot);
	for(int j=0; j<tot; j++)
	{
		xorSum[j].resize(eleSize);
		for(int i=0; i<eleSize; i++)
		{
			xorSum[j][i] = a[j][i] + cBits[j][i];
			if(cBits[j][i]==1)
				xorSum[j][i] -= field->GetElement(2) * a[j][i];
		}
	}
	/*if(flag_print)
	  {
	  cout<<"Xor sum:"<<endl;
	  helper->openShare(eleSize,xorSum,d1);
	  for(int i=0; i<eleSize; i++)
	  cout<<d1[i]<<",";
	  cout<<endl;
	  }*/
	//step 2: compute prefix or
	vector<vector<FieldType> > prefixOr,E;
	//vector<vector<ZZ_p> > prefixOr,E;
	if(flag_print)
		cout<<"Entering PrefixOr"<<endl;
	compPrefixOr(xorSum,prefixOr);
	/*if(flag_print)
	  {
	  cout<<"Exiting compGiven::prefixor"<<endl;
	  helper->openShare(eleSize,prefixOr,d1);
	  for(int i=0; i<eleSize; i++)
	  cout<<d1[i]<<",";
	  cout<<endl;
	  }*/
	E.resize(tot);
	dest.resize(tot);
	for(int j=0; j<tot; j++)
	{
		E[j].resize(eleSize);
		E[j][0] = prefixOr[j][0];
		for(int i=1; i<eleSize; i++)
			E[j][i] = prefixOr[j][i] - prefixOr[j][i-1];
	}
	for(int j=0; j<tot; j++)
	{
		dest[j] = field->GetElement(0);
		for(int i=0; i<eleSize; i++)
			if(cBits[j][i]==0)
				dest[j] = dest[j] + E[j][i];
	}
	/*if(flag_print)
	  {
	  vector<ZZ_p> d0(1),d1;
	  d0[0]=dest;
	  helper->openShare(1,d0,d1);
	  if(d1[0]!=0 && d1[0]!=1)
	  {
	  cout<<"Invalid compGiven:"<<d1[0]<<endl;
	  abort();
	  }
	  cout<<"Exiting CompGiven"<<endl;
	  }*/
}
	template<class FieldType>
void CompareGate<FieldType>::compGivenNoPrefixOr(vector<vector<FieldType> > &x, vector<FieldType> &r, vector<FieldType> &dest)
	/*{
	  }
	  template<>
	  inline void CompareGate<ZZ_p>::compGivenNoPrefixOr(vector<vector<ZZ_p> > &x, vector<ZZ_p> &r, vector<ZZ_p> &dest)
	 */
{
	if(x.size()!=r.size())
	{
		cout<<"INVALID PARM for compGiven"<<endl;
		abort();
	}
	//vector<FieldType> d0,d1;
	//vector<ZZ_p> d0,d1;
	int tot = x.size();
	if(flag_print)
	{
		cout<<"Entering compGiven"<<endl;
	}
	vector<vector<FieldType> > rBits;
	//vector<vector<ZZ_p> > rBits;
	rBits.resize(tot);
	for(int j=0; j<tot; j++)
	{
		getBits(r[j],rBits[j]);
		/*
		long tmp = transElement(r[j]);
		rBits[j].resize(eleSize);
		for(int i=eleSize-1; i>=0; i--)
		{
			rBits[j][i] = field->GetElement(tmp % 2);
			tmp /=2;
		}*/
	}
	if(flag_print)
		cout<<"compGivenNew:step 1"<<endl;
	//step 1: compute b and c
	vector<vector<FieldType> > b,c;
	//vector<vector<ZZ_p> > b,c;
	b.resize(tot);
	c.resize(tot);
	for(int j=0; j<tot; j++)
	{
		b[j].resize(eleSize);
		c[j].resize(eleSize);
		for(int i=0; i<eleSize; i++)
		{
			b[j][i] = field->GetElement(1) - rBits[j][i] - x[j][i] + field->GetElement(2) * rBits[j][i] * x[j][i];
			if(rBits[j][i] == 0) c[j][i] = x[j][i];
			else c[j][i] = 0;
		}
	}
	//step 2: compute r
	computeAnswer(b,c,dest);
}
	template<class FieldType>
void CompareGate<FieldType>::computeAnswer(vector<vector<FieldType> > &b,vector<vector<FieldType> > &c, vector<FieldType> &res)
	//{
	//}
	//tempate<>
	//inline void CompareGate<ZZ_p>::computeAnswer(vector<vector<ZZ_p> > &b,vector<vector<ZZ_p> > &c, vector<ZZ_p> &res)
{
	//step 1: compute C and B (by dividing half)
	int tot = b.size();
	if(b.size()!=c.size())
	{
		cout<<"Invalid Param in compAnswer"<<endl;
		abort();
	}
	if(res.size()<b.size())
		res.resize(b.size());
	vector<FieldType> veca,vecb,vecres;
	int _cnt=0;
	//all length should be even, since we have length 32 here
	for(int l=0; l<tot; l++)
	{
		// if(flag_print)
		// 	cout<<"No "<<l<<", length:"<<b[l].size()<<endl;
		if(b[l].size()==1) //the answer is c_1
			res[l] = c[l][0];
		else if(b[l].size()>1)
		{
			for(int i=0; i<b[l].size(); i+=2)
			{
				veca.push_back(b[l][i]);
				vecb.push_back(c[l][i+1]);
				veca.push_back(b[l][i]);
				vecb.push_back(b[l][i+1]);
			}
		}
		else //size = 0 ?
		{
			cout<<"0-sized array detected!"<<endl;
			abort();
		}
	}
	helper->DNMultVec(veca,vecb,vecres,1);
        int nMults = vecres.size();
        for (int i=0; i<nMults; i++) {
          chkAVec.push_back(vector<FieldType>(1, veca[i]));
          chkBVec.push_back(vector<FieldType>(1, vecb[i]));
        }
	_append(chkCVec,vecres);
	//restore back to get B and C
	vector<vector<FieldType> > B,C;
	//vector<vector<ZZ_p> > B,C;
	for(int l=0; l<tot; l++)
		if(b[l].size()>1)
		{
			B.resize(B.size()+1);
			C.resize(C.size()+1);
			for(int i=0; i<b[l].size(); i+=2)
			{
				B[B.size()-1].push_back(vecres[_cnt+1]);
				C[C.size()-1].push_back(vecres[_cnt] + c[l][i]);
				_cnt+=2;
			}
		}
	//call sub to get answer
	if(B.size()>0)
		computeAnswer(B,C,veca);
	//restore answer
	for(int l=0,_t=0; l<tot; l++)
		if(b[l].size()>1)
		{
			res[l] = veca[_t];
			_t++;
		}
}
/*
	template<class FieldType>
void CompareGate<FieldType>::compP(vector<vector<FieldType> > &a, vector<FieldType> &dest)
	{
	  }
	  template<>
	  inline void CompareGate<ZZ_p>::compP(vector<vector<ZZ_p> > &a, vector<ZZ_p> &dest)
	 
{
	if(flag_print)
		cout<<"Entering compP"<<endl;
	vector<FieldType> cBits;
	//vector<ZZ_p> cBits;
	cBits.resize(eleSize);
	ZZ tmp((1ll<<61)-1);
	for(int i=eleSize-1; i>=0; i--)
	{
		cBits[i] = field->GetElement(tmp % 2);
		tmp /=2;
	}
	//step 1: compute bitwise xor
	int _cnt = a.size();
	if(dest.size() < _cnt)
		dest.resize(_cnt);
	vector<vector<FieldType> > xorSum(_cnt);
	//vector<vector<ZZ_p> > xorSum(_cnt);
	for(int l=0; l<_cnt; l++)
	{
		for(int i=0; i<eleSize; i++)
		{
			xorSum[l].resize(eleSize);
			xorSum[l][i] = a[l][i] + cBits[i];
			if(cBits[i]==1)
				xorSum[l][i] -= field->GetElement(2) * a[l][i];
		}
	}
	//step 2: compute prefix or
	vector<vector<FieldType> > prefixOr;
	vector<FieldType> E;
	//vector<vector<ZZ_p> > prefixOr;
	//vector<ZZ_p> E;
	if(flag_print)
		cout<<"Entering PrefixOr"<<endl;
	compPrefixOr(xorSum,prefixOr);
	if(flag_print)
		cout<<"Exiting compP::prefixor"<<endl;
	if(flag_print)
	  {
	  cout<<"Checking validation of prefixor"<<endl;
	  vector<ZZ_p> tri;
	  helper->openShare(1,prefixOr[0],tri);
	  cout<<"Passed"<<endl;
	  }
	E.resize(eleSize);
	for(int l=0; l<_cnt; l++)
	{
		E[0] = prefixOr[l][0];
		for(int i=1; i<eleSize; i++)
			E[i] = prefixOr[l][i] - prefixOr[l][i-1];
		dest[l] = 0;
		for(int i=0; i<eleSize; i++)
			if(cBits[i]==0)
				dest[l] = dest[l] + E[i];
	}
	if(flag_print)
	  {
	  cout<<"Checking dest"<<endl;
	  vector<ZZ_p> t1(1),t2(1);
	  t1[0] = dest;
	  helper->openShare(1,t1,t2);
	  if(t2[0]!=0 && t2[0]!=1)
	  {
	  cout<<"Invalid Comp:"<<t2[0]<<endl;
	  abort();
	  }
	  cout<<"Exiting compP"<<endl;
	  }
	//step 1: compute bitwise xor
}*/
//compare bit sharing (a<b)
// bits are given from MSB to LSB

// compare the function (x < p/2)
	template<class FieldType>
void CompareGate<FieldType>::compHalfP(vector<FieldType> &X,vector<FieldType> &dest)
	/*{
	  }

	//there are several pairs to be compared [ ? < p/2]
	template<>
	inline void CompareGate<ZZ_p>::compHalfP(vector<ZZ_p> &X, vector<ZZ_p> &dest)
	 */
{
	if(flag_print)
		cout<<"Now at compHalfP"<<endl;
	//we want to compute ([2x]_0)
	int tot = X.size();
	//step 1: generate [r]_B and [r]_p
	vector<FieldType> r;//d0,d1;
	vector<vector<FieldType> > rBits;
	//vector<ZZ_p> r;
	//vector<vector<ZZ_p> > rBits;
	//vector<ZZ_p> d0,d1;
	if(flag_print)
		cout<<"Getting random share r"<<endl;
	getRandomBitShare(tot,r,rBits);
	if(flag_print)
		cout<<"compHalfP::Successfully get shares"<<endl;
	/*if(flag_print)
	  {
	  cout<<"checking valid of bit share"<<endl;
	  vector<ZZ_p> tmp;
	  vector<ZZ> tmpres;
	  for(int i=0; i<tot; i++)
	  {
	  helper->openShare(eleSize,rBits[i],tmp);
	  ZZ _t(0);
	  for(int j=0; j<eleSize; j++)
	  {
	  cout<<"bit "<<j<<":"<<tmp[j]<<endl;
	  _t = _t * 2 + conv<ZZ>(tmp[j]);
	  }
	  if(_t > (1ll<<61))
	  {
	  cout<<"too large generated!"<<endl;
	  abort();
	  }
	  tmpres.push_back(_t);
	  }
	  helper->openShare(tot,r,tmp);
	  for(int i=0; i<tot; i++)
	  if((tmp[i]-conv<ZZ_p>(tmpres[i])) != ZZ_p(0))
	  {
	  cout<<"find difference:"<<tmp[i]<<","<<tmpres[i]<<endl;
	  }
	  }*/
	//step 2: compute [c] = [x] + [r], and reveal c
	vector<FieldType> resvec,tmpvec;
	//vector<ZZ_p> resvec,tmpvec;
	for(int i=0; i<tot; i++)
	{
		FieldType res = field->GetElement(2) * X[i] + r[i];
		//ZZ_p res = 2 * X[i] + r[i];
		resvec.push_back(res);
		/*ZZ_p res = x + r[0], c;
		  vector<ZZ_p> resvec(1),tmpvec(1);
		//resvec.resize(1); 
		resvec[0]=res;*/
	}
	helper->openTShares(tot, true, resvec,tmpvec);
	if(flag_print)
	{
		cout<<"compHalfP::Successfully opened shares"<<endl;
		cout<<"value of c:"<<tmpvec[0]<<endl;
	}
	//c = tmpvec[0];
	//step 3: compute [c <_b r]
	//vector<ZZ_p> comp_res;
	//vector<FieldType> comp_res;
	//compGiven(rBits,tmpvec,comp_res);
	compGivenNoPrefixOr(rBits,tmpvec,resvec);
	if(flag_print)
		cout<<"compGiven finished"<<endl;
	if(resvec.size()!=tot)
	{
		cout<<"Invalid compGiven length"<<endl;
		abort();
	}
	//step 4: compute the overall result
	//vector<ZZ_p> rescr0;
	vector<FieldType> rescr0;
	for(int i=0; i<tot; i++)
	{
		if(getMod(tmpvec[i],2)==0)
			rescr0.push_back(rBits[i][this->eleSize-1]);
		else rescr0.push_back(field->GetElement(1) - rBits[i][this->eleSize-1]);
	}
	//vector<ZZ_p> tmpres;
	//vector<FieldType> tmpres;
	helper->DNMultVec(rescr0,resvec,tmpvec,1);
        int nMults = tmpvec.size();
        for (int i=0; i<nMults; i++) {
          chkAVec.push_back(vector<FieldType>(1, rescr0[i]));
          chkBVec.push_back(vector<FieldType>(1, resvec[i]));
        }
	_append(chkCVec,tmpvec);
	if(flag_print)
		cout<<"compHalfP::Single Mult finished"<<endl;

	if(dest.size()<tot)
		dest.resize(tot);
	for(int i=0; i<tot; i++)
	{
		//ZZ_p tmpr = 2 * tmpres[i];
		dest[i] = field->GetElement(1) - (resvec[i] + rescr0[i] - field->GetElement(2) * tmpvec[i]);
		//dest = 1 - dest;
	}
	/*if(flag_print)
	  {
	  cout<<"compHalfP, result checking:"<<endl;
	  vector<ZZ_p> d1;
	  helper->openShare(tot,dest,d1);
	  cout<<d1[0]<<endl;
	  if(d1[0] != 0 && d1[0] != 1)
	  {
	  cout<<"Invalid CompHalfP Result"<<endl;
	  abort();
	  }
	  }*/

}
	template<class FieldType>
void CompareGate<FieldType>::getRandomBitShare(int num,vector<FieldType> &res,vector<vector<FieldType> > &bits)
{
	if(flag_print)
		cout<<"Now at GetRandomBit with "<<num<<endl;
	//for safety we generate 2 * num
	bits.resize(num);
	//vector<ZZ_p> tmp;
	//tmp.resize(num);
	res.resize(num);
	int suc;
	for(suc = 0; suc < num; suc++)
	{
		bits[suc].resize(eleSize);
		res[suc]=0;
		//generate from MSB to LSB
		//first 3 bits: 0
		if(_bitSharesValue.size() < _bitShareOffset + eleSize)
		{
			cout<<"Not enough bit share!"<<endl;
			generateBitShares(num-suc);
		}
		for(int j=0; j<blank_size; j++)
		{
			bits[suc][j] = 0;
			//bits[suc][j] = _zeroShares[_zeroShareOffset++];
			//res[suc] = res[suc] * FieldType(2) + bits[suc][j];
			/*if(_zeroShareOffset >= _zeroShares.size())
			{
				cout<<"not enough zero shares"<<endl;
				abort();
			}*/
		}
		for(int j=blank_size; j<eleSize && _bitShareOffset < _bitSharesValue.size(); j++)
		{
			bits[suc][j] = _bitSharesValue[_bitShareOffset++];
			res[suc] = res[suc] * FieldType(2) + bits[suc][j];
		}
		if(_bitSharesValue.size() < _bitShareOffset) //not enough?
		{
			cout<<"This can't happen!"<<endl;
			abort();
			if(flag_print)
				cout<<"Not enough bit share! Regenerating....."<<endl;
			generateBitShares(num-suc);
		}
	}
	/*if(flag_print)
	{
		cout<<"gen bit done, checking"<<endl;
		vector<FieldType> tmp;
		helper->openShare(num,res,tmp);
		for(int i=0; i<num; i++)
			cout<<"generated a number "<<tmp[i]<<endl;
		for(int i=0; i<num; i++)
		{
			vector<FieldType> tmp2;
			helper->openShare(eleSize,bits[i],tmp2);
			FieldType cur(0);
			for(int j=0; j<eleSize; j++)
				cur = cur * FieldType(2) + tmp2[j];
			if(cur != tmp[i])
				cout<<"difference in share:"<<cur<<","<<tmp[i]<<endl;
		}
	}*/


	/*compP(bits,chk_suc);
	  if(flag_print)
	  cout<<"genBit:compP ended"<<endl;
	  if(flag_print)
	  {
	  cout<<"Checking consistency of bit share"<<endl;
	  vector<FieldType> r1(1),r2;
	  r1[0]=res[suc];
	  helper->openShare(1,r1,r2);
	  cout<<"Check passed"<<endl;
	  }
	//TODO: should we consider multi-thread here?
	//no need to compare ??
	//IF NO NEED TO COMPARE, JUST COMMENT THE FOLLOWING LINES BY /*
	//if(flag_print)
	//	cout<<"Exiting CompP"<<endl;
	vector<FieldType> vc2;
	helper->openShare(num,chk_suc,vc2);*/
	//TODO: can we use '=' like this?
	//notice: compGiven gives [p < c], here we need to change direction
	/*
	   int tcnt=0;
	   for(int l=0; l<num && tcnt<num; l++)
	//if(vc2[l]==FieldType(0)) //valid share 
	{
	//copy
	for(int j=0; j<eleSize; j++)
	bits[tcnt][j] = bits[l][j];
	res[tcnt] = res[l];
	tcnt++;
	}
	if(tcnt<num) 
	{
	if(flag_print)
	cout<<"First try failed!"<<endl;
	vector<vector<FieldType> > _bits;
	vector<FieldType> _res;
	getRandomBitShare(num - tcnt, _res, _bits);
	for(int l=0,pos=tcnt; l<num-tcnt; l++,pos++)
	{
	res[pos] = _res[l];
	for(int j=0; j<eleSize; j++)
	bits[pos][j] = _bits[l][j];
	}
	}
	else if(flag_print)
	cout<<"First try passed!"<<endl;
	 */
	if(flag_print)
	{
		cout<<"Used #bit:"<<_bitShareOffset<<"/"<<_bitSharesValue.size()<<endl;
	}
	//res.resize(num);
	//bits.resize(num);
}


	template<class FieldType>
void CompareGate<FieldType>::getBitDecomp(FieldType &a, vector<FieldType> &dest)
{
}

	template<>
inline void CompareGate<ZZ_p>::getBitDecomp(ZZ_p &a, vector<ZZ_p> &dest)
{
	/*	//step 1: generate bit
		vector<ZZ_p> rBits,tmp,Ctmp,Cres;
		ZZ_p p = _ZZPRIME;
		ZZ_p r,c;
		getRandomBitShare(1,r,rBits);
		if(dest.size()<eleSize)
		dest.resize(eleSize);
	//step 2: reveal c
	c = a-r;
	Ctmp[0]=c;
	helper->openShare(1,Ctmp,Cres);
	c = Cres[0];
	if(c==0)
	{
	//dest = r, by coincidence
	for(int i=0; i<eleSize; i++)
	dest[i] = rBits[i];
	return;
	}
	//step 3: compute q
	ZZ_p t0 = p - c - 1,q;
	compGiven(rBits,t0,q);
	//now q is the share of compGiven result
	//step 4: restore bit decomp of g = 2^l + c - qp
	vector<ZZ_p> f1Bits,f2Bits,gBits;*/
}


template<class FieldType>
void CompareGate<FieldType>::
TruncPRSecure(vector<FieldType> &aShares,vector<FieldType> &res) {
  // a perfect secure (instead of statistical secure) version:
  // - make r a random field element instead of a ranged integer
  // - add secure comparison for overflow
  auto _t0 = high_resolution_clock::now();
  int tot = aShares.size();
  FieldType _zero = field->GetElement(0);
  FieldType _one = field->GetElement(1);
  FieldType _two = field->GetElement(2);
  vector<FieldType> rShares, r1Shares(tot, _zero);
  vector<vector<FieldType> > rBitShares; // MSB first; 128th bit fixed to zero
  getRandomBitShare(tot,rShares,rBitShares);
  // compute tail 2^m bits of r
  for (int i=0; i<tot; i++) {
    for (int j=eleSize-_m; j<eleSize; j++) {
      r1Shares[i] = r1Shares[i] * _two + rBitShares[i][j];
    }
  }

  // shift all negative fixed pt number to positive and pad w/ r
  FieldType khalf(1ull<<(_k-1)); // TODO: change this if _k > 64
  khalf= khalf * khalf;          // similar for 1ull<<_m if _m > 64
  vector<FieldType> cShares(tot, _zero), cClear;
  for (int i=0; i<tot; i++) {
    cShares[i] = (khalf + aShares[i]) + rShares[i];
  }
  helper->openTShares(tot, true, cShares, cClear);

  // compute overflow bit (r < c) by secure comparison
  vector<FieldType> overflowShares(tot, _zero);
  compGivenNoPrefixOr(rBitShares, cClear, overflowShares);
  // and compute result
  FieldType inv2ToM = _one / FieldType(1ull<<_m);
  res.resize(tot);
  for (int i=0; i<tot; i++) {
    // no overflow, c' = c1 = c mod 2^m
    FieldType c1 = field->GetElement(getMod(cClear[i], 1ull<<_m));
    // overflow, c' = c2 = (c + q) mod 2^m = (c + 2^k - 1) mod 2^m
    FieldType c2 = field->GetElement(getMod(cClear[i] - _one, 1ull<<_m));
    FieldType cpShare = (overflowShares[i]*c2 + (_one-overflowShares[i])*c1);
    res[i]= (aShares[i] - cpShare + r1Shares[i]) * inv2ToM;
  }
  auto _t1 = high_resolution_clock::now();
  TruncPR_t += duration_cast<microseconds>(_t1 - _t0).count();
}


	template<class FieldType>
void CompareGate<FieldType>::TruncPR(vector<FieldType> &a,vector<FieldType> &res)
	/*{
	  }

	  template<>
	  inline void CompareGate<ZZ_p>::TruncPR(vector<ZZ_p> &a,vector<ZZ_p> &res)
	 */
{
	auto _t0 = high_resolution_clock::now();
	int tot = a.size();
	//step 1: generate r and r'
	int m = _m; // # of decimal digits
	int k = _k; // # of range
	int kappa = _kappa;
	FieldType khalf(1ull<<(k-1));
	khalf= khalf * khalf;
	// FieldType twom(1ull<<_m);
        FieldType inv2ToM = FieldType(1) / FieldType(1ull<<_m);
	vector<FieldType> r,r1,r2,r2open;
	vector<vector<FieldType> > rBits;
	//vector<ZZ_p> r,r1,r2,r2open;
	//vector<vector<ZZ_p> > rBits;
	/*if(flag_print)
	{
		cout<<"TruncPR: checking valid"<<endl;
		helper->openShare(tot,a,r);
		cout<<"passed"<<endl;
	}*/
	getRandomBitShare(tot,r,rBits);
	if(flag_print)
		cout<<"TruncPR: get bit share done"<<endl;
	//step 2: add everything up and open them
	//notice the overall process will [not] exceed 2^120
	for(int i=0; i<tot; i++)
	{
		//first use the lower m digits to restore r1
		r1.push_back(field->GetElement(0));
		for(int j=eleSize - m; j < eleSize; j++)
			r1[i] = r1[i] * field->GetElement(2) + rBits[i][j];
		//now restore r2
		r2.push_back(field->GetElement(0));
		for(int j = eleSize - 2 * (k-1) - kappa; j < eleSize; j++)
			r2[i] = r2[i] * field->GetElement(2) + rBits[i][j];
		r2[i] = r2[i] + a[i] + khalf;
		//notice: r2 should be smaller than about 2^120, overwhelmed by a[i]
	}
	helper->openTShares(tot, true, r2,r2open);
	if(flag_print)
		cout<<"TruncPR: r2 passed"<<endl;
	if(res.size()<tot)
		res.resize(tot);
	for(int i=0; i<tot; i++)
	{
		r2open[i] = field->GetElement(getMod(r2open[i],1ull<<m));
		// res[i] = (a[i] - r2open[i] + r1[i]) * twom;
                res[i] = (a[i] - r2open[i] + r1[i]) * inv2ToM;
	}
	auto _t1 = high_resolution_clock::now();
	TruncPR_t += duration_cast<microseconds>(_t1 - _t0).count();
	/*
	   int tot = a.size();
	   if(res.size()<tot)
	   res.resize(tot);
	   for(int i=0; i<tot; i++)
	   res[i]=a[i];
	 */
}
	template<class FieldType> 
void CompareGate<FieldType>::doubleVecMult(vector<FieldType> &a,vector<FieldType> &b,vector<FieldType> &res, int groupSize)
{
	vector<FieldType> tmp;
	//step 1: normal product
	helper->DNMultVec(a,b,tmp,groupSize);
        if (groupSize == 1) {
          int nMults = tmp.size();
          for (int i=0; i<nMults; i++) {
            chkAVec.push_back(vector<FieldType>(1, a[i]));
            chkBVec.push_back(vector<FieldType>(1, b[i]));
          }
          _append(chkCVec, tmp);
        } else {
          auto aStart = a.begin();
          auto bStart = b.begin();
          int nGroups = (a.size() + groupSize -1) / groupSize;
          for (int i=0; i<nGroups-1; i++) {
            chkAVec.push_back(vector<FieldType>(aStart, aStart+groupSize));
            chkBVec.push_back(vector<FieldType>(bStart, bStart+groupSize));
            aStart += groupSize;
            bStart += groupSize;
          }
          chkAVec.push_back(vector<FieldType>(aStart, a.end()));
          chkBVec.push_back(vector<FieldType>(bStart, b.end()));
          _append(chkCVec, tmp);
        }
        //step 2: do trunc
	// TruncPR(tmp,res);
        TruncPRSecure(tmp,res);
}

	template<class FieldType>
void CompareGate<FieldType>::SoftThres(vector<FieldType> &thres, vector<FieldType> &a, vector<FieldType> &res)
	/*{
	  }

	  template<>
	  inline void CompareGate<ZZ_p>::SoftThres(vector<ZZ_p> &thres,vector<ZZ_p> &a,vector<ZZ_p> &res)
	 */
{
	//thres is [open value]
	int tot = a.size();
	vector<FieldType> b,thresPos,thresNeg;	
	FieldType pad(1ull<<(_k-1));
	pad = pad * pad; //field->GetElement(1ll<<20);
	//thres = thres + (1ll << (_k-1)); //automaically moded by q
	for(int i=0; i<tot; i++)
	{
		b.push_back(a[i]+pad);
		thresPos.push_back(thres[i] + pad);
		thresNeg.push_back(pad - thres[i]);
	}
	vector<FieldType> res1,res2,tmp3,add3;
	//test 1: a[i] ?> thres
	//compGiven(a,thresPos,res1);
	//compGiven(a,thresNeg,res2);
	//for(int i=0; i<dim; i++)
	//	res2[i] = 1-res2[i];
	/*vector<vector<FieldType> > bits1,bits2;
	bits1.resize(tot);
	bits2.resize(tot);
	for(int i=0; i<tot; i++)
	{
		getBits(thresPos[i], bits1[i]);
		getBits(thresNeg[i], bits2[i]);
	}
	compGivenNoPrefixOr(bits1, b, res1);
	compGivenNoPrefixOr(bits2, b, res2);
	{
		vector<FieldType> ttt;
		helper->openShare(tot,res2,ttt);
		for(int i=0; i<tot; i++)
		{
			cout<<"compare result:"<<i<<","<<ttt[i]<<endl;
		}
	}
	for(int i=0; i<tot; i++)
		res2[i] = FieldType(1) - res2[i];*/
	compRandom(b, thresPos,res1);
	compRandom(thresNeg, b,res2);

	/*if(flag_print)
	{
		cout<<"SoftThres: checking result of comprandom"<<endl;
		vector<FieldType> _t0,_t1;
		helper->openShare(res1.size(),res1,_t0);
		helper->openShare(res2.size(),res2,_t1);
		cout<<"comprandom passed"<<endl;
		for(int i=0; i<tot; i++)
			cout<<"result of "<<i<<":"<<_t0[i]<<","<<_t1[i]<<endl;
	}*/
	//notice res1 and res2 are 0 ot 1 (without 2^-k)
	helper->DNMultVec(res1,res2,tmp3,1);
        int nMults = tmp3.size();
        for (int i=0; i<nMults; i++) {
          chkAVec.push_back(vector<FieldType>(1, res1[i]));
          chkBVec.push_back(vector<FieldType>(1, res2[i]));
        }
	_append(chkCVec,tmp3);
	if(flag_print)
		cout<<"end of mult1"<<endl;
	helper->DNMultVec(tmp3, a, add3, 1);
        nMults = add3.size();
        for (int i=0; i<nMults; i++) {
          chkAVec.push_back(vector<FieldType>(1, tmp3[i]));
          chkBVec.push_back(vector<FieldType>(1, a[i]));
        }
	_append(chkCVec,add3);
	if(flag_print)
		cout<<"end of mult2"<<endl;
	/*for(int i=0; i<tot; i++)
	{
		tmp1.push_back(thres[i] * (field->GetElement(1)-res1[i]));
		tmp2.push_back(thres[i] * (field->GetElement(1)-res2[i]));
	}*/
	//TruncPR(tmp1,add1);
	//TruncPR(tmp2,add2);
	if(res.size()<tot)
		res.resize(tot);
	for(int i=0; i<tot; i++)
	{
		FieldType t1 = thres[i] - thres[i] * res1[i];
		FieldType t2 = thres[i] - thres[i] * res2[i];
		res[i]= a[i] + t2 - add3[i] - t1;
	}
}
	template<class FieldType>
void CompareGate<FieldType>::doubleInverse(FieldType a,FieldType &res)
	/*{
	  }

	  template<>
	  inline void CompareGate<ZZ_p>::doubleInverse(ZZ_p a,ZZ_p &res)
	 */
{
	long bas(1ll<<(2*_m));
	long r = bas / transElement(a);
	res = r;
	/*if(flag_print)
	{
		cout<<bas<<","<<a<<","<<res<<endl;
	}*/
}

template<>
inline void CompareGate<ZpMersenne127Element>::doubleInverse(ZpMersenne127Element a, ZpMersenne127Element &res)
{
	__uint128_t bas(((__uint128_t)1)<<(2*_m));
	//if(flag_print) cout<<"bas at doubleInverse:"<<ZpMersenne127Element(bas)<<endl;
	res = ZpMersenne127Element(bas / a.elem);
}

template<class FieldType>
void CompareGate<FieldType>::
runLasso(int iter,FieldType lambda, FieldType rho,
         vector<vector<FieldType> >& Ai, vector<FieldType> &bi, vector<FieldType> &res)
{
	// TODO: implement actual inputphase
	helper->preparationPhase(0, 0, 0, 0);
	const int dim = Ai.size(); //Ai: dim * dim matrix, bi: dim * 1 vector
	int N = helper->getN();
	int T = helper->getT();
	int fieldByteSize = eleSize / 8; //field->getElementSizeInBytes();	
	FieldType invN,invrho;
	FieldType ZN(N);
	doubleInverse(ZN * field->GetElement(1ull<<_m),invN);
	doubleInverse(rho * field->GetElement(1ull<<_m),invrho);
	if(flag_print)
		cout<<"InvN:"<<N<<","<<invN<<endl;

	//each party are holding #dim shares, shareOfW[i][j] means the share of the w_i[j]
	//shareofAi[i][j] means the share of Ai[j] 
	vector<vector<FieldType> > shareOfW(N),shareOfU(N),shareOfB(N);
	vector<vector<vector<FieldType> > > shareOfA(N);
	vector<FieldType> shareOfZ;
	// -- prepare the shares for the input
	auto _t01 = high_resolution_clock::now();
        // send input: flat(Ai) || bi || u (= 0) || w (= 0)
        int sendSize = dim*dim + dim*3;
        vector<vector<FieldType>> sendBufsElements;
        vector<FieldType> sendElems(sendSize, field->GetElement(0));
        for (int l1 = 0; l1<dim; l1++) {
          for (int l2 = 0; l2<dim; l2++) {
            sendElems[ l1*dim + l2 ] = Ai[l1][l2]; //
          }
          sendElems[dim * dim + l1] = bi[l1];
        }
        helper->makeTShares(sendElems, sendBufsElements);
        vector<vector<byte>> sendBufsBytes(N, vector<byte>(sendSize*fieldByteSize));
        for (int i=0; i<N; i++) {
          field->elementVectorToByteVector(sendBufsElements[i], sendBufsBytes[i]);
        }
        vector<vector<byte>> recBufBytes(N, vector<byte>(sendSize*fieldByteSize));
        helper->_comm.allToAll(sendBufsBytes, recBufBytes);
	// -- convert received bytes to shares
        vector<vector<FieldType>> recBufElements(N, vector<FieldType>(sendSize));
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < sendSize; j++) {
			recBufElements[i][j] =
				field->bytesToElement(recBufBytes[i].data() + (j * fieldByteSize));
		}
	}	
	if(flag_print)
		cout<<"received"<<endl;


	//receive shares of other parties

	for(int i=0; i<N; i++) //each
	{
		int _tot = 0;
		shareOfA[i].resize(dim);
		//Ai
		for(int l1=0; l1<dim; l1++)
		{
			shareOfA[i][l1].resize(dim);
			for(int l2=0; l2<dim; l2++)
				shareOfA[i][l1][l2] = recBufElements[i][_tot++];
		}
		//Bi
		shareOfB[i].resize(dim);
		for(int l=0; l<dim; l++)
			shareOfB[i][l] = recBufElements[i][_tot++];
		//Wi
		shareOfW[i].resize(dim);
		for(int l=0; l<dim; l++)
			shareOfW[i][l] = recBufElements[i][_tot++];
		//Ui
		shareOfU[i].resize(dim);
		for(int l=0; l<dim; l++)
			shareOfU[i][l] = recBufElements[i][_tot++];
	}
        if(flag_print) {
          cout << "sending zero shares" << endl;
        }
        //Z shares: just set to 0
        shareOfZ.resize(dim, *(field->GetZero()));

        // inputPhase

	//start iteration.
        cout << "starting iteration" << endl;
        for (int curIter = 0; curIter < iter; curIter++) {
          cout<<"Lasso: iteration " << curIter <<endl;

          runOffline();

          // iteration == segment
          vector<FieldType> tmp,tmp1,tmp2;
          for(int i=0; i<N; i++) { //repeat with each w[i]
            for(int j=0; j<dim; j++) {
              tmp.push_back(rho * (shareOfZ[j] - shareOfU[i][j]));
            }
          }
          for(int i=0,_c=0; i<N; i++)
            for(int j=0; j<dim; j++,_c++)
              tmp[_c] = tmp[_c] + shareOfB[i][j];

          tmp1.resize(N*dim*dim);
          tmp2.resize(N*dim*dim);
          for(int i=0,_c=0; i<N; i++) {
            for(int j=0; j<dim; j++)
              for(int l=0; l<dim; l++,_c++) {
                tmp1[_c] = tmp[l + i * dim];
                tmp2[_c] = shareOfA[i][j][l];
              }
          }
          vector<FieldType> tmp3;
          doubleVecMult(tmp1,tmp2,tmp3, dim);
          for(int i=0,_c=0; i<N; i++) {
            for(int j=0; j<dim; j++) { //fill W[i][j]
              shareOfW[i][j] = tmp3[_c++];
            }
          }
          tmp1.resize(dim);
          for(int j=0; j<dim; j++) {
            tmp1[j] = field->GetElement(0);
            //_t0.push_back(field->GetElement(0));
            for(int i=0; i<N; i++)
              tmp1[j] = tmp1[j] + shareOfW[i][j] + shareOfU[i][j];
            tmp1[j] = tmp1[j] * invN;
          }
          tmp2.resize(dim);
          TruncPRSecure(tmp1,tmp2);
          FieldType _tmp = lambda * invN;
          _tmp = getDiv(_tmp,1ll<<_m);
          _tmp = _tmp * invrho;
          _tmp = getDiv(_tmp,1ll<<_m);
          tmp.resize(dim);
          for(int j=0; j<dim; j++)
            tmp[j] = _tmp;
          //_t12.push_back(_tmp);
          SoftThres(tmp, tmp2, shareOfZ);

          for(int i=0; i<N; i++) {
            for(int j=0; j<dim; j++)
              shareOfU[i][j] = shareOfU[i][j] + shareOfW[i][j] - shareOfZ[j];
          }

          verificationPhase();
          cout << "verification passed for iteration " << curIter << endl;
          if ((!helper->verifySharing()) ||
              (!helper->tag()) ) {
            abort();
          }
        }
	helper->openTShares(dim, true, shareOfZ,res);
}

template <class FieldType> void CompareGate<FieldType>::readLassoInputs()
{
	if(flag_print)
		cout<<"Reading lasso"<<endl;
	ifstream myfile;
	myfile.open(inputsFile);
	int n; myfile>>n;
	if(flag_print)
		cout<<"Read "<<n<<endl;
	//read Ai
	_Ai.resize(n);
	for(int i=0; i<n; i++)
		for(int j=0; j<n; j++)
		{
			string tmp;
			myfile>>tmp;
			FieldType cur(0);
			int l=0;
			if(tmp[0]=='-') l++;
			for( ; l<tmp.length(); l++)
				cur = cur * field->GetElement(10) + field->GetElement((long)(tmp[l] - 48));
			if(tmp[0]=='-') cur = field->GetElement(0) - cur;
			_Ai[i].push_back(cur);
		}
	//read Bi
	for(int i=0; i<n; i++)
	{
		string tmp;
		myfile>>tmp;

		FieldType cur(0);
		int l=0;
		if(tmp[0]=='-') l++;
		for( ; l<tmp.length(); l++)
			cur = cur * field->GetElement(10) + field->GetElement((long)(tmp[l] - 48));
		if(tmp[0]=='-') cur = field->GetElement(0) - cur;
		_bi.push_back(cur);
	}
	myfile.close();
	if(flag_print)
		cout<<"Reading lasso ended"<<endl;
}

template <class FieldType>
void CompareGate<FieldType>::
inputPhase() {

}

template <class FieldType>
void CompareGate<FieldType>::
run() {
  int tottme = 0;
  readLassoInputs();
  vector<FieldType> res;
  for (iteration = 0; iteration < times; iteration++) {
    res.clear();
    auto t1start = high_resolution_clock::now();
    //rho = 10, lambda = 0.1
    runLasso(n_iter, field->GetElement((1ull<<(_m))/10), FieldType(10), _Ai, _bi, res);

    auto t2end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(t2end - t1start).count();
    tottme += duration;
    
  }
  writeOutput(res, "");
  cout<<"Total time for " << times << " iterations: "
      << tottme << " ms" <<endl;
}

template <class FieldType> void CompareGate<FieldType>::runOffline() {
  //if(flag_print)
  cout<<"runOffline()"<<endl;
  int dim = _Ai.size();
  int cnt = 40 * dim * dim  * n_iter * eleSize / 10;
  cnt = cnt * 2 / 80;
  cnt = cnt / 10;         // per iteration cost should roughly be 1/10

  // TODO: change for gas dataset
  int cnt_bit = 6 * dim;
  if (helper->preparationPhase(n_iter*5,
                               cnt_bit*eleSize+n_iter*2,
                               n_iter*60,
                               cnt_bit*eleSize*3) == false) {
    cout << "preparationPhase faild" << endl;
    abort();
  }
  cout<<"generating bit"<<endl;
  generateBitShares(cnt_bit);
  cout<<"bit generation done"<<endl;
}

template <class FieldType> void CompareGate<FieldType>::runOnline() {

	vector<FieldType> res;
	runLasso(n_iter, field->GetElement((1ull<<(_m))/10), FieldType(10), _Ai, _bi, res);
        
}

	template<class FieldType>
void CompareGate<FieldType>::verificationPhase()
{
  // assuming semi-honest, hence no final check performed
  int nParties = helper->getN();
  int K = helper->getK();
  for (int king=0; king<nParties; king++) {
    vector<FieldType> aShares, bShares, transcript, ultimateTranscript(7);
    helper->delinearizationHelper(king, chkAVec, chkBVec, chkCVec,
                                  aShares, bShares, transcript);
    while(helper->dimensionReduction(king, aShares, bShares, transcript) >= K);
    helper->randomization(king, aShares, bShares, transcript, ultimateTranscript);
    // TODO: do actual check and debug
    // if (!helper->checkSingleMult(king, ultimateTranscript)) {
    //   abort();
    // }
  }
  chkAVec.clear();
  chkBVec.clear();
  chkCVec.clear();
}

//only party 0 outputs
	template<class FieldType>
void CompareGate<FieldType>::writeOutput(vector<FieldType> &res,string f)
{
	if(m_partyID!=0) return;
	ofstream ouf;
	ouf.open("output" + f + ".txt");
	for(int i=0; i<res.size(); i++)
		ouf<<res[i]; //no need for endl
	ouf.close();
}
#endif /* COMPARISON_H_ */
