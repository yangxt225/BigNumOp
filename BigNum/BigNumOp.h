/*****************************************************************************************************************************************************/
#ifndef _BIG_NUM_H
#define _BIG_NUM_H

#ifdef __cplusplus
extern "C" {
#endif

//存大数的数组,可根据大数位数定义,以免影响速度
#define BINUM_MAXLEN 128

//大数存储结构
typedef struct 
{
	//大数在0x10000000进制下的长度
	unsigned m_nLength;

	//用数组记录大数在0x10000000进制下每一位的值
    unsigned long m_ulValue[BINUM_MAXLEN];
}BigNumber;

/************************************************大数基本操作与运算*****************************************************************/

//init,赋值运算,赋值为大数
void initByInt(BigNumber* BigNum,unsigned long long A);

//init,赋值运算,赋值为普通整数
void initByBigNum(BigNumber* BigNumDes,const BigNumber* BigNumSrc);

//Add,加运算,求大数与大数的和
BigNumber addByBigNum(const BigNumber* BigNumSrc1,const BigNumber* BigNumSrc2);

//Add,加运算,求大数与普通整数的和
BigNumber addByInt(const BigNumber* BigNumSrc1,unsigned long long A);

//Sub,减运算,求大数与大数的差
BigNumber subByBigNum(const BigNumber* BigNumSrc1,const BigNumber* BigNumSrc2);

//Sub,减动算,求大数与普通整数的差
BigNumber subByInt(const BigNumber* BigNumSrc1,unsigned long A);

//Mul,乘运算,求大数与大数整数的积
BigNumber mulByBigNum(const BigNumber* BigNumSrc1,const BigNumber* BigNumSrc2);

//Mul,乘运算,求大数与普通整数的积
BigNumber mulByInt(const BigNumber* BigNumSrc1,unsigned long A);

//Div,除运算,求大数与大数的商
BigNumber divByBigNum(const BigNumber* BigNumSrc1,const BigNumber* BigNumSrc2);

//Div,除运算,求大数与普通整数的商
BigNumber divByInt(const BigNumber* BigNumSrc1,unsigned long A);

//Mod,模运算,求大数与大数的模
BigNumber modByBigNum(const BigNumber* BigNumSrc1,const BigNumber* BigNumSrc2);

//Mod,模运算,求大数与普通整数的模
unsigned long modByInt(const BigNumber* BigNumSrc1,unsigned long A);

//Cmp,比较运算 1:> 0:= -1:<
int Cmp(const BigNumber* BigNumSrc1,const BigNumber* BigNumSrc2);

//如果BigNum==0返回1，否则返回0
int CmpZero(BigNumber* BigNum);

//将大整数转化为二进制数
void BigNumToIndex(BigNumber* BN,int b[],int* count);

//(BN1 - BN2) % N 的处理
BigNumber SubModProc(BigNumber* BN1,BigNumber* BN2,BigNumber* N);

//Euc,欧几里德算法求解同余方程,求不定方程ax-by=1的最小整数解,返回值：X,满足：NX mod A=1
BigNumber Euc(const BigNumber* BigNumSrc1,const BigNumber* BigNumSrc2);
