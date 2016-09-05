/*****************************************************************************************************************************************************/
#ifndef _BIG_NUM_H
#define _BIG_NUM_H

#ifdef __cplusplus
extern "C" {
#endif

//�����������,�ɸ��ݴ���λ������,����Ӱ���ٶ�
#define BINUM_MAXLEN 128

//�����洢�ṹ
typedef struct 
{
	//������0x10000000�����µĳ���
	unsigned m_nLength;

	//�������¼������0x10000000������ÿһλ��ֵ
    unsigned long m_ulValue[BINUM_MAXLEN];
}BigNumber;

/************************************************������������������*****************************************************************/

//init,��ֵ����,��ֵΪ����
void initByInt(BigNumber* BigNum,unsigned long long A);

//init,��ֵ����,��ֵΪ��ͨ����
void initByBigNum(BigNumber* BigNumDes,const BigNumber* BigNumSrc);

//Add,������,�����������ĺ�
BigNumber addByBigNum(const BigNumber* BigNumSrc1,const BigNumber* BigNumSrc2);

//Add,������,���������ͨ�����ĺ�
BigNumber addByInt(const BigNumber* BigNumSrc1,unsigned long long A);

//Sub,������,�����������Ĳ�
BigNumber subByBigNum(const BigNumber* BigNumSrc1,const BigNumber* BigNumSrc2);

//Sub,������,���������ͨ�����Ĳ�
BigNumber subByInt(const BigNumber* BigNumSrc1,unsigned long A);

//Mul,������,���������������Ļ�
BigNumber mulByBigNum(const BigNumber* BigNumSrc1,const BigNumber* BigNumSrc2);

//Mul,������,���������ͨ�����Ļ�
BigNumber mulByInt(const BigNumber* BigNumSrc1,unsigned long A);

//Div,������,��������������
BigNumber divByBigNum(const BigNumber* BigNumSrc1,const BigNumber* BigNumSrc2);

//Div,������,���������ͨ��������
BigNumber divByInt(const BigNumber* BigNumSrc1,unsigned long A);

//Mod,ģ����,������������ģ
BigNumber modByBigNum(const BigNumber* BigNumSrc1,const BigNumber* BigNumSrc2);

//Mod,ģ����,���������ͨ������ģ
unsigned long modByInt(const BigNumber* BigNumSrc1,unsigned long A);

//Cmp,�Ƚ����� 1:> 0:= -1:<
int Cmp(const BigNumber* BigNumSrc1,const BigNumber* BigNumSrc2);

//���BigNum==0����1�����򷵻�0
int CmpZero(BigNumber* BigNum);

//��������ת��Ϊ��������
void BigNumToIndex(BigNumber* BN,int b[],int* count);

//(BN1 - BN2) % N �Ĵ���
BigNumber SubModProc(BigNumber* BN1,BigNumber* BN2,BigNumber* N);

//Euc,ŷ������㷨���ͬ�෽��,�󲻶�����ax-by=1����С������,����ֵ��X,���㣺NX mod A=1
BigNumber Euc(const BigNumber* BigNumSrc1,const BigNumber* BigNumSrc2);
