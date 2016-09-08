/*****************************************************************************************************************************************************/
#include <windows.h>
#include "BigNumOp.h"

/************************************************��̬����ں���*****************************************************************/
BOOL APIENTRY DllMain(HINSTANCE hInstDLL, DWORD fdwReason, LPVOID lpvReserved)
{
    switch (fdwReason)
    {
        case DLL_PROCESS_ATTACH:
            break;
        case DLL_THREAD_ATTACH:
            break;
        case DLL_THREAD_DETACH:
            break;
        case DLL_PROCESS_DETACH:
            break;
    }
    return TRUE;
}

/************************************************������������������*****************************************************************/
/*
���������������RSA�ӽ����㷨�У����ڲ�����Կ�Լ��ӽ�������.
*/
/*
���������ʵ�ַ�����Ҫ�����¼��֣�
	1)���ַ�����ʾ��������������ʮ�����ַ������ʾ��Ȼ���ա���ʽ���㡱��˼����м��㡣���ַ����Ƚ�������⣬���Ǽ���Ч�ʺܵ͡�
	2)���������ɶ����������д���ʹ�ø���λ������߼�������ʵ�ִ�������㡣�÷�����Ƹ��ӣ��ɶ��Խϲ�������Ե��ԡ�
	3)��������ʾ��һ��n�������顣n��ȡֵԽ������Ĵ�СԽС�������������������ʱ�估�ռ临�Ӷȣ�����㷨��Ч�ʡ�
	��32λϵͳ�У�n����ȡ2^32����ʱÿһλ��ȡֵ��Χ��0~0xffffffff��
 
�������Ե�3���ַ���������������RSA���漰�Ĵ���ͨ��������0������Ϊ�˼����⣬����������������д�����������0�ġ�
*/
//init,��ֵ����,��ֵΪ���� initByInt
void initByInt(BigNumber* BigNum,unsigned long long A)
{
	// ���� 2^32λ������(������ÿһ��"λ"ֻ�ܴ��2^32��С������).
	if(A>0xffffffff)
    {
		// ����A��������long long, m_nLengthΪ2.
        BigNum->m_nLength = 2;
        BigNum->m_ulValue[1] = (unsigned long)(A>>32);
        BigNum->m_ulValue[0] = (unsigned long)A;
    }
    else
    {
        BigNum->m_nLength = 1;
        BigNum->m_ulValue[0] = (unsigned long)A;
    }

	//memset(BigNum->m_ulValue + BigNum->m_nLength,0,(BINUM_MAXLEN - BigNum->m_nLength)*sizeof(unsigned long));
}

//init,��ֵ����,��ֵΪ��ͨ���� initByBigNum
void initByBigNum(BigNumber* BigNumDes,const BigNumber* BigNumSrc)
{
	unsigned int i;

	BigNumDes->m_nLength = BigNumSrc->m_nLength;
	
	for (i = 0;i < BigNumSrc->m_nLength;i++)
		BigNumDes->m_ulValue[i] = BigNumSrc->m_ulValue[i];
}

//Add,������,�����������ĺ� addByBigNum
BigNumber addByBigNum(const BigNumber* BigNumSrc1,const BigNumber* BigNumSrc2)
{
    BigNumber X;
	unsigned i;
	unsigned carry;
	// long long���ͣ����Լ򻯴����λ����
	unsigned long long sum;
	unsigned long num1;

	carry = 0;
	sum = 0;

    initByBigNum(&X,BigNumSrc1);

    if(X.m_nLength < BigNumSrc2->m_nLength)
		X.m_nLength = BigNumSrc2->m_nLength;

    for(i = 0;i<X.m_nLength;i++)
    {
		if (i >= BigNumSrc2->m_nLength)
			sum = 0;
		else
			sum = BigNumSrc2->m_ulValue[i];

		if (i >= BigNumSrc1->m_nLength)
			num1 = 0;
		else
			num1 = X.m_ulValue[i];

		sum = sum+num1+carry;
        X.m_ulValue[i]=(unsigned long)sum;
        carry=(unsigned)(sum>>32);
    }
	// �ӷ� ��λcarryֻ����0��1.
    X.m_ulValue[X.m_nLength]=carry;
    X.m_nLength+=carry;

    return X;
}

//Add,������,���������ͨ�����ĺ� addByInt
BigNumber addByInt(const BigNumber* BigNumSrc1,unsigned long long A)
{
    BigNumber X;
	unsigned long long sum;
	unsigned i;

    initByBigNum(&X,BigNumSrc1);
    
    sum = X.m_ulValue[0];
	sum += A;
    X.m_ulValue[0] = (unsigned long)sum;

    if(sum > 0xffffffff)
    {
        i = 1;

        while(X.m_ulValue[i]==0xffffffff)
		{
			X.m_ulValue[i]=0;
			i++;

			if (i >= X.m_nLength)
				X.m_ulValue[i]=0;
		}

        X.m_ulValue[i]++;

		if(BigNumSrc1->m_nLength == i)
			X.m_nLength++;
    }

    return X;
}

//Sub,������,�����������Ĳ� subByBigNum
BigNumber subByBigNum(const BigNumber* BigNumSrc1,const BigNumber* BigNumSrc2)
{
	unsigned carry;
	unsigned long numt;
    unsigned long long num;
	unsigned i;
    BigNumber X;

	carry = 0;

	//initByInt(&X,0);
    initByBigNum(&X,BigNumSrc1);

	// ������ С�� ����, ������ֱ����0.
    if(Cmp(&X,BigNumSrc2)<=0)
	{
		initByInt(&X,0);
		return X;
	}

    for(i = 0;i < BigNumSrc1->m_nLength;i++)
    {
		if (i >= BigNumSrc2->m_nLength)
			numt = 0;
		else
			numt = BigNumSrc2->m_ulValue[i];

        if((BigNumSrc1->m_ulValue[i] > numt) || 
			((BigNumSrc1->m_ulValue[i] == numt) && (carry==0)))
        {
			// ���ý�λ
            X.m_ulValue[i] = BigNumSrc1->m_ulValue[i] - carry - numt;
            carry=0;
        }
        else
        {
			// ��1λ"0x100000000LL".
            num = 0x100000000LL + BigNumSrc1->m_ulValue[i];
            X.m_ulValue[i]=(unsigned long)(num - carry - numt);
            carry=1;
        }
    }

    while(X.m_ulValue[X.m_nLength-1]==0)
		X.m_nLength--;
    
	return X;
}

//Sub,������,���������ͨ�����Ĳ� subByInt
BigNumber subByInt(const BigNumber* BigNumSrc1,unsigned long A)
{
    BigNumber X;
	unsigned long long num;
	int i;

    initByBigNum(&X,BigNumSrc1);

    if(X.m_ulValue[0] >= A)
	{
		X.m_ulValue[0] -= A;
		return X;
	}

    if(X.m_nLength == 1)
	{
		initByInt(&X,0);
		return X;
	}

	num = 0x100000000LL;
    num += X.m_ulValue[0];
    X.m_ulValue[0] = (unsigned long)(num - A);

    i = 1;
    while(X.m_ulValue[i] == 0)
	{
		X.m_ulValue[i] = 0xffffffff;
		i++;
	}

    X.m_ulValue[i]--;
    if(X.m_ulValue[i]==0)
		X.m_nLength--;
    
	return X;
}

//Mul,������,���������������Ļ� mulByBigNum
BigNumber mulByBigNum(const BigNumber* BigNumSrc1,const BigNumber* BigNumSrc2)
{
	BigNumber X;
	unsigned long long sum;
	unsigned long long mul;
	unsigned long long carry;
	unsigned i;
	unsigned j;

	mul=0;
	carry=0;
	initByInt(&X,0);

	// ��һ������Ϊ0,ֱ�ӷ���0.
	if ((BigNumSrc1->m_nLength == 1 && BigNumSrc1->m_ulValue[0] == 0) || (BigNumSrc2->m_nLength == 1 && BigNumSrc2->m_ulValue[0] == 0))
		return X;

	// ��������ͨ�����Ļ�
    if(BigNumSrc2->m_nLength==1)
		return mulByInt(BigNumSrc1,BigNumSrc2->m_ulValue[0]);

	X.m_nLength = BigNumSrc1->m_nLength + BigNumSrc2->m_nLength - 1;

	/* ����ÿһλ����ֵ�� ���㷽ʽ:
		����i=4,��Ϊ���� �ڻ��ĵ�5λ����ֵ.
	*/
    for(i = 0;i < X.m_nLength;i++)
	{
		sum = carry;
		carry = 0;

		for(j = 0;j < BigNumSrc2->m_nLength;j++)
		{
            if(((i - j) >= 0) && ((i - j) < BigNumSrc1->m_nLength))
			{
				mul = BigNumSrc1->m_ulValue[i - j];
				mul *= BigNumSrc2->m_ulValue[j];
			    carry += mul>>32;
				mul = mul&0xffffffff;
				sum += mul;
			}
        }

		carry+=sum>>32;
		X.m_ulValue[i]=(unsigned long)sum;
	}

	if(carry)
	{
		X.m_nLength++;
		X.m_ulValue[X.m_nLength-1] = (unsigned long)carry;
	}

    return X;
}

//Mul,������,���������ͨ�����Ļ� mulByInt
BigNumber mulByInt(const BigNumber* BigNumSrc1,unsigned long A)
{
    BigNumber X;
    unsigned long long mul;
    unsigned long carry;
	unsigned i;

	carry = 0;

    initByBigNum(&X,BigNumSrc1);

    for(i = 0;i < BigNumSrc1->m_nLength;i++)
    {
        mul = BigNumSrc1->m_ulValue[i];
        mul = mul*A + carry;
        X.m_ulValue[i] = (unsigned long)mul;
        carry = (unsigned long)(mul>>32);
    }

    if(carry)
	{
		X.m_nLength++;
		X.m_ulValue[X.m_nLength-1] = carry;
	}

    return X;
}

//Div,������,�������������� divByBigNum
BigNumber divByBigNum(const BigNumber* BigNumSrc1,const BigNumber* BigNumSrc2)
{
	BigNumber Temp1;
	BigNumber X;
	BigNumber Y; 
	BigNumber Z;
    unsigned i;
	unsigned len;
    unsigned long long num;
	unsigned long long div;

	initByInt(&X,0);

	if(BigNumSrc2->m_nLength == 1)
		return divByInt(BigNumSrc1,BigNumSrc2->m_ulValue[0]);
    
    initByBigNum(&Y,BigNumSrc1);

    while(Cmp(&Y,BigNumSrc2) >= 0)
    {
		div = Y.m_ulValue[Y.m_nLength-1];
		num = BigNumSrc2->m_ulValue[BigNumSrc2->m_nLength - 1];
		len = Y.m_nLength - BigNumSrc2->m_nLength;
		
		// �������ͳ����� ���λ��λ��һ��������Ϊ1��ȡ�����ʶ���ʱ���ÿ������λ֮�������λ��.
		if((div == num) && (len == 0))
		{
			X = addByInt(&X,1);
			break;
		}

		if((div <= num) && len)
		{
			// ������С�ڳ�������λ. 
			len--;
			div = (div<<32) + Y.m_ulValue[Y.m_nLength-2];
		}	

		// ��num+1������ֹnumΪ0, ��Ϊdiv>num,���ԡ�num+1�����ᵼ���̵ı仯.
		div = div/(num+1);
		// Z���� ÿһ�γ�������.
		initByInt(&Z,div);
		if(len)
		{
			// ��Z����λ���������Ա����潫ÿһλ�����������.
			Z.m_nLength += len;
			for(i = Z.m_nLength-1;i >= len;i--)
				Z.m_ulValue[i] = Z.m_ulValue[i-len];

			for(i = 0;i < len;i++)
				Z.m_ulValue[i]=0;
		}
		// ����Z��������ע�⣬Z�Ѿ�����λ��������.
		X = addByBigNum(&X,&Z);
		// ��ȡ������������
		Temp1 = mulByBigNum(BigNumSrc2,&Z);
		// YΪ������������Y������һ�ֵĳ�������.
		Y = subByBigNum(&Y,&Temp1);
    }

    return X;
}

//Div,������,���������ͨ�������� divByInt
BigNumber divByInt(const BigNumber* BigNumSrc1,unsigned long A)
{
	int i;
    BigNumber X;
	unsigned long long div;
	unsigned long long mul;
    unsigned long carry;

	carry = 0;

    initByBigNum(&X,BigNumSrc1);

    if(X.m_nLength==1)
	{
		X.m_ulValue[0]=X.m_ulValue[0]/A;
		return X;
	}
 
	// ��������ͨ�������̣������ֵ��λ��������߼�1.
    for(i = X.m_nLength - 1;i >= 0;i--)
    {
        div = carry;
		// div<<32,����������Ϊ��һ�γ����ĸ�λ.
        div = (div<<32) + X.m_ulValue[i];
        X.m_ulValue[i] = (unsigned long)(div/A);
        mul = (div/A)*A;
		// ��������
        carry = (unsigned long)(div-mul);
    }

    if(X.m_ulValue[X.m_nLength-1] == 0)
		X.m_nLength--;
    
	return X;
}

//Mod,ģ����,������������ģ modByBigNum
BigNumber modByBigNum(const BigNumber* BigNumSrc1,const BigNumber* BigNumSrc2)
{
    BigNumber X;
	BigNumber Y;
	unsigned long long div;
	unsigned long long num;
	unsigned i;
	unsigned len;

    initByBigNum(&X,BigNumSrc1);
    
	while(Cmp(&X,BigNumSrc2) >= 0)
    {
		div = X.m_ulValue[X.m_nLength - 1];
		num = BigNumSrc2->m_ulValue[BigNumSrc2->m_nLength - 1];
		len = X.m_nLength -  BigNumSrc2->m_nLength;
		
		// ���λ��λ��һ�������ÿ�������λ��ֱ��ȡ��.
		if((div == num) && (len == 0))
		{
			X = subByBigNum(&X,BigNumSrc2);
			break;
		}

		if((div <= num) && len)
		{
			len--;
			div = (div<<32) + X.m_ulValue[X.m_nLength-2];
		}

		div = div/(num+1);
		// ���� BigNumSrc2��������.
		initByInt(&Y,div);
		Y = mulByBigNum(BigNumSrc2,&Y);

		if(len)
		{
			// λ������
			Y.m_nLength += len;
			for(i = Y.m_nLength - 1;i >= len;i--)
				Y.m_ulValue[i] = Y.m_ulValue[i-len];

			for(i = 0;i < len;i++)
				Y.m_ulValue[i]=0;
		}
		// ����������һ��ȡ��
		X = subByBigNum(&X,&Y);
    }

    return X;
}
 
//Mod,ģ����,���������ͨ������ģ modByInt
unsigned long modByInt(const BigNumber* BigNumSrc1,unsigned long A)
{
	int i;
	unsigned long long div;
    unsigned long carry;

	carry = 0;

    if(BigNumSrc1->m_nLength == 1)
		return(BigNumSrc1->m_ulValue[0]%A);
    
    for(i = BigNumSrc1->m_nLength - 1;i >= 0;i--)
    {
        div = BigNumSrc1->m_ulValue[i];
		div += carry*0x100000000LL;
        carry = (unsigned long)(div%A);
    }

    return carry;
}

//Cmp,�Ƚ�����
int Cmp(const BigNumber* BigNumSrc1,const BigNumber* BigNumSrc2)
{
	int i;
	
    if(BigNumSrc1->m_nLength > BigNumSrc2->m_nLength)
		return 1;
	
    if(BigNumSrc1->m_nLength < BigNumSrc2->m_nLength)
		return -1;
	
    for(i = BigNumSrc1->m_nLength-1;i >= 0;i--)
    {
        if(BigNumSrc1->m_ulValue[i] > BigNumSrc2->m_ulValue[i])
			return 1;
		
        if(BigNumSrc1->m_ulValue[i] < BigNumSrc2->m_ulValue[i])
			return -1;
    }
	
    return 0;
}

//���BigNum==0����1�����򷵻�0. CmpZero
int CmpZero(BigNumber* BigNum)
{
	if (BigNum->m_nLength == 1 && BigNum->m_ulValue[0] == 0)
	{
		return 1;
	}
	else
		return 0;
}

//��������ת��Ϊ��������
void BigNumToIndex(BigNumber* BN, int b[], int* count)
{
	BigNumber BigNum;
	
	*count = 0;
	initByBigNum(&BigNum, BN);
	while (!(BigNum.m_nLength == 1 && BigNum.m_ulValue[0] == 0))
	{
		b[*count]=modByInt(&BigNum,2);
		*count += 1;
		BigNum = divByInt(&BigNum,2);
	}
}

//(BN1 - BN2) % N �Ĵ���
BigNumber SubModProc(BigNumber* BN1,BigNumber* BN2,BigNumber* N)
{
	BigNumber BNTempn,BNTemp,BNTemp1,BNTemp2;

	// BN1 >= BN2.
	if (Cmp(BN1,BN2) >= 0)
		BNTemp = subByBigNum(BN1,BN2);
	else
	{
		// BN1 < BN2.
		BNTemp1 = subByBigNum(BN2,BN1);

		if (Cmp(&BNTemp1,N))
		{
			BNTemp2 = divByBigNum(&BNTemp1,N);
			BNTempn = addByInt(&BNTemp2,1);
			BNTemp2 = mulByBigNum(N,&BNTempn);
			BNTemp = subByBigNum(&BNTemp2,&BNTemp1);
		}
		else
			BNTemp = subByBigNum(N,&BNTemp1);
	}
	// ȡ��
	BNTemp1 = modByBigNum(&BNTemp,N);

	return BNTemp1;
}

/****************************************************************************************
һ��RSA��Կ���ɵĲ���:
��һ�������ѡ����������ȵ�����p��q(ʵ��Ӧ���У�����������Խ�󣬾�Խ���ƽ�);
�ڶ���������p��q�ĳ˻�n(n�Ķ����Ƴ��Ⱦ�����Կ����);
������������n��ŷ��������(n);(���ݺ���������(n) = (p-1)(q-1));
���Ĳ������ѡ��һ������e��������1< e < ��(n)����e���(n) ����;
���岽������e���ڦ�(n)��ģ��Ԫ��d.��ν"ģ��Ԫ��"����ָ��һ������d������ʹ��ed����(n)��������Ϊ1����ed �� 1 (mod ��(n));
����������n��e��װ�ɹ�Կ��n��d��װ��˽Կ��

����岽�е�ģ��Ԫ�أ�����ʹ�������ŷ������㷨��շת���������
Euc,ŷ������㷨���ͬ�෽��,�󲻶�����ax-by=1����С������
����ֵ��X,���㣺BigNumSrc1*X mod A=1
****************************************************************************************/
BigNumber Euc(const BigNumber* BigNumSrc1,const BigNumber* BigNumSrc2)
{
	int i;
	BigNumber M,E,X,Y,I,J;
    int x,y;

	initByBigNum(&M,BigNumSrc2);
	initByBigNum(&E,BigNumSrc1);
	initByInt(&X,0);
	initByInt(&Y,1);
	x = y = 1;

	i = 0;
	while((E.m_nLength!=1) || (E.m_ulValue[0]!=0))
	{
		i++;
		I = divByBigNum(&M,&E);
		J = modByBigNum(&M,&E);
		initByBigNum(&M,&E);
		initByBigNum(&E,&J);
		initByBigNum(&J,&Y);
		Y = mulByBigNum(&Y,&I);
		
		if (x == y)
		{
		    if(Cmp(&X,&Y) >= 0)
			{
				Y = subByBigNum(&X,&Y);
			}
			else
			{
				Y = subByBigNum(&Y,&X);
				y = 0;
			}
		}
		else
		{
			Y = addByBigNum(&X,&Y);

			x = 1 - x;
			y = 1 - y;
		}

		initByBigNum(&X,&J);
	}

	if (x == 0)
	{
		X = subByBigNum(BigNumSrc2,&X);
	}
	
	return X;
}

/****************************************************************************************
ģ�ظ�ƽ���㷨���ǵݹ飩ʵ��˼·��b^n mod m����
(1) a��ʼ��Ϊ1����ָ��nת��Ϊ�����ƣ��Զ����Ƶ�ÿһλ�������ж�;
ע��������λ���ԴӸ�λ����λ�����λ����λ��˳�� ���жϣ������ǴӸ�λ��ʼ�жϡ�
(2) b = b*b (mod m);
(3) ���Ϊ1,��a = a*b (mod m); ����a����;
������ϸ���Ͳμ���http://blog.csdn.net/yxtxiaotian/article/details/52464496
�Ǵ�������λ��ʼ�жϵ�C������Բμ���https://github.com/yangxt225/powerMod

****************************************************************************************
// ����output���������ڷ�װ��̬���ǵ���ȡ�ý����
RsaTransTest,����ƽ���㷨������ģ����,��˷���ģ
N = BigNumN,A = BigNumSrc1,B = BigNumSrc2
����ֵ��X=N^A MOD B
****************************************************************************************/
int RsaTransTest(const BigNumber* BigNumN,const BigNumber* BigNumSrc1,const BigNumber* BigNumSrc2, BigNumber* output)
{
    BigNumber X,Y,Temp;
	int i,j,k;
	unsigned long n;
	unsigned long num;

	// ָ���Ķ�����ʵ����Чλ
	k = BigNumSrc1->m_nLength*32 - 32; 		
	num = BigNumSrc1->m_ulValue[BigNumSrc1->m_nLength - 1];
	while(num)
	{
		num = num>>1;
		k++;
	}
	
	// X = N
	initByBigNum(&X,BigNumN);
	
	/* i��ʼ��Ϊk-2,��һ:�����Ķ������±���Ҫ��1, 
	���:��������߶�����λ(�����������Ƶĵ�k-1λ)һ��Ϊ1.
	iΪʲô��ʼ��Ϊk-2����
	*/
	for(i = k-2;i >= 0;i--)
	{
		/*�������λ���� b=b*b Mod m;
		�൱�ڣ�X = X*X Mod BigNumSrc2;
		*/ 
		Y = mulByInt(&X,X.m_ulValue[X.m_nLength-1]);
		Y = modByBigNum(&Y,BigNumSrc2);

		// ��ȥ���λ(n=1)���ٶԴ�����ÿһλ��ִ������ b*b Mod m
        for(n = 1;n < X.m_nLength;n++)
		{   
			/*��Y������32��������Ϊ32���ƣ�,Ҳ������1λ��
			��Ϊ����������򣬴Ӹ�λ�̳ж���ֵ(Yֵ)���ڵ�λ��ʱ����Ҫ��չΪ(���Ƶı���)32����
			*/ 
			for(j = Y.m_nLength;j > 0;j--)
				Y.m_ulValue[j] = Y.m_ulValue[j - 1];
			Y.m_ulValue[0] = 0;
			Y.m_nLength++;
			// 
			Temp = mulByInt(&X,X.m_ulValue[X.m_nLength-n-1]);
			Y = addByBigNum(&Y,&Temp);
			Y = modByBigNum(&Y,BigNumSrc2);
		}
		// ִ�и�ֵ b=b*b Mod m
		initByBigNum(&X,&Y);
		
		/* if�����жϣ�
			"i>>5" i��ʾ������λ����ţ�����5λ(����32)�õ�N,��ʾ����ֵ�ĵ�Nλ(Dec).
			"i&31" ��ȡ���Ĵ���ֵ��Nλ����ֵ����ȡ���������ʽ�ĵ�"i&31"λ.
			�ж��Ƿ�Ϊ1������ǣ�ִ������.
		*/
		if((BigNumSrc1->m_ulValue[i>>5]>>(i&31)) & 1)
		{
			/* ִ������ a=a*b Mod m;
			�൱�ڣ�X = BigNumN*X Mod BigNumSrc2;
			*/
			// ��BigNumN��ˣ���
			Y = mulByInt(BigNumN,X.m_ulValue[X.m_nLength-1]);
			Y = modByBigNum(&Y,BigNumSrc2);

            for(n = 1;n < X.m_nLength;n++)
			{          
			    for(j = Y.m_nLength;j > 0;j--)
					Y.m_ulValue[j] = Y.m_ulValue[j-1];
			    Y.m_ulValue[0] = 0;
			    Y.m_nLength++;
				// 
				Temp = mulByInt(BigNumN,X.m_ulValue[X.m_nLength-n-1]);
				Y = addByBigNum(&Y,&Temp);
				Y = modByBigNum(&Y,BigNumSrc2);
			}
			// ���Ǵ����X ����
		    initByBigNum(&X,&Y);
		}
	}
	
	initByBigNum(output, &X);
    return 0;
}





