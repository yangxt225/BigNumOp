/*****************************************************************************************************************************************************/
#include <windows.h>
#include "BigNumOp.h"

/************************************************动态库入口函数*****************************************************************/
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

/************************************************大数基本操作与运算*****************************************************************/
/*
大数运算可以用在RSA加解密算法中，用于产生密钥对及加解密运算.
*/
/*
大数运算的实现方法主要有以下几种：
	1)用字符串表示大数。将大数用十进制字符数组表示，然后按照“竖式计算”的思想进行计算。这种方法比较容易理解，但是计算效率很低。
	2)将大数看成二进制流进行处理。使用各种位运算和逻辑操作来实现打算的运算。该方法设计复杂，可读性较差，而且难以调试。
	3)将大数表示成一个n进制数组。n的取值越大，数组的大小越小，这样可以缩短运算的时间及空间复杂度，提高算法的效率。
	在32位系统中，n可以取2^32，这时每一位的取值范围是0~0xffffffff。
 
下面就针对第3）种方法进行描述。在RSA中涉及的大数通常都大于0，所以为了简化问题，假设运算过程中所有大数均都大于0的。
*/
//init,赋值运算,赋值为大数 initByInt
void initByInt(BigNumber* BigNum,unsigned long long A)
{
	// 大于 2^32位的数据(大数的每一个"位"只能存放2^32大小的数据).
	if(A>0xffffffff)
    {
		// 这里A的类型是long long, m_nLength为2.
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

//init,赋值运算,赋值为普通整数 initByBigNum
void initByBigNum(BigNumber* BigNumDes,const BigNumber* BigNumSrc)
{
	unsigned int i;

	BigNumDes->m_nLength = BigNumSrc->m_nLength;
	
	for (i = 0;i < BigNumSrc->m_nLength;i++)
		BigNumDes->m_ulValue[i] = BigNumSrc->m_ulValue[i];
}

//Add,加运算,求大数与大数的和 addByBigNum
BigNumber addByBigNum(const BigNumber* BigNumSrc1,const BigNumber* BigNumSrc2)
{
    BigNumber X;
	unsigned i;
	unsigned carry;
	// long long类型，可以简化处理进位问题
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
	// 加法 进位carry只会是0或1.
    X.m_ulValue[X.m_nLength]=carry;
    X.m_nLength+=carry;

    return X;
}

//Add,加运算,求大数与普通整数的和 addByInt
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

//Sub,减运算,求大数与大数的差 subByBigNum
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

	// 被减数 小于 减数, 被减数直接置0.
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
			// 不用借位
            X.m_ulValue[i] = BigNumSrc1->m_ulValue[i] - carry - numt;
            carry=0;
        }
        else
        {
			// 借1位"0x100000000LL".
            num = 0x100000000LL + BigNumSrc1->m_ulValue[i];
            X.m_ulValue[i]=(unsigned long)(num - carry - numt);
            carry=1;
        }
    }

    while(X.m_ulValue[X.m_nLength-1]==0)
		X.m_nLength--;
    
	return X;
}

//Sub,减运算,求大数与普通整数的差 subByInt
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

//Mul,乘运算,求大数与大数整数的积 mulByBigNum
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

	// 有一个因子为0,直接返回0.
	if ((BigNumSrc1->m_nLength == 1 && BigNumSrc1->m_ulValue[0] == 0) || (BigNumSrc2->m_nLength == 1 && BigNumSrc2->m_ulValue[0] == 0))
		return X;

	// 大数与普通整数的积
    if(BigNumSrc2->m_nLength==1)
		return mulByInt(BigNumSrc1,BigNumSrc2->m_ulValue[0]);

	X.m_nLength = BigNumSrc1->m_nLength + BigNumSrc2->m_nLength - 1;

	/* 积的每一位的数值的 计算方式:
		比如i=4,则为计算 在积的第5位的数值.
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

//Mul,乘运算,求大数与普通整数的积 mulByInt
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

//Div,除运算,求大数与大数的商 divByBigNum
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
		
		// 被除数和除数的 最高位和位数一样，则商为1（取整，故而此时不用考虑最高位之外的其他位）.
		if((div == num) && (len == 0))
		{
			X = addByInt(&X,1);
			break;
		}

		if((div <= num) && len)
		{
			// 被除数小于除数，借位. 
			len--;
			div = (div<<32) + Y.m_ulValue[Y.m_nLength-2];
		}	

		// “num+1”：防止num为0, 因为div>num,所以“num+1”不会导致商的变化.
		div = div/(num+1);
		// Z保存 每一次除法的商.
		initByInt(&Z,div);
		if(len)
		{
			// 对Z进行位数调整，以便下面将每一位的商组合起来.
			Z.m_nLength += len;
			for(i = Z.m_nLength-1;i >= len;i--)
				Z.m_ulValue[i] = Z.m_ulValue[i-len];

			for(i = 0;i < len;i++)
				Z.m_ulValue[i]=0;
		}
		// 将商Z加起来（注意，Z已经作了位数调整）.
		X = addByBigNum(&X,&Z);
		// 获取除数的整数倍
		Temp1 = mulByBigNum(BigNumSrc2,&Z);
		// Y为余数，对余数Y进入下一轮的除法运算.
		Y = subByBigNum(&Y,&Temp1);
    }

    return X;
}

//Div,除运算,求大数与普通整数的商 divByInt
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
 
	// 大数与普通整数的商，结果数值的位数不变或者减1.
    for(i = X.m_nLength - 1;i >= 0;i--)
    {
        div = carry;
		// div<<32,余数部分作为下一次除法的高位.
        div = (div<<32) + X.m_ulValue[i];
        X.m_ulValue[i] = (unsigned long)(div/A);
        mul = (div/A)*A;
		// 余数部分
        carry = (unsigned long)(div-mul);
    }

    if(X.m_ulValue[X.m_nLength-1] == 0)
		X.m_nLength--;
    
	return X;
}

//Mod,模运算,求大数与大数的模 modByBigNum
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
		
		// 最高位和位数一样，不用考虑其他位，直接取余.
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
		// 计算 BigNumSrc2的整数倍.
		initByInt(&Y,div);
		Y = mulByBigNum(BigNumSrc2,&Y);

		if(len)
		{
			// 位数调整
			Y.m_nLength += len;
			for(i = Y.m_nLength - 1;i >= len;i--)
				Y.m_ulValue[i] = Y.m_ulValue[i-len];

			for(i = 0;i < len;i++)
				Y.m_ulValue[i]=0;
		}
		// 余数进入下一次取余
		X = subByBigNum(&X,&Y);
    }

    return X;
}
 
//Mod,模运算,求大数与普通整数的模 modByInt
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

//Cmp,比较运算
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

//如果BigNum==0返回1，否则返回0. CmpZero
int CmpZero(BigNumber* BigNum)
{
	if (BigNum->m_nLength == 1 && BigNum->m_ulValue[0] == 0)
	{
		return 1;
	}
	else
		return 0;
}

//将大整数转化为二进制数
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

//(BN1 - BN2) % N 的处理
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
	// 取余
	BNTemp1 = modByBigNum(&BNTemp,N);

	return BNTemp1;
}

/****************************************************************************************
一、RSA密钥生成的步骤:
第一步，随机选择两个不相等的质数p和q(实际应用中，这两个质数越大，就越难破解);
第二步，计算p和q的乘积n(n的二进制长度就是密钥长度);
第三步，计算n的欧拉函数φ(n);(根据函数：　φ(n) = (p-1)(q-1));
第四步，随机选择一个整数e，条件是1< e < φ(n)，且e与φ(n) 互质;
第五步，计算e对于φ(n)的模反元素d.所谓"模反元素"就是指有一个整数d，可以使得ed被φ(n)除的余数为1。即ed ≡ 1 (mod φ(n));
第六步，将n和e封装成公钥，n和d封装成私钥。

求第五步中的模反元素，可以使用下面的欧几里得算法（辗转相除法）：
Euc,欧几里德算法求解同余方程,求不定方程ax-by=1的最小整数解
返回值：X,满足：BigNumSrc1*X mod A=1
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
模重复平方算法（非递归）实现思路（b^n mod m）：
(1) a初始化为1，将指数n转换为二进制，对二进制的每一位，进行判断;
注：二进制位可以从高位到低位，或低位到高位的顺序 来判断，本例是从高位开始判断。
(2) b = b*b (mod m);
(3) 如果为1,则a = a*b (mod m); 否则，a不变;
更多详细解释参见：http://blog.csdn.net/yxtxiaotian/article/details/52464496
非大数，低位开始判断的C代码可以参见：https://github.com/yangxt225/powerMod

****************************************************************************************
// 传入output参数，用于封装动态库是调用取得结果。
RsaTransTest,反复平方算法进行幂模运算,求乘方的模
N = BigNumN,A = BigNumSrc1,B = BigNumSrc2
返回值：X=N^A MOD B
****************************************************************************************/
int RsaTransTest(const BigNumber* BigNumN,const BigNumber* BigNumSrc1,const BigNumber* BigNumSrc2, BigNumber* output)
{
    BigNumber X,Y,Temp;
	int i,j,k;
	unsigned long n;
	unsigned long num;

	// 指数的二进制实际有效位
	k = BigNumSrc1->m_nLength*32 - 32; 		
	num = BigNumSrc1->m_ulValue[BigNumSrc1->m_nLength - 1];
	while(num)
	{
		num = num>>1;
		k++;
	}
	
	// X = N
	initByBigNum(&X,BigNumN);
	
	/* i初始化为k-2,其一:大数的二进制下标需要减1, 
	其二:大数的最高二进制位(即大数二进制的第k-1位)一定为1.
	i为什么初始化为k-2？？
	*/
	for(i = k-2;i >= 0;i--)
	{
		/*大数最高位运算 b=b*b Mod m;
		相当于：X = X*X Mod BigNumSrc2;
		*/ 
		Y = mulByInt(&X,X.m_ulValue[X.m_nLength-1]);
		Y = modByBigNum(&Y,BigNumSrc2);

		// 除去最高位(n=1)，再对大数的每一位，执行运算 b*b Mod m
        for(n = 1;n < X.m_nLength;n++)
		{   
			/*对Y，扩大32倍（大数为32进制）,也即增多1位。
			因为根据运算规则，从高位继承而来值(Y值)，在低位的时候需要扩展为(进制的倍数)32倍。
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
		// 执行赋值 b=b*b Mod m
		initByBigNum(&X,&Y);
		
		/* if条件判断：
			"i>>5" i表示二进制位的序号，右移5位(除以32)得到N,表示大数值的第N位(Dec).
			"i&31" 对取到的大数值第N位的数值，再取其二进制形式的第"i&31"位.
			判断是否为1，如果是，执行运算.
		*/
		if((BigNumSrc1->m_ulValue[i>>5]>>(i&31)) & 1)
		{
			/* 执行运算 a=a*b Mod m;
			相当于：X = BigNumN*X Mod BigNumSrc2;
			*/
			// 和BigNumN相乘？？
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
			// 又是存放在X ？？
		    initByBigNum(&X,&Y);
		}
	}
	
	initByBigNum(output, &X);
    return 0;
}





