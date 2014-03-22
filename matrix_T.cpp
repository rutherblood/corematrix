//Implementation by Rutherblood
//Contact Author: <thecaptcha.314159@gmail.com> <github.com/rutherblood>
//*_T.* denotes code yet to be fully-tested or under the testing phase
//mxMatrix: Implements Addition, Subtracion, Scalar Mult, Matrix Mult, Determinant (if m==n),Adjoint(if m==n), Inverse(if m==n)
//Functionality To Be Added (FTBA): add exception handling for mem-allocation failures on mxData
//FTBA: Add exception for n==m==0
//FTBA: Add Co-factor finding function (DONE)
//FTBA: Optimize!

#include "matrix_T.h"

//Matrix Addition

void mxMatrix::mxAdd(mxMatrix A,mxMatrix B,mxMatrix C) //C should be pre-initialized appropriately
{
		if((A.rows==B.rows==C.rows) && (A.columns==B.columns==C.columns))
		{
			for(int i=0;i<C.rows;i++)
			{
				for(int j=0;j<C.columns;j++)
					C.mxData[i][j]=A.mxData[i][j]+B.mxData[i][j];
			}
		}
}

mxMatrix mxMatrix::mxAdd(mxMatrix A,mxMatrix B)
{
	if((A.rows==B.rows) && (A.columns==B.columns))
	{
		mxMatrix C=mxMatrix(A.rows,A.columns);
		for(int i=0;i<C.rows;i++)
		{
			for(int j=0;j<C.columns;j++)
				C.mxData[i][j]=A.mxData[i][j]+B.mxData[i][j];
		}
		return C;
	}
}

void mxMatrix::mxAdd(mxMatrix B)
{
	if((rows==B.rows) && (columns==B.columns))
	{
		for(int i=0;i<rows;i++)
		{
			for(int j=0;j<columns;j++)
				mxData[i][j]+=B.mxData[i][j];
		}
	}
}


//Scalar Multiplication

void mxMatrix::mxScMult(mxMatrix A,float scalar,mxMatrix C) //C should be pre-initialized appropriately
{
	if((A.rows==C.rows) && (A.columns==C.columns))
	{
		for(int i=0;i<C.rows;i++)
		{
			for(int j=0;j<C.columns;j++)
				C.mxData[i][j]=A.mxData[i][j]*scalar;
		}
	}
}

mxMatrix mxMatrix::mxScMult(mxMatrix A,float scalar)
{
	mxMatrix C(A.rows,A.columns);
	for(int i=0;i<C.rows;i++)
	{
		for(int j=0;j<C.columns;j++)
			C.mxData[i][j]=A.mxData[i][j]*scalar;
	}
	return C;
}

void mxMatrix::mxScMult(float scalar)
{
	for(int i=0;i<rows;i++)
	{
		for(int j=0;j<columns;j++)
			mxData[i][j]*=scalar;
	}
}


//Matrix Multiplication

void mxMatrix::mxMxMult(mxMatrix A,mxMatrix B,mxMatrix C)  //C should be pre-initialized appropriately
{
	if(A.columns==B.rows && A.rows==C.rows && B.columns==C.columns)
	{
		float temVaria;
		for(int i=0;i<C.rows;i++)
		{
			for(int j=0;j<C.columns;j++)
			{
				for(int k=0,temVaria=0;k<A.columns;k++)
					temVaria+=A.mxData[i][k]*B.mxData[k][j];

				C.mxData[i][j]=temVaria;
			}
		}
	}
}

mxMatrix mxMatrix::mxMxMult(mxMatrix A,mxMatrix B)
{
	if(A.columns==B.rows)
	{
		mxMatrix C(A.rows,B.columns);
		float temVaria;
		for(int i=0;i<C.rows;i++)
		{
			for(int j=0;j<C.columns;j++)
			{
				for(int k=0,temVaria=0;k<A.columns;k++)
					temVaria+=A.mxData[i][k]*B.mxData[k][j];

				C.mxData[i][j]=temVaria;
			}
		}
		return C;
	}
}

mxMatrix mxMatrix::mxMxMult(mxMatrix B)
{
	if(columns==B.rows)
	{
		mxMatrix C(rows,B.columns);
		float temVaria;
		for(int i=0;i<C.rows;i++)
		{
			for(int j=0;j<C.columns;j++)
			{
				for(int k=0,temVaria=0;k<columns;k++)
					temVaria+=mxData[i][k]*B.mxData[k][j];

				C.mxData[i][j]=temVaria;
			}
		}
		return C;
	}
}


//Matrix Transpose
mxMatrix mxMatrix::mxMxTranspose(mxMatrix A)
{
	mxMatrix C(A.columns,A.rows);
	for(int i=0;i<C.rows;i++)
	{
		for(int j=0;j<C.columns;j++)
		{
			C.mxData[i][j]=A.mxData[j][i];
		}
	}
	return C;
}

mxMatrix mxMatrix::mxMxTranspose(void)
{
	mxMatrix C(columns,rows);
	for(int i=0;i<C.rows;i++)
	{
		for(int j=0;j<C.columns;j++)
		{
			C.mxData[i][j]=mxData[j][i];
		}
	}
	return C;
}



//Co-Factor
mxMatrix mxMatrix::mxMxCoFact(mxMatrix A,int aRow, int bColumn)
{
	if(aRow<A.rows && bColumn<A.columns)
	{
		mxMatrix C(A.rows-1,A.columns-1);
		for(int i=0,s=0;i<A.rows;i++)
		{
			if(i==aRow) continue;
			for(int j=0,t=0;j<A.columns;j++)
			{
				if(j==bColumn) continue;
				C.mxData[s][t]=A.mxData[i][j];
				t+=1;
			}
			s+=1;
		}
		return C;
	}

	else return mxMatrix(0,0);
}

mxMatrix mxMatrix::mxMxCoFact(int aRow, int bColumn)
{
	if(aRow<rows && bColumn<columns)
	{
		mxMatrix C(rows-1,columns-1);
		for(int i=0,s=0;i<rows;i++)
		{
			if(i==aRow) continue;
			for(int j=0,t=0;j<columns;j++)
			{
				if(j==bColumn) continue;
				C.mxData[s][t]=mxData[i][j];
				t+=1;
			}
			s+=1;
		}
		return C;
	}

	else return mxMatrix(0,0);
}



//Determinant (if n==m)
//A[1][1]*determinant(cofactor(A,1,1))-A[1][2]*determinant(cofactor(A,2,2))+....

float mxMatrix::mxDeterminant(mxMatrix A)
{
	if(A.rows!=A.columns)
	{
		return 0;
	}
	else if(A.rows==A.columns==1)
	{
		return A.mxData[0][0];  //Fixed 1 1 bug
	}
	else
	{
		float summage=0;
		for(int j=0;j<A.columns;j++)
		{
			mxMatrix temVaria=mxMxCoFact(A,0,j);
			summage+=(pow(-1,j)*A.mxData[0][j]*mxDeterminant(temVaria));
		}
		return summage;
	}
}

float mxMatrix::mxDeterminant()
{
	if(rows!=columns)
	{
		return 0;
	}
	else if(rows==columns==1)
	{
		return mxData[0][0]; //Fixed 1 1 bug
	}
	else
	{
		float summage=0;
		for(int j=0;j<columns;j++)
		{
			mxMatrix temVaria=mxMxCoFact(*this,0,j);
			summage+=(pow(-1,j)*mxData[0][j]*mxDeterminant(temVaria));
		}
		return summage;
	}
}


//Adjoint

mxMatrix mxMatrix::mxAdj(mxMatrix A)
{
	if(A.rows!=A.columns)
	{
		return mxMatrix(0,0);
	}
	else
	{
		mxMatrix C(A.rows,A.columns);
		for(int i=0;i<A.rows;i++)
		{
			for(int j=0;j<A.columns;j++)
			{
				mxMatrix temVaria=mxMxCoFact(A,i,j);
				C.mxData[i][j]=pow(-1,(i+j))*mxDeterminant(temVaria);
			}
		}
		return mxMxTranspose(C);                                              //Transpose bug fixed
	}
}

mxMatrix mxMatrix::mxAdj()
{
	if(rows!=columns)
	{
		return mxMatrix(0,0);
	}
	else
	{
		mxMatrix C(rows,columns);
		for(int i=0;i<rows;i++)
		{
			for(int j=0;j<columns;j++)
			{
				mxMatrix temVaria=mxMxCoFact(*this,i,j);
				C.mxData[i][j]=pow(-1,(i+j))*mxDeterminant(temVaria);
			}
		}
		return C.mxMxTranspose();                                           //Transpose bug fixed
	}
}


//Matrix Inversing

mxMatrix mxMatrix::mxInv(mxMatrix A)
{
	float invDetA=A.mxDeterminant();
	invDetA=1.0/invDetA;
	mxMatrix invA=A.mxAdj();
	invA.mxScMult(invDetA);
	return invA;
}

mxMatrix mxMatrix::mxInv(void)
{
	float invDetA=mxDeterminant();
	invDetA=1.0/invDetA;
	mxMatrix invA=mxAdj();
	invA.mxScMult(invDetA);
	return invA;
}


