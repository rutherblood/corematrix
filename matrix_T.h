//Implementation by Rutherblood
//Contact Author: <thecaptcha.314159@gmail.com> <github.com/rutherblood>
//*_T.* denotes code yet to be tested
//mxMatrix: Implements Addition, Subtracion, Scalar Mult, Matrix Mult, Transpose, Determinant (if m==n), Adjoint(if m==n), Inverse(if m==n)
// Functionality To Be Added: add exception handling for mem-allocation failures on mxData
//FTBA: Add exception for n==m==0
//FTBA: Add Co-factor finding function (DONE)
//FTBA: Optimize!

#include <math.h>

class mxMatrix
{
	private:
		float** mxData;
		int rows, columns;

	public:
		//Constructors
		mxMatrix(float** placeholderData,int rowsX,int columnsY) //placeholderData must be appropriately initialized
		{
			rows=rowsX; columns=columnsY;
			mxData=new float*[rows];
			for(int i=0;i<rows;i++)
			{
				mxData[i]=new float[columns];
				for(int j=0;j<columns;j++)
				{
					mxData[i][j]=placeholderData[i][j];
				}
			}

		}

		mxMatrix(int rowsX,int columnsY)
		{
			rows=rowsX; columns=columnsY;
			mxData=new float*[rows];
			for(int i=0;i<rows;i++)
			{
				mxData[i]=new float[columns];
			}
		}

		//Destructor
		~mxMatrix()
		{
			for(int i=0;i<rows;i++)
			{
				delete[] mxData[i];
			}
			delete mxData;
		}

		//Matrix-operation function prototypes: Defined in the matrix_T.cpp file
		//Matrix Addition
		static void mxAdd(mxMatrix,mxMatrix,mxMatrix);
		static mxMatrix mxAdd(mxMatrix,mxMatrix);
		void mxAdd(mxMatrix);

		//Scalar Multiplication
		static void mxScMult(mxMatrix,float,mxMatrix);
		static mxMatrix mxScMult(mxMatrix,float);
		void mxScMult(float);

		//Matrix Multiplication
		static void mxMxMult(mxMatrix,mxMatrix,mxMatrix);
		static mxMatrix mxMxMult(mxMatrix,mxMatrix);
		mxMatrix mxMxMult(mxMatrix);

		//Matrix Transpose
		static mxMatrix mxMxTranspose(mxMatrix);
		mxMatrix mxMxTranspose(void);

		//Co-factor
		static mxMatrix mxMxCoFact(mxMatrix,int, int);
		mxMatrix mxMxCoFact(int, int);

		//Determinant (if n==m)
		static float mxDeterminant(mxMatrix);
		float mxDeterminant(void);

		//Adjoint (if n==m)
		static mxMatrix mxAdj(mxMatrix);
		mxMatrix mxAdj(void);

		//Inverse (if n==m)
		static mxMatrix mxInv(mxMatrix);
		mxMatrix mxInv(void);
};
