//Implementation by Rutherblood
//Contact Author: <thecaptcha.314159@gmail.com> <github.com/rutherblood>
//*_T.* denotes code yet to be tested
//mxMatrix: Implements Addition, Subtracion, Scalar Mult, Matrix Mult, Transpose, Determinant (if m==n), Adjoint(if m==n), Inverse(if m==n)
// Functionality To Be Added: add exception handling for mem-allocation failures on mxData
// Add exception & error handling
//FTBA: Optimize!



	pub fn pow(base:f32,exp:int)->f32
	{
	    let mut temVaria=1.0;
	    if exp==0
	    {
	        return 1.0;
	    }
	    else
	    {
	        for i in range(0,exp)
	        {
	            temVaria*=base;
	        }
	        return temVaria;
	    }
	}


	pub struct mxMatrix{
	    rows:int,
	    columns:int,
	    mxData:~[~[f32]]
	}

	//static fns
	//  fn mxAdd(A:&mxMatrix,B:&mxMatrix)->mxMatrix;
	//  fn mxScMult(A:&mxMatrix,scalar:f32)->mxMatrix;
	//  fn mxMxMult(A:&mxMatrix,B:&mxMatrix)->mxMatrix;
	//  fn mxMxTranspose(A:&mxMatrix)->mxMatrix;
	//  fn mxMxCoFact(A:&mxMatrix,aRow:int,bColumn:int)->mxMatrix;
	//  fn mxDeterminant(A:&mxMatrix)->f32;
	//  fn mxAdj(A:&mxMatrix)->mxMatrix;
	//  fn mxInv(A:&mxMatrix)->mxMatrix;


	impl mxMatrix{
	    pub fn new(rowsX:int,columnsY:int)->mxMatrix
	    {
	        let mut temp=mxMatrix{rows:rowsX,columns:columnsY,mxData:~[]};
	        for i in range(0,rowsX)
	        {
	            temp.mxData.push(~[]);
	            for j in range(0,columnsY)
	            {
	                temp.mxData[i].push(0.0);
	            }
	        }
	        return temp;
	    }

	    pub fn initData(&mut self,dataArg:~[~[f32]])
	    {
	        self.mxData=dataArg;
	    }

	    pub fn mxAdd(A:&mxMatrix,B:&mxMatrix)->mxMatrix
	    {
	        if (A.rows==B.rows) && (A.columns==B.columns)
	        {
	            let mut C=mxMatrix::new(A.rows,A.columns);
	            for i in range(0,C.rows)
	            {
	                for j in range(0,C.columns)
	                {
	                    C.mxData[i][j]=(A.mxData[i][j]+B.mxData[i][j]);
	                }
	            }
	            return C;
	        }
	        else{return mxMatrix::new(0,0);}
	    }

	    pub fn mxScMult(A:&mxMatrix,scalar:f32)->mxMatrix
	    {
	        let mut C=mxMatrix::new(A.rows,A.columns);
	        for i in range(0,C.rows)
	        {
	            for j in range(0,C.columns)
	            {
	                C.mxData[i][j]=A.mxData[i][j]*scalar;
	            }
	        }
	        return C;
	    }

	    pub fn mxMxMult(A:&mxMatrix,B:&mxMatrix)->mxMatrix
	    {
	        if A.columns==B.rows
	        {
	            let mut C=mxMatrix::new(A.rows,B.columns);
	            let mut temVaria=0.0;
	            for i in range(0,C.rows)
	            {
	                for j in range(0,C.columns)
	                {
	                    temVaria=0.0;
	                    for k in range(0,A.columns)
	                    {
	                        temVaria+=A.mxData[i][k]*B.mxData[k][j];
	                    }
	                    C.mxData[i][j]=temVaria;
	                }
	            }
	            return C;
	        }
	        else{return mxMatrix::new(0,0);}
	    }

	    pub fn mxMxTranspose(A:&mxMatrix)->mxMatrix
	    {
	        let mut C=mxMatrix::new(A.columns,A.rows);
	        for i in range(0,C.rows)
	        {
	            for j in range(0,C.columns)
	            {
	                C.mxData[i][j]=A.mxData[j][i];
	            }
	        }
	        return C;
	    }

	    pub fn mxMxCoFact(A:&mxMatrix,aRow:int,bColumn:int)->mxMatrix
	    {
	        if aRow<A.rows && bColumn<A.columns
	        {
	            let mut C=mxMatrix::new(A.rows-1,A.columns-1);
	            let mut s=0;
	            for i in range(0,A.rows)
	            {
	                if i==aRow{continue;}
	                let mut t=0;
	                for j in range(0,A.columns)
	                {
	                    if j==bColumn{continue;}
	                    C.mxData[s][t]=A.mxData[i][j];
	                    t+=1;
	                }
	                s+=1;
	            }
	            return C;
	        }
	        else{return mxMatrix::new(0,0);}
	    }

	    pub fn mxDeterminant(A:&mxMatrix)->f32
	    {
	        //let temVaria=mxMatrix::mxMxCoFact(A,0,0);
	        if A.rows!=A.columns
	        {
	            return 0.0;
	        }
	        else if A.rows==1 && A.columns==1
	        {
	            return A.mxData[0][0];
	        }
	        else
	        {
	            let mut summage=0.0;
	            for j in range(0,A.columns)
	            {
	                let temVaria=mxMatrix::mxMxCoFact(A,0,j);
	                summage+=A.mxData[0][j]*mxMatrix::mxDeterminant(&temVaria)*pow(-1.0,j);
	            }
	            return summage;
	        }
	    }

	    pub fn mxAdj(A:&mxMatrix)->mxMatrix
	    {
	        if A.rows!=A.columns
	        {
	            return mxMatrix::new(0,0);
	        }
	        else
	        {
	            let mut C=mxMatrix::new(A.rows,A.columns);
	            for i in range(0,C.rows)
	            {
	                for j in range(0,C.columns)
	                {
	                    let temVaria=mxMatrix::mxMxCoFact(A,i,j);
	                    C.mxData[i][j]=pow(-1.0,(i+j))*mxMatrix::mxDeterminant(&temVaria);
	                }
	            }
	            return mxMatrix::mxMxTranspose(&C);
	        }
	    }

	    pub fn mxInv(A:&mxMatrix)->mxMatrix
	    {
	        let mut invDetA=mxMatrix::mxDeterminant(A);
	        invDetA=1.0/invDetA;
	        let mut invA=mxMatrix::mxAdj(A);
	        invA=mxMatrix::mxScMult(&invA,invDetA);
	        return invA;
	    }
	}

	pub trait mxTrait{
	    fn mxAdd(&mut self,B:&mxMatrix);
	    fn mxScMult(&mut self,scalar:f32);
	    fn mxMxMult(&self,B:&mxMatrix)->mxMatrix;
	    fn mxMxTranspose(&self)->mxMatrix;
	    fn mxMxCoFact(&self,aRow:int,bColumn:int)->mxMatrix;
	    fn mxDeterminant(&self)->f32;
	    fn mxAdj(&self)->mxMatrix;
	    fn mxInv(&self)->mxMatrix;
	}


	impl mxTrait for mxMatrix
	{
	    fn mxAdd(&mut self,B:&mxMatrix)
	    {
	        if (self.rows==B.rows) && (self.columns==B.columns)
	        {
	            for i in range(0,self.rows)
	            {
	                for j in range(0,self.columns)
	                {
	                    self.mxData[i][j]+=B.mxData[i][j];
	                }
	            }
	        }
	    }

	    fn mxScMult(&mut self,scalar:f32)
	    {
	        for i in range(0,self.rows)
	        {
	            for j in range(0,self.columns)
	            {
	                self.mxData[i][j]*=scalar;
	            }
	        }
	    }

	    fn mxMxMult(&self,B:&mxMatrix)->mxMatrix
	    {
	        if self.columns==B.rows
	        {
	            let mut C=mxMatrix::new(self.rows,B.columns);
	            let mut temVaria=0.0;
	            for i in range(0,C.rows)
	            {
	                for j in range(0,C.columns)
	                {
	                    temVaria=0.0;
	                    for k in range(0,self.columns)
	                    {
	                        temVaria+=self.mxData[i][k]*B.mxData[k][j];
	                    }
	                    C.mxData[i][j]=temVaria;
	                }
	            }
	            return C;
	        }
	        else{return mxMatrix::new(0,0);}
	    }

	    fn mxMxTranspose(&self)->mxMatrix
	    {
	        let mut C=mxMatrix::new(self.columns,self.rows);
	        for i in range(0,C.rows)
	        {
	            for j in range(0,C.columns)
	            {
	                C.mxData[i][j]=self.mxData[j][i];
	            }
	        }
	        return C;
	    }

	    fn mxMxCoFact(&self,aRow:int,bColumn:int)->mxMatrix
	    {
	        if aRow<self.rows && bColumn<self.columns
	        {
	            let mut C=mxMatrix::new(self.rows-1,self.columns-1);
	            let mut s=0;
	            for i in range(0,self.rows)
	            {
	                if i==aRow{continue;}
	                let mut t=0;
	                for j in range(0,self.columns)
	                {
	                    if j==bColumn{continue;}
	                    C.mxData[s][t]=self.mxData[i][j];
	                    t+=1;
	                }
	                s+=1;
	            }
	            return C;
	        }
	        else{return mxMatrix::new(0,0);}
	    }

	    fn mxDeterminant(&self)->f32
	    {
	        if self.rows!=self.columns
	        {
	            return 0.0;
	        }
	        else if self.rows==1 && self.columns==1
	        {
	            return self.mxData[0][0];
	        }
	        else
	        {
	            let mut summage=0.0;
	            for j in range(0,self.columns)
	            {
	                let temVaria=self.mxMxCoFact(0,j);
	                summage+=self.mxData[0][j]*temVaria.mxDeterminant()*pow(-1.0,j);
	            }
	            return summage;
	        }
	    }

	    fn mxAdj(&self)->mxMatrix
	    {
	        if self.rows!=self.columns
	        {
	            return mxMatrix::new(0,0);
	        }
	        else
	        {
	            let mut C=mxMatrix::new(self.rows,self.columns);
	            for i in range(0,C.rows)
	            {
	                for j in range(0,C.columns)
	                {
	                    let temVaria=self.mxMxCoFact(i,j);
	                    C.mxData[i][j]=pow(-1.0,(i+j))*temVaria.mxDeterminant();
	                }
	            }
	            return C.mxMxTranspose();
	        }
	    }

	    fn mxInv(&self)->mxMatrix
	    {
	        let mut invDetA=self.mxDeterminant();
	        invDetA=1.0/invDetA;
	        let mut invA=self.mxAdj();
	        invA.mxScMult(invDetA);
	        return invA;
	    }
	}