
//Test
use matrix::mxMatrix;
use matrix::mxTrait;
//use super;
mod matrix;

fn main()
{
    let mut C=mxMatrix::new(2,2);C.initData(~[~[5.0,6.0],~[7.0,8.0]]);
    let mut J=mxMatrix::new(2,2);J.initData(~[~[0.0,55.0],~[71.0,3.0]]);
    //let mut f=mxMatrix::mxAdj(&C); //print!("{} ",C.columns);
    //let f=mxMatrix::mxAdd(&C,&J);
    C.mxAdd(&J);
    displayMatrix(&C);
    //f=mxMatrix::
    //print!("{} ",f);
}

fn displayMatrix(temp:&mxMatrix)
{
    for i in range(0,2)
    {
        for j in range(0,2)
        {
            print!("{} ",temp.mxData[i][j]);
        }
        println!("");
    }
}