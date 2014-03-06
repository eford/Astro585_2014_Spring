extern "C"   // ensure function name to be left alone
{

    __global__ void normal_pdf_gpu(const double *x, double *y, unsigned int n)
    {
	// assumes a 2-d grid of 1-d blocks
	unsigned int i = (blockIdx.y * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;
        if(i<n)  y[i] = exp(-0.5*x[i]*x[i])*rsqrt(2.0*M_PI);
    }

}


