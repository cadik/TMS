#ifdef cl_khr_fp64
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#else
#ifdef cl_amd_fp64
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#endif
#endif

//______________________________________________________________________________
__kernel void calc_error(__global double* grad, __global double* err,
                         const uint rows, const uint cols)
{
	const uint2 id = {get_global_id(1), get_global_id(0)};
	if (id.y >= rows || id.x >= cols)
		return;
	const uint i = id.y * cols + id.x;

	err[i] = grad[3 * i] - grad[3 * i + 1] +
	         (id.x + 1 < cols) ? grad[3 * (i + 1) + 1] : 0. -
	         (id.y + 1 < rows) ? grad[3 * (i + cols)] : 0.;
}

//______________________________________________________________________________
__kernel void correct_grad(__global double* grad, __global double* err,
                           const double s, const uchar mode,
                           const uint rows, const uint cols)
{
	const uint2 id = {get_global_id(1), get_global_id(0)};
	if (id.y >= rows || id.x >= cols)
		return;
	const uint i = id.y * cols + id.x;

	const double e = err[i] * .25 * s;

	switch (mode) {
		case 0:
		grad[3 * i] -= e;
		break;
		case 1:
		grad[3 * i + 1] += e;
		break;
		case 2:
		if (id.x + 1 < cols)
			grad[3 * (i + 1) + 1] -= e;
		break;
		case 3:
		if (id.y + 1 < rows)
			grad[3 * (i + cols)] += e;
		default:
	}
}

//______________________________________________________________________________
__kernel void reduce(__global double* in, __local double* tmp,
                     __global double* maximas, const uint n)
{
	const size_t gid = get_global_id(0),
	             lid = get_local_id(0),
	             wgs = get_local_size(0);

	tmp[lid] = gid < n ? fabs(in[gid]) : 0.;
	tmp[lid + wgs] = gid + wgs < n ? fabs(in[gid + wgs]) : 0.;

	for (ushort offset = wgs; offset; offset >>= 1) {
		barrier(CLK_LOCAL_MEM_FENCE);
		if (lid < offset)
			tmp[lid] = max(tmp[lid], tmp[lid + offset]);
	}

	barrier(CLK_LOCAL_MEM_FENCE);
	if (!lid)
		maximas[get_group_id(0)] = tmp[0];
}
/*
		if (i > 0)
			pDst_image[3 * tmp_y] = pDst_image[3 * tmp_y + 1] =
			pDst_image[3*tmp_y+2] = (pDst_image[3 * (tmp_y - xmax)] +
			                        pG_image[3 * (tmp_y - xmax) + 1]);
			//neboli: OUTPUT_BW[0][y] = OUTPUT_BW[0][y-1] + Grad_Y[0][y-1];

		for (j = 1; j < xmax; ++j) {
			tmp_ind = j + tmp_y;
			pDst_image[3 * tmp_ind] = pDst_image[3 * tmp_ind + 1] =
			pDst_image[3 * tmp_ind + 2] = (pDst_image[3 * (tmp_ind - 1)] +
			                              pG_image[3 * (tmp_ind - 1)]);
			//neboli: OUTPUT_BW[x][y] = OUTPUT_BW[x-1][y] + Grad_X[x-1][y]; 
		}
*/
//______________________________________________________________________________
__kernel void integrate(__global double* grad, __global double* out,
                        const uint rows, const uint cols)
{
	const size_t y = get_global_id(0);
	if (y >= rows)
		return;

	if (!y)
		out[0] = out[1] = out[2] = 0.;

	const size_t i = y * cols;
	if (y > 1)
		out[3 * i] = out[3 * i + 1] = out[3 * i + 2] =
		out[3 * (i - cols)] + grad[3 * (i - cols) + 1];

	for (size_t x = 1; x < cols; ++x) {
		const size_t j = i + x;
		out[3 * j] = out[3 * j + 1] = out[3 * j + 2] =
		out[3 * (j - 1)] + grad[3 * (j - 1)];
	}
}
