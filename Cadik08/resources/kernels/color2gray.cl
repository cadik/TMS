#ifdef cl_amd_fp64
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#endif

#ifdef cl_amd_printf
#pragma OPENCL EXTENSION cl_amd_printf : enable
#endif

//______________________________________________________________________________
__kernel void calc_error(__global const double* const grad,
                         __global double* const err,
                         const uint rows, const uint cols)
{
	const uint2 gid = {get_global_id(1), get_global_id(0)};
	if (gid.y >= rows || gid.x >= cols)
		return;

	const uint i = gid.y * cols + gid.x;
	err[i] = grad[3 * i] - grad[3 * i + 1] +
	         grad[3 * (i + 1) + 1] - grad[3 * (i + cols)];
}

//______________________________________________________________________________
__kernel void correct_grad(__global double* const grad,
                           __global const double* const err, const double s,
                           const uint rows, const uint cols, const uint n)
{
	const uint2 gid = {get_global_id(1), get_global_id(0)};
	if (gid.y >= rows || gid.x >= cols)
		return;

	const uint i = gid.y * cols + gid.x;
	const double e = .25 * s * err[i];
	if (n == 0)
		grad[3 * i] -= e;
	else if (n == 1)
		grad[3 * i + 1] += e;
	else if (n == 2) {
		if (gid.y < rows && gid.x + 1 < cols)
			grad[3 * (i + 1)] -= e;
	}
	else if (n == 3) {
		if ((gid.y + 1 < rows && gid.x < cols))
			grad[3 * (i + cols)] += e;
	}
}

//______________________________________________________________________________
__kernel void reduce_min(__global const double* const in,
                         __local double* const tmp,
                         __global double* const minimas, const uint n)
{
	const size_t gid = get_global_id(0),
	             lid = get_local_id(0),
	             wgs = get_local_size(0),
	             wgid = get_group_id(0);

	const uint stride = wgid * wgs;
	tmp[lid] = stride + gid < n ? in[stride + gid] : 0.;
	tmp[lid + wgs] = stride + gid + wgs < n ? in[stride + gid + wgs] :
	                                          0.;

	for (ushort offset = wgs; offset; offset >>= 1) {
		barrier(CLK_LOCAL_MEM_FENCE);
		if (lid < offset)
			tmp[lid] = min(tmp[lid], tmp[lid + offset]);
	}

	barrier(CLK_LOCAL_MEM_FENCE);
	if (!lid)
		minimas[wgid] = tmp[0];
}

//______________________________________________________________________________
__kernel void reduce_absmax(__global const double* const in,
                            __local double* const tmp,
                            __global double* const maximas, const uint n)
{
	const size_t gid = get_global_id(0),
	             lid = get_local_id(0),
	             wgs = get_local_size(0),
	             wgid = get_group_id(0);

	const uint stride = wgid * wgs;
	tmp[lid] = stride + gid < n ? fabs(in[stride + gid]) : 0.;
	tmp[lid + wgs] = stride + gid + wgs < n ? fabs(in[stride + gid + wgs]) :
	                                          0.;

	for (ushort offset = wgs; offset; offset >>= 1) {
		barrier(CLK_LOCAL_MEM_FENCE);
		if (lid < offset)
			tmp[lid] = max(tmp[lid], tmp[lid + offset]);
	}

	barrier(CLK_LOCAL_MEM_FENCE);
	if (!lid)
		maximas[wgid] = tmp[0];
}

//______________________________________________________________________________
__kernel void reduce_max(__global const double* const in,
                         __local double* const tmp,
                         __global double* const maximas, const uint n)
{
	const size_t gid = get_global_id(0),
	             lid = get_local_id(0),
	             wgs = get_local_size(0),
	             wgid = get_group_id(0);

	const uint stride = wgid * wgs;
	tmp[lid] = stride + gid < n ? in[stride + gid] : 0.;
	tmp[lid + wgs] = stride + gid + wgs < n ? in[stride + gid + wgs] :
	                                          0.;

	for (ushort offset = wgs; offset; offset >>= 1) {
		barrier(CLK_LOCAL_MEM_FENCE);
		if (lid < offset)
			tmp[lid] = max(tmp[lid], tmp[lid + offset]);
	}

	barrier(CLK_LOCAL_MEM_FENCE);
	if (!lid)
		maximas[wgid] = tmp[0];
}

//______________________________________________________________________________
__kernel void integrate2x(__global double* grad, __global double* out,
                          const uint rows, const uint cols)
{
	const size_t y = get_global_id(0);
	if (y >= rows)
		return;

	const size_t i = y * cols;
	for (size_t x = 1; x < cols; ++x) {
		const size_t j = i + x;
		out[3 * j] =
		out[3 * j + 1] =
		out[3 * j + 2] = out[3 * (j - 1)] +
		                 grad[3 * (j - 1)];
	}
}

//______________________________________________________________________________
__kernel void trans_range(__global double* const in,
                          const double old_min, const double old_max,
                          const double new_min, const double new_max,
                          const uint n)
{
	const size_t gid = get_global_id(0);

	if (gid >= n)
		return;

	const size_t i = 3 * gid;
	in[i] =
	in[i + 1] =
	in[i + 2] = ((in[i] - old_min) * (new_max - new_min)) /
	            (old_max - old_min) + new_min;
}

// MAXIMUM ERROR SELECTION code: rather slow
//______________________________________________________________________________
/*__kernel void correct_gradi(__global double* const grad,
                            __global const double* const err,
                            const double s, const uint uv,
                            const uint rows, const uint cols)
{
	const uint2 id = {get_global_id(1), get_global_id(0)};
	if (id.y >= rows || id.x >= cols)
		return;
	const uint i = id.y * cols + id.x;

	if (i == uv) {
		const double e = err[i] * .25 * s;

		grad[3 * i] -= e;
		grad[3 * i + 1] += e;
		if (id.x + 1 < cols)
			grad[3 * (i + 1) + 1] -= e;
		if (id.y + 1 < rows)
			grad[3 * (i + cols)] += e;
	}
}

//______________________________________________________________________________
__kernel void reduce_maxi(__global const double* const in,
                          __global const uint* const uv,
                          const uint use_uv,
                          __local double* const tmp,
                          __local uint* const inds,
                          __global double* const maximas,
                          __global uint* const coords,
                          const uint n)
{
	const size_t gid = get_global_id(0),
	             lid = get_local_id(0),
	             wgs = get_local_size(0),
	             wgid = get_group_id(0);

	const uint stride = wgid * wgs;
	tmp[lid] = stride + gid < n ? fabs(in[stride + gid]) : 0.;
	tmp[lid + wgs] = stride + gid + wgs < n ? fabs(in[stride + gid + wgs]) :
	                                          0.;
	if (use_uv) {
		inds[lid] = stride + gid < n ? uv[stride + gid] : 0;
		inds[lid + wgs] = stride + gid + wgs < n ? uv[stride + gid + wgs] :
		                                           0;
	}
	else {
		inds[lid] = stride + gid;
		inds[lid + wgs] = stride + gid + wgs;
	}

	for (ushort offset = wgs; offset; offset >>= 1) {
		barrier(CLK_LOCAL_MEM_FENCE);
		if (lid < offset)
			if (tmp[lid] < tmp[lid + offset]) {
				tmp[lid] = tmp[lid + offset];
				inds[lid] = inds[lid + offset];
			}
	}

	barrier(CLK_LOCAL_MEM_FENCE);
	if (!lid) {
		maximas[wgid] = tmp[0];
		coords[wgid] = inds[0];
	}
}*/


// INTERLEAVED GRADIENT CORRECTION
//______________________________________________________________________________
/*__kernel void correct_grad(__global double* const grad,
                           __local double* const tmp,
                           __global double* const err,
                           const double s,
                           const uint rows, const uint cols)
{
	const uint2 gid = {get_global_id(1), get_global_id(0)},
	            lid = {get_local_id(1), get_local_id(0)},
	            wgs = {get_local_size(1), get_local_size(0)};

	const bool inrange = gid.y < rows && gid.x < cols;
	const uint i = gid.y * cols + gid.x,
	           j = lid.y * wgs.x + lid.x;

	tmp[2 * j] = inrange ? grad[3 * i] : 0.;
	tmp[2 * j + 1] = inrange ? grad[3 * i + 1] : 0.;
	if (lid.x == wgs.x - 1) {
		tmp[2 * (j + 1)] = (gid.y < rows &&
		                    gid.x + 1 < cols) ? grad[3 * (i + 1)] : 0.;
		tmp[2 * (j + 1) + 1] = (gid.y < rows &&
		                        gid.x + 1 < cols) ? grad[3 * (i + 1) + 1] : 0.;
	}
	if (lid.y == wgs.y - 1) {
		tmp[2 * (j + wgs.x)] = (gid.y + 1 < rows &&
		                        gid.x < cols) ? grad[3 * (i + cols)] : 0.;
		tmp[2 * (j + wgs.x) + 1] = (gid.y + 1 < rows &&
		                            gid.x < cols) ? grad[3 * (i + cols) + 1] : 0.;
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	if (!(lid.y % 2) && !(lid.x % 2)) {
		double e = tmp[2 * j] - tmp[2 * j + 1] +
		           tmp[2 * (j + 1) + 1] - tmp[2 * (j + wgs.x)];
		e *= .25 * s;

		tmp[2 * j] -= e;
		tmp[2 * j + 1] += e;
		tmp[2 * (j + 1) + 1] -= e;
		tmp[2 * (j + wgs.x)] += e;

		if (inrange)
			err[i] = e;
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	if ((lid.y % 2 && !(lid.x % 2)) || (lid.x % 2 && !(lid.y % 2))) {
		double e = tmp[2 * j] - tmp[2 * j + 1] +
		           tmp[2 * (j + 1) + 1] - tmp[2 * (j + wgs.x)];
		e *= .25 * s;

		tmp[2 * j] -= e;
		tmp[2 * j + 1] += e;
		tmp[2 * (j + 1) + 1] -= e;
		tmp[2 * (j + wgs.x)] += e;

		if (inrange)
			err[i] = e;
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	if (lid.y % 2 && lid.x % 2) {
		double e = tmp[2 * j] - tmp[2 * j + 1] +
		           tmp[2 * (j + 1) + 1] - tmp[2 * (j + wgs.x)];
		e *= .25 * s;

		tmp[2 * j] -= e;
		tmp[2 * j + 1] += e;
		tmp[2 * (j + 1) + 1] -= e;
		tmp[2 * (j + wgs.x)] += e;

		if (inrange)
			err[i] = e;
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	if (inrange) {
		grad[3 * i] = tmp[2 * j];
		grad[3 * i + 1] = tmp[2 * j + 1];
	}
}*/
