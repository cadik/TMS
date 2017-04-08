#ifdef cl_amd_fp64
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#endif

#include "morton.cl"

//______________________________________________________________________________
// chessboard version
/*enum { CHESS_EVEN, CHESS_ODD, CHESS_DIAG };

__kernel void correct_grad_chess(__global double* const g,
                                 __global double* const err,
                                 const double s,
                                 const uint rows, const uint cols,
                                 const uchar mode)
{
	const uint2 gid = {get_global_id(1), get_global_id(0)};
	if (gid.y >= rows || gid.x >= cols)
		return;
	const uint i = gid.y * cols + gid.x;

	if ((mode == CHESS_EVEN && !(gid.y % 2) && !(gid.x % 2)) ||
	    (mode == CHESS_ODD && (gid.y % 2 && !(gid.x % 2)) ||
	     (gid.x % 2 && !(gid.y % 2)) ||
	    (mode == CHESS_DIAG && gid.y % 2 && gid.x % 2))) {
		const double e = g[3 * i] - g[3 * i + 1] +
		                 ((gid.x + 1 < cols) ? g[3 * (i + 1) + 1] :
		                                       0.) -
		                 ((gid.y + 1 < rows) ? g[3 * (i + cols)] :
		                                       0.);
		const double r = .25 * s * e;

		g[3 * i] -= r;
		g[3 * i + 1] += r;
		if (gid.x + 1 < cols)
			g[3 * (i + 1) + 1] -= r;
		if (gid.y + 1 < rows)
			g[3 * (i + cols)] += r;

		err[i] = e;
	}
}*/

//______________________________________________________________________________
// used to evaluate coarser quadtree levels
__kernel void avg_grad(__global double2* const rooti,
                       const uint aligni, const uint sizei)
{
	const size_t gid = get_global_id(0);
	if (gid >= sizei)
		return;

	rooti[aligni + gid] = (rooti[aligni + sizei + 4 * gid] +
	                       rooti[aligni + sizei + 4 * gid + 1] +
	                       rooti[aligni + sizei + 4 * gid + 2] +
	                       rooti[aligni + sizei + 4 * gid + 3]) /
	                      4;
}

//______________________________________________________________________________
__kernel void tag_active(__global const uchar* const flag,
                         __global const uint* const offsets,
                         __global const morton* const zs,
                         __global morton* const parent_index,
                         const uint n_a_i)
{
	const size_t gid = get_global_id(0);
	if (gid < n_a_i && flag[gid])
		parent_index[offsets[gid]] = zs[gid];
}

//______________________________________________________________________________
// for all active loops in a quadtree level i = n_a_i
__kernel void calc_err(__global const double2* const root, const uint align,
                       const uint size,
                       __global const morton* const parent_index,
                       __global morton* const zs,
                       __global uchar* const flag, const double eps,
                       const uint n_a_i)
{
	const size_t gid = get_global_id(0);
	if (gid >= n_a_i)
		return;

	const uint side = convert_uint(sqrt(convert_float(size)));
	const morton z = (parent_index[gid / 4] << 2) | (gid % 4);
	const uint2 pos = morton2d_decode(z);

	if (z >= size) {
		printf(" [%d, %lu] ", size, z);
		flag[gid] = 0;
		zs[gid] = 0;
		return;
	}

	const double e = root[align + z].x - root[align + z].y +
	                 (pos.x + 1 < side ? root[align + morton2d_encode(pos.x + 1, pos.y)].y : 0.) -
	                 (pos.y + 1 < side ? root[align + morton2d_encode(pos.x, pos.y + 1)].x : 0.);

	flag[gid] = fabs(e) > eps ? 1 : 0;
	zs[gid] = z;
}

enum { CHESS_EVEN, CHESS_ODD, CHESS_DIAG };

//______________________________________________________________________________
__kernel void correct_grad(__global const double2* const root0,
                           const uint align0,
                           __global const morton* const parent_index,
                           __global double* const err, const double eps,
                           const double s, const uint size, const uchar mode,
                           const uint n_a_i)
{
	const size_t gid = get_global_id(0);
	if (gid >= n_a_i)
		return;

	const uint side = convert_uint(sqrt(convert_float(size)));
	const morton z = (parent_index[gid / 4] << 2) | (gid % 4);
	const uint2 pos = morton2d_decode(z);

	if (pos.x >= side || pos.y >= side) {
		printf("%lu ", z);
		return;
	}

	if ((mode == CHESS_EVEN && !(pos.y % 2) && !(pos.x % 2)) ||
	    (mode == CHESS_ODD && (pos.y % 2 && !(pos.x % 2)) ||
	     (pos.x % 2 && !(pos.y % 2)) ||
	    (mode == CHESS_DIAG && pos.y % 2 && pos.x % 2))) {
		const morton z1 = morton2d_encode(pos.x + 1, pos.y),
		             z2 = morton2d_encode(pos.x, pos.y + 1);
		const double e = root0[align0 + z].x - root0[align0 + z].y +
		                 (pos.x + 1 < side ? root0[align0 + z1].y : 0.) -
		                 (pos.y + 1 < side ? root0[align0 + z2].x : 0.);

		const double r = .25 * s * e;
		root0[align0 + z].x -= r;
		root0[align0 + z].y += r;
		if (pos.x + 1 < side)
			root0[align0 + z1].y -= r;
		if (pos.y + 1 < side)
			root0[align0 + z2].x += r;

		err[gid] = e;
	}
}

//______________________________________________________________________________
/*__kernel void correct_grad(__global double* const grad,
                           const double eps, const double s,
                           const uint rows, const uint cols)
{
	const uint2 gid = {get_global_id(1), get_global_id(0)};
	if (gid.y >= rows || gid.x >= cols)
		return;
	const uint i = gid.y * cols + gid.x;

	double tmp[4] = {grad[3 * i], grad[3 * i + 1],
	                 gid.x + 1 < cols ? grad[3 * (i + 1) + 1] : 0.,
	                 gid.y + 1 < rows ? grad[3 * (i + cols)] : 0.};

	double e;
	do {
		e = tmp[0] - tmp[1] + tmp[2] - tmp[3];
		const double c = .25 * s * e;
		tmp[0] -= c;
		tmp[1] += c;
		if (gid.x + 1 < cols)
			tmp[2] -= c;
		if (gid.y + 1 < rows)
			tmp[3] += c;
	} while (fabs(e) > eps);

	grad[3 * i] = tmp[0];
	grad[3 * i + 1] = tmp[1];
}*/

//______________________________________________________________________________
/*__kernel void calc_error(__global const double* const grad,
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
}*/

//______________________________________________________________________________
__kernel void reduce(__global const double* const in,
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

//______________________________________________________________________________
__kernel void nullify(__global uint* const buf, const uint n)
{
	const size_t gid = get_global_id(0);
	if (gid < n)
		buf[gid] = 0;
}

//______________________________________________________________________________
__kernel void mutate(__global const uchar* const a, __global uint* const b,
                     const uint n)
{
	const size_t gid = get_global_id(0);
	if (gid < n)
		b[gid] = a[gid];
}

//______________________________________________________________________________
__kernel void scan(__global uint* const a, __global uint* const sums,
                   __local uint* const tmp, const uint n)
{
	const size_t gid = get_global_id(0),
	             lid = get_local_id(0),
	             wgs = get_local_size(0);

	tmp[2 * lid] = 2 * gid < n ? a[2 * gid] : 0;
	tmp[2 * lid + 1] = 2 * gid + 1 < n ? a[2 * gid + 1] : 0;

	for (ushort d = wgs, ofs = 1; d; d >>= 1, ofs <<= 1) {
		barrier(CLK_LOCAL_MEM_FENCE);
		if (lid < d)
			tmp[(2 * lid + 2) *
			    ofs - 1] += tmp[(2 * lid + 1) * ofs - 1];
	}

	for (ushort ofs = wgs; ofs; ofs >>= 1) {
		barrier(CLK_LOCAL_MEM_FENCE);
		const size_t i = 2 * (lid + 1) * ofs - 1;
		if (i + ofs < 2 * wgs)
			tmp[i + ofs] += tmp[i];
	}

	barrier(CLK_LOCAL_MEM_FENCE);
	if (2 * gid < n)
		a[2 * gid] = tmp[2 * lid];
	if (2 * gid + 1 < n)
		a[2 * gid + 1] = tmp[2 * lid + 1];
	if (!lid)
		sums[get_group_id(0)] = tmp[2 * wgs - 1];
}

//______________________________________________________________________________
__kernel void prescan(__global uint* const a, __global uint* const sums,
                      __local uint* const tmp, const uint n)
{
	const size_t gid = get_global_id(0),
	             lid = get_local_id(0),
	             wgs = get_local_size(0);

	tmp[2 * lid] = 2 * gid < n ? a[2 * gid] : 0;
	tmp[2 * lid + 1] = 2 * gid + 1 < n ? a[2 * gid + 1] : 0;

	for (ushort d = wgs, ofs = 1; d; d >>= 1, ofs <<= 1) {
		barrier(CLK_LOCAL_MEM_FENCE);
		if (lid < d)
			tmp[(2 * lid + 2) * ofs  - 1] += tmp[(2 * lid + 1) *
			                                     ofs  - 1];
	}

	barrier(CLK_LOCAL_MEM_FENCE);
	if (!lid) {
		sums[get_group_id(0)] = tmp[2 * wgs - 1];
		tmp[2 * wgs - 1] = 0;
	}

	for (ushort d = 1, ofs = wgs; d <= wgs; d <<= 1, ofs >>= 1) {
		barrier(CLK_LOCAL_MEM_FENCE);
		if (lid < d) {
			const size_t i = (2 * lid + 1) * ofs  - 1,
			             j = (2 * lid + 2) * ofs  - 1;
			const uint t = tmp[i];
			tmp[i] = tmp[j];
			tmp[j] += t;
		}
	}

	barrier(CLK_LOCAL_MEM_FENCE);
	if (2 * gid < n)
		a[2 * gid] = tmp[2 * lid];
	if (2 * gid + 1 < n)
		a[2 * gid + 1] = tmp[2 * lid + 1];
}

//______________________________________________________________________________
__kernel void add(__global uint* const a, __global const uint* const inc,
                  const uint n)
{
	const size_t gid = get_global_id(0);
	if (gid < n)
		a[gid] += inc[gid / (2 * get_local_size(0))];
}
