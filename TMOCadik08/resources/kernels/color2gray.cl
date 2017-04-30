#ifdef cl_amd_fp64
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#elif cl_khr_fp64
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#endif

#include "morton.cl"
#include "core.cl"

//==============================================================================
// naive cyclical version
//______________________________________________________________________________
__kernel void cyc_calc_error(__global const double2* const grad,
                             __global double* const err,
                             const uint rows, const uint cols)
{
	const uint2 gid = {get_global_id(1), get_global_id(0)};
	if (gid.y >= rows || gid.x >= cols)
		return;

	const uint i = gid.y * cols + gid.x;
	err[i] = grad[i].x - grad[i].y +
	         ((gid.x + 1 < cols) ? grad[i + 1].y : 0.) -
	         ((gid.y + 1 < rows) ? grad[i + cols].x : 0.);
}

//______________________________________________________________________________
__kernel void cyc_correct_grad(__global double2* const grad,
                               __global const double* const err, const double s,
                               const uint rows, const uint cols, const uint n)
{
	const uint2 gid = {get_global_id(1), get_global_id(0)};
	if (gid.y >= rows || gid.x >= cols)
		return;

	const uint i = gid.y * cols + gid.x;
	const double r = .25 * s * err[i];
	if (n == 0)
		grad[i].x -= r;
	else if (n == 1)
		grad[i].y += r;
	else if (n == 2) {
		if (gid.y < rows && gid.x + 1 < cols)
			grad[i + 1].y -= r;
	}
	else if (n == 3) {
		if ((gid.y + 1 < rows && gid.x < cols))
			grad[i + cols].x += r;
	}
}

//==============================================================================
// chessboard version
enum { CHESS_EVEN, CHESS_ODD, CHESS_DIAG };

//______________________________________________________________________________
__kernel void cb_correct_grad(__global double2* const g,
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
		const double e = g[i].x - g[i].y +
		                 ((gid.x + 1 < cols) ? g[i + 1].y :
		                                       0.) -
		                 ((gid.y + 1 < rows) ? g[i + cols].x :
		                                       0.);
		const double r = .25 * s * e;

		g[i].x -= r;
		g[i].y += r;
		if (gid.x + 1 < cols)
			g[i + 1].y -= r;
		if (gid.y + 1 < rows)
			g[i + cols].x += r;

		err[i] = e;
	}
}

//______________________________________________________________________________
// local memory version
__kernel void cb_correct_grad_l(__global double2* const g,
                                __local double2* const tmp,
                                __global double* const err,
                                const double s,
                                const uint rows, const uint cols)
{
	const uint2 gid = {get_global_id(1), get_global_id(0)},
	            lid = {get_local_id(1), get_local_id(0)},
	            wgs = {get_local_size(1), get_local_size(0)};
	const uint dim = wgs.x + 1;

	const size_t i = gid.y * cols + gid.x,
	             l = lid.y * dim + lid.x;

	tmp[l] = (gid.y < rows) && (gid.x < cols) ? g[i] : (double2) (0., 0.);
	if ((lid.x == (wgs.x - 1)))
		tmp[l + 1] = ((gid.x + 1) < cols && gid.y < rows) ? g[i + 1] : (double2) (0., 0.);
	if ((lid.y == (wgs.y - 1)))
		tmp[l + dim] = ((gid.y + 1) < rows && gid.x < cols) ? g[i + cols] : (double2) (0., 0.);
	// no need to setup the tmp[(wgs.y + 1) * (wgs.x + 1) - 1] gradient, it is never used
	barrier(CLK_LOCAL_MEM_FENCE);

	if (!(lid.y % 2) && !(lid.x % 2)) {
		const double e = tmp[l].x - tmp[l].y +
		                 tmp[l + 1].y - tmp[l + dim].x;
		const double r = .25 * s * e;

		tmp[l].x -= r;
		tmp[l].y += r;
		tmp[l + 1].y -= r;
		tmp[l + dim].x += r;

		if (gid.y < rows && gid.x < cols)
			err[i] = e;
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	if (((lid.y % 2) && !(lid.x % 2)) || ((lid.x % 2) && !(lid.y % 2))) {
		const double e = tmp[l].x - tmp[l].y +
		                 tmp[l + 1].y - tmp[l + dim].x;
		const double r = .25 * s * e;

		tmp[l].x -= r;
		tmp[l].y += r;
		tmp[l + 1].y -= r;
		tmp[l + dim].x += r;

		if (gid.y < rows && gid.x < cols)
			err[i] = e;
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	if ((lid.y % 2) && (lid.x % 2)) {
		const double e = tmp[l].x - tmp[l].y +
		                 tmp[l + 1].y - tmp[l + dim].x;
		const double r = .25 * s * e;

		tmp[l].x -= r;
		tmp[l].y += r;
		tmp[l + 1].y -= r;
		tmp[l + dim].x += r;

		if (gid.y < rows && gid.x < cols)
			err[i] = e;
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	if (gid.y < rows && gid.x < cols)
		g[i] = tmp[l];
}

//==============================================================================
// hierarchical correction
//______________________________________________________________________________
// used to evaluate coarser quadtree levels
__kernel void hie_avg_grad(__global double2* const rooti,
                           const uint aligni, const uint sizei)
{
	const size_t gid = get_global_id(0);
	if (gid >= sizei)
		return;

	rooti[aligni + gid] = (rooti[aligni + sizei + (gid << 2)] +
	                       rooti[aligni + sizei + (gid << 2) | 1] +
	                       rooti[aligni + sizei + (gid << 2) | 2] +
	                       rooti[aligni + sizei + (gid << 2) | 3]);
}

//______________________________________________________________________________
__kernel void hie_tag_active(__global const uchar* const flag,
                             __global const uint* const offsets,
                             __global const morton* const zs,
                             __global morton* const parent_index,
                             const uint n_a_i)
{
	const size_t gid = get_global_id(0);
	if (gid < n_a_i && flag[gid])
		parent_index[offsets[gid]] = zs[gid];// (parent_index_old[gid / 4] << 2) | (gid % 4)
}

//______________________________________________________________________________
// for all active loops in a quadtree level i = n_a_i
__kernel void hie_calc_err(__global const double2* const root, const uint align,
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

//______________________________________________________________________________
__kernel void hie_correct_grad(__global const double2* const root0,
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

	if (pos.x >= side || pos.y >= side)
		return;

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

//==============================================================================
// double integration 
//_____________________________________________________________________________
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
