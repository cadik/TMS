#ifndef TMOCADIK08_CORE_CL
#define TMOCADIK08_CORE_CL

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

#endif
