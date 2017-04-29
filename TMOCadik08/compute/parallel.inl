#include <string>

namespace cl
{

//______________________________________________________________________________
inline parallel::parallel(const cl_device_type type) :
                          host{type},
                          gpu{host, type},
                          env{{CL_CONTEXT_PLATFORM,
                               reinterpret_cast<cl_context_properties>
                               ((const cl_platform_id&) host), 0}, gpu},
                          lms{gpu.info<cl_ulong>(CL_DEVICE_LOCAL_MEM_SIZE)},
                          wgs{gpu.info<size_t>(CL_DEVICE_MAX_WORK_GROUP_SIZE)},
                          maba{gpu.info<cl_uint>(CL_DEVICE_MEM_BASE_ADDR_ALIGN) /
                               8}
{
}

//______________________________________________________________________________
inline command_queue parallel::create_command_queue(const cl_command_queue_properties opts) const
{
	return {gpu, env, opts};
}

//______________________________________________________________________________
inline buffer parallel::create_buffer(const cl_mem_flags opts, const size_t size,
                                      const void* const ptr) const
{
	return {env, opts, size, ptr};
}

//______________________________________________________________________________
inline program parallel::create_program(const char* source,
                                        const char* const opts) const
{
	return {env, source, gpu, opts};
}

//______________________________________________________________________________
inline size_t parallel::get_wgs() const
{
	return wgs;
}

//______________________________________________________________________________
inline unsigned parallel::get_maba() const
{
	return maba;
}

//______________________________________________________________________________
inline std::string parallel::info() const
{
	const char* const dev{gpu.info<char*>(CL_DEVICE_NAME)};
	std::string info{dev};
	delete[] dev;
	info += "\n|_ local mem. size: " + std::to_string(lms) + " B\n" +
	        "|_ max. work group size: " + std::to_string(wgs) + '\n' +
	        "|_ device address bits: " + std::to_string(gpu.info<size_t>(CL_DEVICE_ADDRESS_BITS)) + " b\n" +
	        "|_ device mem. base addr. align: " + std::to_string(maba) + " B\n"; 
	return info;
}

}// namespace cl
