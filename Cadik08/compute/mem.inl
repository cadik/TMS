namespace cl
{

inline mem::mem() : object{}
{
}
	
inline mem::mem(const cl_mem& id) : object{id}
{
}

inline void mem::retain() const
{
	clRetainMemObject(id);
}

inline void mem::release() const
{
	clReleaseMemObject(id);
}

}// namespace cl
