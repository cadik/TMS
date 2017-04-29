namespace cl
{

//______________________________________________________________________________
inline mem::mem() : object{}
{
}

//______________________________________________________________________________
inline mem::mem(const cl_mem& id) : object{id}
{
}

//______________________________________________________________________________
inline void mem::retain() const
{
	clRetainMemObject(id);
}

//______________________________________________________________________________
inline void mem::release() const
{
	clReleaseMemObject(id);
}

}// namespace cl
