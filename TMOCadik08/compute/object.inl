namespace cl
{

//______________________________________________________________________________
template <class D, typename T>
inline object<D, T>::object() : id{nullptr}
{
}

//______________________________________________________________________________
template <class D, typename T>
inline object<D, T>::object(const object& img) : id{img.id}
{
	retain();
}

//______________________________________________________________________________
template <class D, typename T>
inline object<D, T>::object(const T& id) : id{id}
{
}

//______________________________________________________________________________
template <class D, typename T>
inline object<D, T>::~object()
{
	release();
}

//______________________________________________________________________________
template <class D, typename T>
inline object<D, T>::operator const T&() const
{
	return id;
}

//______________________________________________________________________________
template <class D, typename T>
inline object<D, T>::operator T&()
{
	return id;
}

//______________________________________________________________________________
template <class D, typename T>
inline object<D, T>::operator const T*() const
{
	return &id;
}

//______________________________________________________________________________
template <class D, typename T>
inline object<D, T>::operator T*()
{
	return &id;
}

//______________________________________________________________________________
template <class D, typename T>
inline object<D, T>& object<D, T>::operator=(const object& rhs)
{
	release();
	id = rhs.id;
	retain();
	return *this;
}

//______________________________________________________________________________
template <class D, typename T>
inline void object<D, T>::retain() const
{
	if (id)
		static_cast<const D*>(this)->retain();
}

//______________________________________________________________________________
template <class D, typename T>
inline void object<D, T>::release() const
{
	if (id)
		static_cast<const D*>(this)->release();
}

}// namespace cl
