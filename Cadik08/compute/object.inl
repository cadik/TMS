namespace cl
{

template <class D, typename T>
inline object<D, T>::object() : id{nullptr}
{
}

template <class D, typename T>
inline object<D, T>::object(const object& img) : id{img.id}
{
	retain();
}

template <class D, typename T>
inline object<D, T>::object(const T& id) : id{id}
{
}

template <class D, typename T>
inline object<D, T>::~object()
{
	release();
}

template <class D, typename T>
inline object<D, T>::operator const T&() const
{
	return id;
}

template <class D, typename T>
inline object<D, T>::operator T&()
{
	return id;
}

template <class D, typename T>
inline object<D, T>::operator const T*() const
{
	return &id;
}

template <class D, typename T>
inline object<D, T>::operator T*()
{
	return &id;
}

template <class D, typename T>
inline object<D, T>& object<D, T>::operator=(const object& rhs)
{
	release(); // XXX ?
	id = rhs.id;
	retain();
	return *this;
}

template <class D, typename T>
inline void object<D, T>::retain() const
{
	if (id)
		static_cast<const D*>(this)->retain();
}

template <class D, typename T>
inline void object<D, T>::release() const
{
	if (id)
		static_cast<const D*>(this)->release();
}

}// namespace cl
