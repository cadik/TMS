namespace com
{

inline text_loader::text_loader(const std::string path) : loader{path}
{
}

inline std::string text_loader::operator()()
{
	input.seekg(0);
	const std::string content{{std::istreambuf_iterator<char>(input)}, std::istreambuf_iterator<char>()};
	return content;
}

}// namespace com
