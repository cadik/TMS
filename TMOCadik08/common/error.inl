namespace com
{

inline error::error(const std::string file, const unsigned line, const char* const func, const std::string msg) : std::runtime_error{file + ": " + std::to_string(line) + ": " + func + ": " + msg}
{
}

}// namespace com
