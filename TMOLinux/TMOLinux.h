class TMO;

int CloseLibrary(wchar_t* filename, TMO** pOperators);
int OpenLibrary(wchar_t* filename, TMO** pOperators);
int EnumLibraries(wchar_t** filenames, int maxcount);
int OperatorCount(wchar_t* filename);
