Name:           armadillo
Version:        12.6.x
Release:        1%{?dist}
Summary:        C++ library for linear algebra and scientific computing with MATLAB-like syntax

License:        ASL 2.0
URL:            http://arma.sourceforge.net/
Source:         http://sourceforge.net/projects/arma/files/%{name}-%{version}.tar.xz

BuildRequires:  gcc-c++, cmake, lapack-devel, arpack-devel
%{!?openblas_arches:%global openblas_arches x86_64 %{ix86} armv7hl %{power64} aarch64}
%ifarch %{openblas_arches}
BuildRequires:  openblas-devel
%endif
BuildRequires:  SuperLU-devel


%description
Armadillo is a high quality C++ library for linear algebra and scientific computing,
aiming towards a good balance between speed and ease of use.
Useful for algorithm development directly in C++,
and/or quick conversion of research code into production environments.
The syntax (API) is deliberately similar to Matlab.
The library provides efficient classes for vectors, matrices and cubes,
as well as 200+ associated functions (eg. contiguous and non-contiguous
submatrix views).  Various matrix decompositions are provided through
integration with LAPACK, or one of its high performance drop-in replacements
(eg. OpenBLAS, Intel MKL, Apple Accelerate framework, etc).
A sophisticated expression evaluator (via C++ template meta-programming)
automatically combines several operations (at compile time) to increase
speed and efficiency.
The library can be used for machine learning, pattern recognition,
computer vision, signal processing, bioinformatics, statistics, finance, etc.


%package devel
Summary:        Development headers and documentation for the Armadillo C++ library
Requires:       %{name} = %{version}-%{release}
Requires:       lapack-devel,  arpack-devel, libstdc++-devel
%ifarch %{openblas_arches}
Requires:       openblas-devel
%endif
Requires:       SuperLU-devel


%description devel
This package contains files necessary for development using the
Armadillo C++ library. It contains header files, example programs,
and user documentation (API reference guide).


%prep
%setup -q

# convert DOS end-of-line to UNIX end-of-line

for file in README.md; do
  sed 's/\r//' $file >$file.new && \
  touch -r $file $file.new && \
  mv $file.new $file
done

%build
%{cmake}
%{__make} VERBOSE=1 %{?_smp_mflags}


%install
rm -rf $RPM_BUILD_ROOT
%{__make} install DESTDIR=$RPM_BUILD_ROOT
rm -f examples/Makefile.cmake
rm -f examples/example1_win64.sln
rm -f examples/example1_win64.vcxproj
rm -f examples/example1_win64.README.txt
rm -rf examples/lib_win64


%post -p /sbin/ldconfig

%postun -p /sbin/ldconfig


%files
%{_libdir}/*.so.*
%license LICENSE.txt NOTICE.txt

%files devel
%{_libdir}/*.so
%{_libdir}/pkgconfig/%{name}.pc
%{_includedir}/armadillo
%{_includedir}/armadillo_bits/
%{_datadir}/Armadillo/
%doc README.md index.html docs.html
%doc examples armadillo_icon.png
%doc armadillo_joss_2016.pdf
%doc armadillo_lncs_2018.pdf
%doc armadillo_nicta_2010.pdf
%doc armadillo_solver_2020.pdf
%doc armadillo_spcs_2017.pdf
%doc armadillo_rcpp_2014.pdf
%doc mex_interface
