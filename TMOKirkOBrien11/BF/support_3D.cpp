/*! \file
  \verbatim

     Copyright (c) 2006, Sylvain Paris and Frédo Durand

     Permission is hereby granted, free of charge, to any person
     obtaining a copy of this software and associated documentation
     files (the "Software"), to deal in the Software without
     restriction, including without limitation the rights to use, copy,
     modify, merge, publish, distribute, sublicense, and/or sell copies
     of the Software, and to permit persons to whom the Software is
     furnished to do so, subject to the following conditions:

     The above copyright notice and this permission notice shall be
     included in all copies or substantial portions of the Software.

     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
     MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
     NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
     WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
     OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
     DEALINGS IN THE SOFTWARE.

   \endverbatim  
*/


#include <cstring>
#include "support_3D.h"

namespace FFT{
  
  // Static values initialization
  
  unsigned              Support_3D::fftw_flags       = FFTW_EXHAUSTIVE;
  bool                  Support_3D::wisdom_file_set  = false;
  std::string           Support_3D::wisdom_file      = std::string();
  Support_3D::size_type Support_3D::support_number   = 0;
  bool                  Support_3D::auto_save_wisdom = true;
  bool                  Support_3D::wisdom_loaded    = false;
}
