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



#ifndef __LINEAR_BF__
#define __LINEAR_BF__

#include <cmath>

#include "array.h"
#include "fft_3D.h"
#include "math_tools.h"

#ifdef CHRONO
#include <iostream>
#include "chrono.h"
#endif

namespace Image_filter{

  typedef double       real_type;
  typedef unsigned int size_type;

  template <typename Array>
  void linear_BF(const Array&    input,
		 const Array&    base,
		 const real_type space_sigma,
		 const real_type range_sigma,
		 const real_type space_sampling,
		 const real_type range_sampling,
		 const bool      early_division,
		 Array* const    result);

  template <typename Array>
  inline void linear_BF(const Array&    input,
			const real_type space_sigma,
			const real_type range_sigma,
			Array* const    result);
  


  
/*
  
  #############################################
  #############################################
  #############################################
  ######                                 ######
  ######   I M P L E M E N T A T I O N   ######
  ######                                 ######
  #############################################
  #############################################
  #############################################
  
*/


  template <typename Array>
  void linear_BF(Array&    input,
		 real_type space_sigma,
		 real_type range_sigma,
		 Array* result){

    linear_BF(input,
	      space_sigma,range_sigma,
	      space_sigma,range_sigma,
	      result);
  }



  template <typename Array>
  void linear_BF(const Array&    input,
		 const real_type space_sigma,
		 const real_type range_sigma,
		 const real_type space_sampling,
		 const real_type range_sampling,
		 Array* const    result){

    using namespace std;

    typedef Array               real_array_2D_type;
    typedef Array_3D<real_type> real_array_3D_type;
    
    const size_type width  = input.x_size();
    const size_type height = input.y_size();

    const real_type sigma_r = range_sigma / range_sampling; 
    const real_type sigma_s = space_sigma / space_sampling;

    const size_type padding_xy = static_cast<size_type>(2.0 * sigma_s) + 1;
    const size_type padding_z  = static_cast<size_type>(2.0 * sigma_r) + 1;

    const real_type input_min   = *min_element(input.begin(),input.end());
    const real_type input_max   = *max_element(input.begin(),input.end());
    const real_type input_delta = input_max - input_min;
    
    const size_type small_width  = static_cast<size_type>((width  - 1) / space_sampling) + 1 + 2*padding_xy; 
    const size_type small_height = static_cast<size_type>((height - 1) / space_sampling) + 1 + 2*padding_xy; 
    const size_type small_depth  = static_cast<size_type>(input_delta / range_sampling) + 1 + 2*padding_z; 

    FFT::Support_3D fft_support(small_width,small_height,small_depth);

#ifdef CHRONO
    Chrono chrono("filter");
    chrono.start();

    Chrono chrono_down("downsampling");
    chrono_down.start();
#endif
    
    real_array_3D_type w(small_width,small_height,small_depth,0.0);
    real_array_3D_type iw(small_width,small_height,small_depth,0.0);

    for(size_type x=0,x_end=width;x<x_end;x++){
      for(size_type y=0,y_end=height;y<y_end;y++){
      
	const size_type small_x = static_cast<size_type>(1.0*x / space_sampling + 0.5) + padding_xy;
	const size_type small_y = static_cast<size_type>(1.0*y / space_sampling + 0.5) + padding_xy;
	const size_type small_z = static_cast<size_type>((input(x,y)-input_min) / range_sampling + 0.5) + padding_z;

	w(small_x,small_y,small_z)  += 1.0;
	iw(small_x,small_y,small_z) += input(x,y);
      } // END OF for y
    } // END OF for x
  
    real_array_3D_type kernel(small_width,small_height,small_depth);

    const size_type half_width  = small_width/2;
    const size_type half_height = small_height/2;
    const size_type half_depth  = small_depth/2;
  
    for(size_type x=0,x_end=small_width;x<x_end;x++){

      const real_type X = static_cast<real_type>(x) - ((x>half_width) ? small_width : 0.0);
    
      for(size_type y=0,y_end=small_height;y<y_end;y++){

	const real_type Y = static_cast<real_type>(y) - ((y>half_height) ? small_height : 0.0);
      
	for(size_type z=0,z_end=small_depth;z<z_end;z++){

	  const real_type Z = static_cast<real_type>(z) - ((z>half_depth) ? small_depth : 0.0);

	  const real_type rr = (X*X + Y*Y) / (sigma_s*sigma_s) + Z*Z / (sigma_r*sigma_r);	
	
	  kernel(x,y,z) = exp(-rr*0.5);
	
	} // END OF for z
      } // END OF for y
    } // END OF for x

#ifdef CHRONO
    chrono_down.stop();

    Chrono chrono_convolution("convolution");
    chrono_convolution.start();
#endif  
  
    FFT::convolve_3D(iw,kernel,&iw,fft_support);
    FFT::convolve_3D(w,kernel,&w,fft_support);

#ifdef CHRONO
    chrono_convolution.stop();

    Chrono chrono_nonlinearities("nonlinearities");
    chrono_nonlinearities.start();
#endif

    result->resize(width,height);
    
    for(size_type x=0,x_end=width;x<x_end;x++){
      for(size_type y=0,y_end=height;y<y_end;y++){
      
	const real_type z = input(x,y) - input_min;

	const real_type IW = Math_tools::trilinear_interpolation(iw,
								 static_cast<real_type>(x)/space_sampling + padding_xy,
								 static_cast<real_type>(y)/space_sampling + padding_xy,
								 z/range_sampling + padding_z);
      
	const real_type W = Math_tools::trilinear_interpolation(w,
								static_cast<real_type>(x)/space_sampling + padding_xy,
								static_cast<real_type>(y)/space_sampling + padding_xy,
								z/range_sampling + padding_z);

	(*result)(x,y) = IW / W;
      
      } // END OF for y
    } // END OF for x
    
#ifdef CHRONO
    chrono_nonlinearities.stop();
    chrono.stop();
    
    cout<<chrono.report()<<endl;
    cout<<chrono_down.report()<<endl;
    cout<<chrono_convolution.report()<<endl;
    cout<<chrono_nonlinearities.report()<<endl;
#endif
    
  }

} // END OF namespace Image_filter



#endif
