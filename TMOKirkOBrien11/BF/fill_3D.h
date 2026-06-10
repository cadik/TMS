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

#ifndef __FFT_FILL_3D__
#define __FFT_FILL_3D__


namespace FFT{

  //! Fill in the real part of the array with function f. Imaginary part is
  //! not changed.
  template <typename Function,typename Complex_data_3D_array>
  void fill_real_part_3D(const Function& f,
			 Complex_data_3D_array* const array);

  //! Fill in the imaginary part of the array with function f. Real part is
  //! not changed.
  template <typename Function,typename Complex_data_3D_array>
  void fill_imaginary_part_3D(const Function& f,
			      Complex_data_3D_array* const array);



  
/*
  ###########################################################
  ###########################################################
  ###########################################################
  ##############                               ##############
  ##############  I M P L E M E N T A T I O N  ##############
  ##############                               ##############
  ###########################################################
  ###########################################################
  ###########################################################
*/



  
  /*!
    Complex_data_3D_array must provide the following functions:
    x_size(), y_size(), z_size() and a (x,y,z) access operator.
    Complex_data_3D_array::value_type is assumed to be compatible
    with the standard complex<double> structure.
  */
  template <typename Function,typename Complex_data_3D_array>
  void fill_real_part_3D(const Function& f,
			 Complex_data_3D_array* const array){

    typedef typename Complex_data_3D_array::value_type complex_type;
    typedef float                                      real_type;
    typedef unsigned int                               size_type;

    const size_type width  = array->x_size();
    const size_type height = array->y_size();
    const size_type depth  = array->z_size();
  
    const size_type half_width  = width/2;
    const size_type half_height = height/2;

    for(size_type x=0;x<width;x++){
      for(size_type y=0;y<height;y++){
	for(size_type z=0;z<depth;z++){

	  const real_type X = static_cast<real_type>(x) -
	    width*((x+half_width)/width);
      
	  const real_type Y = static_cast<real_type>(y) -
	    height*((y+half_height)/height);
	
	  const real_type Z = static_cast<real_type>(z);
 
	  (*array)(x,y,z) = complex_type(f(X,Y,Z),(*array)(x,y,z).imag());
	}
      }
    }
  }

  
    /*!
    Complex_data_3D_array must provide the following functions:
    x_size(), y_size(), z_size() and a (x,y,z) access operator.
    Complex_data_3D_array::value_type is assumed to be compatible
    with the standard complex<double> structure.
  */
  template <typename Function,typename Complex_data_3D_array>
  void fill_imaginary_part_3D(const Function& f,
			      Complex_data_3D_array* const array){

    typedef typename Complex_data_3D_array::value_type complex_type;
    typedef float                                      real_type;
    typedef unsigned int                               size_type;

    const size_type width  = array->x_size();
    const size_type height = array->y_size();
    const size_type depth  = array->z_size();
  
    const size_type half_width  = width/2;
    const size_type half_height = height/2;

    for(size_type x=0;x<width;x++){
      for(size_type y=0;y<height;y++){
	for(size_type z=0;z<depth;z++){

	  const real_type X = static_cast<real_type>(x) -
	    width*((x+half_width)/width);
      
	  const real_type Y = static_cast<real_type>(y) -
	    height*((y+half_height)/height);
	
	  const real_type Z = static_cast<real_type>(z); 
 
	  (*array)(x,y,z) = complex_type((*array)(x,y,z).real(),f(X,Y,Z));
	}
      }
    }
  }

}

#endif
