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

#ifndef __CONVOLUTION_FFT_3D__
#define __CONVOLUTION_FFT_3D__

#include <algorithm>

#include "fftw3.h"

#include "support_3D.h"



namespace FFT{

  
  

  //! Compute h=f*g for f,g and h real functions.
  template<typename Function>
  inline void convolve_3D(const Function& f,
			  const Function& g,
			  Function* const h);

  
  //! Compute h=f*g for f,g and h real functions. Support is given,
  //! goes faster for multiple computations.
  template<typename Function>
  inline void convolve_3D(const Function& f,
			  const Function& g,
			  Function* const h,
			  Support_3D&     support);

  
  //! This function mixes frequency representation for f and space
  //! representation for g to provide a fast process when f is convolved
  //! several times with different g.
  //! Compute h=f*g for f,g and h real functions. Support is given,
  //! goes faster for multiple computations.
  //! f is given through its Fourier transform F in a Support_3D::buffer_type.
  template<typename Function>
  void convolve_3D(const Support_3D::buffer_type& F,
		const Function& g,
		Function* const h,
		Support_3D&     support);

  
  //! This function mixes frequency representation for f and g, space
  //! representation for the result to provide a fast process when f is convolved
  //! several times with different g.
  //! Compute h=f*g for f,g and h real functions. Support is given,
  //! goes faster for multiple computations.
  //! f and are given through their Fourier transform F and G in a Support_3D::buffer_type.
  template<typename Function>
  void convolve_3D(const Support_3D::buffer_type& F,
		   const Support_3D::buffer_type& G,
		   Function* const h,
		   Support_3D&        support);
  
  
  //! Adapt the function so that the point (0,0,0) is at (x,y,z) before
  //! the funtion call.
  template<typename Function>
  void move_center_3D(const Function&    f,
		      const unsigned int x,
		      const unsigned int y,
		      const unsigned int z,
		      Function* const    result);



  


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

  
  template<typename Function>
  void convolve_3D(const Support_3D::buffer_type& F,
		const Function& g,
		Function* const h,
		Support_3D&     support){

    // ################
    // ### FFT of g ###
    // ################

    // Store g data in the fftw structures.    
    support.load_space_data(g);
    
    // Compute the FFT.
    support.space_to_frequency();



    
    // Multiplication "in place" in the frequency space.    
    support.multiply_frequency_data_by_buffer(F,true);


    

    // ###################
    // ### Inverse FFT ###
    // ###################
    
    // Compute inverse FFT transform.
    support.frequency_to_space();

    // Return data into the original structures    
    support.save_space_data(h);
  }


  /*!
    A Function must have members \e x_size(), \e y_size() and \e y_size() to get
    the dimensions of the data. It must have an access operator \e (x,y,z) too.
  */
  template<typename Function>
  void convolve_3D(const Function& f,
		   const Function& g,
		   Function* const h,
		   Support_3D&     support){

    // ################
    // ### FFT of f ###
    // ################

    // We need a temporary buffer.
    Support_3D::buffer_type buffer = support.create_buffer();

    // Store f data in the fftw structures.
    support.load_space_data(f);
    
    // Compute the FFT.
    support.space_to_frequency();
    
    // Copy the result in a safe place.
    support.save_frequency_data_into_buffer(buffer);


    convolve_3D(buffer,g,h,support);
    
    support.destroy_buffer(buffer);
    
  }

  

  template<typename Function>
  void convolve_3D(const Function& f,
		const Function& g,
		Function* const h){

    Support_3D s(f.x_size(),f.y_size(),f.z_size());

    convolve_3D(f,g,h,s);
  }

  
  template<typename Function>
  void convolve_3D(const Support_3D::buffer_type& F,
		   const Support_3D::buffer_type& G,
		   Function* const h,
		   Support_3D&     support){

    // Store G data in the fftw structures.    
    support.load_frequency_data_from_buffer(G);
    
    // Multiplication "in place" in the frequency space.    
    support.multiply_frequency_data_by_buffer(F,true);
    

    // ###################
    // ### Inverse FFT ###
    // ###################
    
    // Compute inverse FFT transform.
    support.frequency_to_space();

    // Return data into the original structures    
    support.save_space_data(h);


  }


  
  template<typename Function>
  void move_center_3D(const Function&    f,
		      const unsigned int x_center,
		      const unsigned int y_center,
		      const unsigned int z_center,
		      Function* const    result){

    typedef unsigned int size_type;

    const size_type width  = f.x_size();
    const size_type height = f.y_size();
    const size_type depth  = f.z_size();

    Function g(width,height);
    
    for(size_type x=0;x<width;x++){
      for(size_type y=0;y<height;y++){
	for(size_type z=0;z<depth;z++){
	(*result)(x,y,z) = f((x + width - x_center) % width,
			     (y + height - y_center) % height,
			     (z + depth - z_center) % depth);
	}
      }
    }
  }
  
}

#endif
