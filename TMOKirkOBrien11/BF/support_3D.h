/*! \file
  \verbatim

     Copyright (c) 2006, Sylvain Paris and Fr�do Durand

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



#ifndef __SUPPORT_FFT_3D__
#define __SUPPORT_FFT_3D__

#include <cmath>

#include <errno.h>
#include <stdio.h>

#include <sstream>
#include <string>

#include <fftw3.h>

#include "msg_stream.h"


namespace FFT{

  /*
     ####################
     # class Support_3D #
     ####################
  */


  
  //! Type for support of FFT calculation in 3D.
  class Support_3D {
    
  public:

    typedef unsigned int  size_type;

    typedef double        real_type;
    typedef fftw_complex  complex_type;
    
    typedef real_type*    space_data_type;
    typedef complex_type* frequency_data_type;
    typedef void*         buffer_type;
    
    
    //! Needs the \e width, \e height and \e depth in the space domain of the real
    //! function with which it will deal.
    inline Support_3D(const size_type w,
		      const size_type h,
		      const size_type d);

    inline ~Support_3D();


    //@{
    //! Give the size of the space domain.
    inline size_type x_space_size() const;
    inline size_type y_space_size() const;
    inline size_type z_space_size() const;
    //@}

    //@{
    //! Give the size of the frequency domain.
    inline size_type x_frequency_size() const;
    inline size_type y_frequency_size() const;
    inline size_type z_frequency_size() const;
    //@}
   
    //! Create a buffer to store data.
    inline buffer_type create_buffer() const;

    //! Destroy the buffer.    
    inline void destroy_buffer(buffer_type buffer) const;



    
    //! Put space data into a 3D array.
    template <typename Real_data_3D_array>
    void save_space_data(Real_data_3D_array* const array) const;

    //! Get space data from a 3D array.
    template <typename Real_data_3D_array>
    void load_space_data(const Real_data_3D_array& array);


    
    
    //! Put frequency data into a 3D array.
    template <typename Complex_data_3D_array>
    void save_frequency_data(Complex_data_3D_array* const array) const;

    //! Get frequency data from a 3D array.
    template <typename Complex_data_3D_array>
    void load_frequency_data(const Complex_data_3D_array& array);


    
    
    //! Put the space data into a buffer.
    inline void save_space_data_into_buffer(buffer_type buffer) const;

    //! Get the space data from a buffer.
    inline void load_space_data_from_buffer(buffer_type buffer);
    



    //! Put the frequency data into a buffer.
    inline void save_frequency_data_into_buffer(buffer_type buffer) const;

    //! Get the frequency data from a buffer.
    inline void load_frequency_data_from_buffer(buffer_type buffer);


    

    //! Compute the FFT.
    inline void space_to_frequency();

    //! Compute the inverse FFT.
    inline void frequency_to_space();



    //! Apply the normalization coefficient (image_size)^(-0.5) to all the values.
    template<typename Iterator>
    inline void half_normalize(Iterator begin,Iterator end) const;

    //! Apply the normalization coefficient (image_size)^(-1) to all the values.
    template<typename Iterator>
    inline void full_normalize(Iterator begin,Iterator end) const;


    //! Compute a fast, one-by-one, in-place, complex-valued multiplication
    //! of the frequency data.
    inline void multiply_frequency_data_by_buffer(buffer_type buffer,
						  const bool  scale = false);



    //! Convert a buffer into space data
    template<typename Real_data_3D_array>
    void buffer_to_space_data(const buffer_type& buffer,
			      Real_data_3D_array* const array) const;


     //! Convert a buffer into frequency data
    template<typename Real_data_3D_array>
    void buffer_to_frequency_data(const buffer_type& buffer,
				  Real_data_3D_array* const array) const;
   

    // ####################
    // # Static functions #
    // ####################


    //! Set the flags for FFTW.
    static inline void set_fftw_flags(const unsigned flags);
    
    //! Set the file name of the wisdom file. "" not to use a file.
    static inline void set_wisdom_file(const std::string file_name);

    //! Set 'autosave' mode.
    static inline void set_auto_save(const bool auto_save);
    
    //! Save the wisdom file; nothing fancy.
    static inline void save_wisdom();
    
  private:

    //@{
    //! FFTW structure.
    fftw_plan forward,backward;
    //@}

    //! Real data (in space domain).
    space_data_type space_data;

    //! Complex data (in frequency domain).
    frequency_data_type frequency_data;

    //@{
    //! Size of the real data.
    
    const size_type width;
    const size_type height;
    const size_type depth;
    //@}

    //@{
    //! Size of the complex data.
    
    const size_type support_w;
    const size_type support_h;
    const size_type support_d;
    const size_type mem_size;
    //@}
    
    static inline size_type index(const size_type x,
				  const size_type y,
				  const size_type z,
				  const size_type w,
				  const size_type h,
				  const size_type d);
   

    // ####################
    // # Static functions #
    // ####################


    //! See set_fftw_flags()
    static unsigned fftw_flags;

    //! Flag to indicate whether a wisdom file has been set or not
    static bool wisdom_file_set;
    
    //! See set_wisdom_file()
    static std::string wisdom_file;

    //! Number of created support.
    static size_type support_number;

    //! If true, the wisdom file is saved each time the last support is deleted.
    //! This induces useless disk access if several supports are created then deleted.
    static bool auto_save_wisdom;

    //! Indicates if the wisdom has been loaded. Avoid reading it several times.
    static bool wisdom_loaded;
    

    //! Must be called in a support constructor.
    static inline void new_support();

    //! Must be called in a support destructor.
    static inline void delete_support();

  };



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


  


  
  

  void Support_3D::set_fftw_flags(const unsigned flags){
    fftw_flags = flags;
  }

  void Support_3D::set_wisdom_file(const std::string file_name){
    wisdom_file     = file_name;
    wisdom_file_set = true;
  }

  void Support_3D::set_auto_save(const bool auto_save){
    auto_save_wisdom = auto_save;
  }


  void Support_3D::save_wisdom(){
    FILE* f = fopen(wisdom_file.c_str(),"w");

    if (f != NULL){ // System says "ok".
      fftw_export_wisdom_to_file(f);
      fflush(f);
      fclose(f);
    }
    else{ // System says "error".
      Message::warning<<"error while opening (write mode) FFTW wisdom file \""
		      <<wisdom_file<<"\"\n"
		      <<"no file saved"<<Message::done;
    }
  }
  

  void Support_3D::new_support(){
    support_number++;

    // Use a default value for the wisdom file if not yet set.
    if(!wisdom_file_set){
      char* proxy = std::getenv("FFTW_WISDOM");
      wisdom_file = (proxy == NULL) ? std::string() : std::string(proxy);
      wisdom_file_set = true;
    }
     
 
    // Job is done if not the first support or if not file support required.
    if ((wisdom_loaded)||(support_number > 1)||(wisdom_file.size() == 0)){
      return;
    }

    
    // If it is the first support, try to load the wisdom data.

    wisdom_loaded = true;

    
    // Opening the wisdom file.
    FILE* f = fopen(wisdom_file.c_str(),"r");
    
    if (f != NULL){ // System says "ok"
      if(fftw_import_wisdom_from_file(f) == 0){ // But FFTW says "error".
	Message::warning<<"error while reading FFTW wisdom file \""
			<<wisdom_file<<"\""<<Message::done;
      }
      fclose(f);
    }
    else if (errno != ENOENT){ // System says "error".
      Message::warning<<"error while opening (read mode) FFTW wisdom file \""
		      <<wisdom_file<<"\""<<Message::done;
    }
    // Wisdom is now initialized.

  }


  
  void Support_3D::delete_support(){
    support_number--;

    // Job is done if not the last support.
    if ((support_number > 0)||(wisdom_file.size() == 0)||(!auto_save_wisdom)){
      return;
    }

    // If it is the last support, try to save the wisdom data.

    save_wisdom();
  }


  
  // ######################
  // # Non-static section #
  // ######################



  /*!
    If \e wisdom_file is specified, uses the optimized process. This may
    take time at the first run but it goes faster during the next ones.
    Informations are read from and written into \e wisdom_file.

    TODO: if the file does not exist, it segfaults. Use the \e touch shell
    command to workaround.
  */
  Support_3D::Support_3D(const size_type w,
			 const size_type h,
			 const size_type d):
    width(w),
    height(h),
    depth(d),
    
    support_w(w),
    support_h(h),
    support_d(d/2+1),
    mem_size(w*h*(d/2+1)*sizeof(complex_type)) {

    new_support();
    
    // Pointer to the space data.
    space_data = reinterpret_cast<space_data_type>(fftw_malloc(mem_size));

    // Same pointer to perform "in place" computation.
    frequency_data = reinterpret_cast<frequency_data_type>(space_data);

    
    // We make the best measurement possible because we use the wisdom system.
    forward  = fftw_plan_dft_r2c_3d(width,height,depth,
				    space_data,frequency_data,
				    fftw_flags);
    
    backward = fftw_plan_dft_c2r_3d(width,height,depth,
				    frequency_data,space_data,
				    fftw_flags);
  }





  
  Support_3D::~Support_3D(){

    delete_support();
    
    // Free memory.
    fftw_destroy_plan(forward);
    fftw_destroy_plan(backward);

    fftw_free(space_data);
  }



  
  Support_3D::size_type Support_3D::x_space_size() const{
    return width;
  }


  
  Support_3D::size_type Support_3D::y_space_size() const{
    return height;
  }


  Support_3D::size_type Support_3D::z_space_size() const{
    return depth;
  }


  Support_3D::size_type Support_3D::x_frequency_size() const{
    return support_w;
  }


  Support_3D::size_type Support_3D::y_frequency_size() const{
    return support_h;
  }


  Support_3D::size_type Support_3D::z_frequency_size() const{
    return support_d;
  }


  
  
  Support_3D::size_type Support_3D::index(const size_type x,
					  const size_type y,
					  const size_type z,
					  const size_type w,
					  const size_type h,
					  const size_type d){
    return z + y*d + x*h*d;
  }


  
  
  Support_3D::buffer_type Support_3D::create_buffer() const{

    return /*fftw_*/malloc(mem_size);
  }

  void Support_3D::destroy_buffer(buffer_type buffer) const{
    
    fftw_free(buffer);
  }

  

  
  /*!
    \e Real_data_3D_array must provide the following functions: resize(),
    and a (x,y,z) access operator.
  */
  template <typename Real_data_3D_array>
  void Support_3D::save_space_data(Real_data_3D_array* const array) const{

    array->resize(width,height,depth);

    for(size_type x=0;x<width;x++){
      for(size_type y=0;y<height;y++){
	for(size_type z=0;z<depth;z++){
	  
	  const size_type i = index(x,y,z,support_w,support_h,2*support_d);

	  (*array)(x,y,z) = space_data[i];
	}
      }
    }
  }

  /*!
    Real_data_3D_array must provide the following functions:
    x_size(), y_size(), z_size() and a (x,y,z) access operator.
  */
  template <typename Real_data_3D_array>
  void Support_3D::load_space_data(const Real_data_3D_array& array){

    if ((array.x_size()!=width)||(array.y_size()!=height)||(array.z_size()!=depth)){
      Message::error<<"Support_3D::load_space_data : size does not match "
		    <<array.x_size()<<'x'<<array.y_size()<<'x'<<array.z_size()
		    <<" instead of "<<width<<'x'<<height<<'x'<<depth<<'\n'
		    <<Message::done;
    }

    for(size_type x=0;x<width;x++){
      for(size_type y=0;y<height;y++){
	for(size_type z=0;z<depth;z++){
	  const size_type i = index(x,y,z,support_w,support_h,2*support_d);

	  space_data[i] = array(x,y,z);
	}
      }
    }

  }



  
  /*!
    \e Complex_data_3D_array must provide the following functions: resize(),
    and a (x,y,z) access operator. Complex_data_3D_array::value_type is assumed
    to be compatible with the standard complex<double> structure.
  */
  template <typename Complex_data_3D_array>
  void Support_3D::save_frequency_data(Complex_data_3D_array* const array) const{

    typedef typename Complex_data_3D_array::value_type cplx_type;
    
    array->resize(support_w,support_h,support_d);

    for(size_type x=0;x<support_w;x++){
      for(size_type y=0;y<support_h;y++){
	for(size_type z=0;z<support_d;z++){
	  const size_type i = index(x,y,z,support_w,support_h,support_d);
	  
	  (*array)(x,y,z) = cplx_type(frequency_data[i][0],frequency_data[i][1]);
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
  template <typename Complex_data_3D_array>
  void Support_3D::load_frequency_data(const Complex_data_3D_array& array){

    if ((array.x_size()!=support_w)||(array.y_size()!=support_h)||(array.z_size()!=support_d)){
      Message::error<<"Support_3D::load_frequency_data : size does not match "
		    <<array.x_size()<<'x'<<array.y_size()<<'x'<<array.z_size()<<" instead of "<<support_w<<'x'<<support_h<<'x'<<support_d<<'\n'
		    <<Message::done;
    }

    for(size_type x=0;x<support_w;x++){
      for(size_type y=0;y<support_h;y++){
	for(size_type z=0;z<support_d;z++){
	  
	  const size_type i = index(x,y,z,support_w,support_h,support_d);

	  frequency_data[i][0] = array(x,y,z).real();
	  frequency_data[i][1] = array(x,y,z).imag();
	}
      }
    }

  }



  
  
  void Support_3D::save_space_data_into_buffer(buffer_type buffer) const{
    memcpy(buffer,space_data,mem_size);
  }

  void Support_3D::load_space_data_from_buffer(buffer_type buffer){
    memcpy(space_data,buffer,mem_size);
  }

  

  
  void Support_3D::save_frequency_data_into_buffer(buffer_type buffer) const{
    memcpy(buffer,frequency_data,mem_size);
  }

  void Support_3D::load_frequency_data_from_buffer(buffer_type buffer){
    memcpy(frequency_data,buffer,mem_size);
  }





  void Support_3D::space_to_frequency(){
    fftw_execute(forward);
  }


  void Support_3D::frequency_to_space(){
    fftw_execute(backward);
  }


  template<typename Iterator>
  void Support_3D::half_normalize(Iterator begin,Iterator end) const{
    
    const real_type scale = 1.0 / std::sqrt(static_cast<real_type>(width*height*depth));
    
    for(Iterator i=begin;i!=end;i++){
      *i *= scale;
    }
  }
  
  

  template<typename Iterator>
  void Support_3D::full_normalize(Iterator begin,Iterator end) const{
    
    const real_type scale = 1.0 / static_cast<real_type>(width*height*depth);
    
    for(Iterator i=begin;i!=end;i++){
      *i *= scale;
    }
  }


  void Support_3D::multiply_frequency_data_by_buffer(buffer_type buffer,
						     const bool  scale){
    
    const size_type size = support_w*support_h*support_d;
    const frequency_data_type freq_buffer =
      reinterpret_cast<frequency_data_type>(buffer);

    const real_type s = scale ? 1.0/(width*height*depth) : 1.0;
    
    for(size_type i=0;i<size;i++){

      const double a_re = frequency_data[i][0];
      const double a_im = frequency_data[i][1];
      
      const double b_re = freq_buffer[i][0];
      const double b_im = freq_buffer[i][1];
      
      frequency_data[i][0] = (a_re*b_re - a_im*b_im)*s;
      frequency_data[i][1] = (a_re*b_im + a_im*b_re)*s;
      
    }
    
  }


  
  template<typename Real_data_3D_array>
  void Support_3D::buffer_to_space_data(const buffer_type& buffer,
					Real_data_3D_array* const array) const{
    
    space_data_type space_buffer = reinterpret_cast<space_data_type>(buffer);

    array->resize(width,height);
    
    for(size_type x=0;x<width;x++){
      for(size_type y=0;y<height;y++){
	for(size_type z=0;z<depth;z++){
	  const size_type i = index(x,y,z,support_w,support_h,2*support_d);

	  (*array)(x,y,z) = space_buffer[i];
	}
      }
    }

  }

  


  template<typename Complex_data_3D_array>
  void Support_3D::buffer_to_frequency_data(const buffer_type& buffer,
					 Complex_data_3D_array* const array) const{
    
    typedef typename Complex_data_3D_array::value_type cplx_type;

    frequency_data_type frequency_buffer = reinterpret_cast<frequency_data_type>(buffer);

    array->resize(support_w,support_h,support_d);

    for(size_type x=0;x<support_w;x++){
      for(size_type y=0;y<support_h;y++){
	for(size_type z=0;z<support_d;z++){
	  const size_type i = index(x,y,z,support_w,support_h,support_d);

	  (*array)(x,y,z) = cplx_type(frequency_buffer[i][0],frequency_buffer[i][1]);
	}
      }
    }
      
  }

}


#endif
