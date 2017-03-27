Tone Mapping Studio (TMS), 
2004 - 2017, Martin Cadik (cadikm@centrum.cz, http://cadik.posvete.cz/)
developers: Ondrej Hajdok, Martin Cadik, Michal Augustyn, Ondrej Fialka, Antonin Lejsek, Petr Bilek, Ondrej Pecina, Pavel Fryz, Martin Molek, Vladimir Vlkovic, Jan Brida

ABOUT 
Tone Mapping Studio is a plugin-based framework for image processing operations, 
in particular HDR tone mapping and color-to-grayscale conversions. 

If you use TMS for your research work, please cite at least one of the following papers:

@article{cadik08evaluation,
  title =      "Evaluation of HDR Tone Mapping Methods Using Essential
               Perceptual Attributes",
  author =     "Martin {\v{C}}ad\'{i}k and Michael Wimmer and Laszlo
               Neumann and Alessandro Artusi",
  year =       "2008",
  volume =     "32",
  issue =      "3",
  pages =      "330-349",
  issn =       "0097-8493",
  journal =    "Computers \& Graphics",
  keywords =   "tone mapping operators, tone mapping evaluation, image
               attributes, high dynamic range",
  URL =        "http://cadik.posvete.cz/tmo/cadik08cag.pdf",
}

@article{ cadik08perceptual,
  author    = {Martin {\v{C}}ad\'{i}k},
  title     = {Perceptual Evaluation of Color-to-Grayscale Image Conversions},
  journal   = {Comput. Graph. Forum},
  volume    = {27},
  number    = {7},
  year      = {2008},
  pages     = {1745-1754},
}

@PhdThesis{cadik08phd,
 author = "Martin {\v{C}}ad\'{i}k",
 title  = "Perceptually Based Image Quality Assessment and Image Transformations",
 school = "Department of Computer Science and Engineering,
           Faculty of Electrical Engineering,
           Czech Technical University in Prague",
 type   = "Ph.D. Thesis",
 year   = "2008",
 month  = "January",
 url    = "http://cadik.posvete.cz/diss",
}
 


The TMS package consists of 
-tmolib - a library for loading, writing and processing of LDR and HDR files in various formats
-tmocmd - a command line tool, which loads an image, runs a specific plugin, and writes the final output
-tmogui - an application, which uses the same plugins as tmocmd, to enable visual inspection of the results, interactive setting of parameters, and other operations
-TMOYourOperatorName/TMOYourOperatorNameLinux - templates for your own plugins
-TMOWard94 - an example plugin implementing a simple tone mapping method [Ward 1994]
-TMOLinux/TMOW32 - linux/Windows specific plugin code


COMPILING

To compile run in this directory 
	cmake . 
	make 
Use cmake version 2.4.8 and higher.

You need at least these libraries (and versions):
OpenExr 1.2.2-4.4       http://www.openexr.com/
tiff4 3.8.2-8           www.libtiff.org/
libjpeg 6 or higher     http://www.ijg.org/
FFTW 3.3.4              http://www.fftw.org/
Boost.MultiArray 1.6    http://www.boost.org/doc/libs/1_56_0/libs/multi_array/doc/index.html
CGAL 3.3                http://www.cgal.org/

and for GUI:
Qt3 3.3.8b              https://download.qt.io/archive/qt/3/


On Debian these are provided by packages:
libopenexr-dev
libtiff4-dev
libjpeg-dev
libqt3-mt-dev
qt3-assistant
qt3-apps-dev
libfftw3-3
libboost-all-dev


INSTALL
After compiling run make install. If you want to specify directory to install, you have to run cmake like this:
cmake -D CMAKE_INSTALL_PREFIX='/where/to/install/tmo/' .

