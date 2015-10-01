test.tml: TMOPlugin.o TMOYourOperatorName.o 
	g++ -shared -fPIC TMOYourOperatorName.o TMOPlugin.o ../TMOCmd/*.o -o ../TMOCmd/test.tml
TMOYourOperatorName.o: TMOYourOperatorName.cpp TMOYourOperatorName.h
	g++ -c -fPIC -I../tmolib/ TMOYourOperatorName.cpp
TMOPlugin.o:	
	g++ -c -fPIC -I../tmolib/ TMOPlugin.cpp
clean: 
	\rm *.o ../TMOCmd/test.tml

