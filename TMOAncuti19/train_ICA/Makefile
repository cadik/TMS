all: train_ica

train_ica: train_ICA.cpp
	g++ train_ICA.cpp -o train_ica -std=c++17 -litpp -llapack -lblas `pkg-config --cflags --libs opencv4`

clean:
	rm -f train_ica
