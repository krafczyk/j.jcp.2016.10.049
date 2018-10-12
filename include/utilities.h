#ifndef UTILITIES_HDR
#define UTILITIES_HDR

template<class T> T** New2DArray(const int NX, const int NY) {
	T** new_array = new T*[NX];
	for(int i=0; i<NX; ++i) {
		new_array[i] = new T[NY];
	}
	return new_array;
}

template<class T> T*** New3DArray(const int NX, const int NY, const int NZ) {
	T*** new_array = new T**[NX];
	for(int i=0; i<NX; ++i) {
		new_array[i] = new T*[NY];
		for(int j=0; j<NY; ++j) {
			new_array[i][j] = new T[NZ];
		}
	}
	return new_array;
}

template<class T> void Delete2DArray(T** array, const int NX, const int NY) {
	for(int i=0; i<NX; ++i) {
		delete[] array[i];
	}
	delete[] array;
}

template<class T> void Delete3DArray(T*** array, const int NX, const int NY, const int NZ) {
	for(int i=0; i<NX; ++i) {
		for(int j=0; j<NY; ++j) {
			delete[] array[i][j];
		}
		delete[] array[i];
	}
	delete[] array;
}

#endif
