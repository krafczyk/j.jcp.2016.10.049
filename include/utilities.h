#ifndef UTILITIES_HDR
#define UTILITIES_HDR

template<class T>
class Array2D {
	public:
		Array2D() {
			this->data = NULL;
			this->_NX = 0;
			this->_NY = 0;
		}
		~Array2D() {
			if (this->data != NULL) {
				delete [] this->data;
				this->data = NULL;
			}
		}

		void init(const int NX, const int NY) {
			this->_NX = NX;
			this->_NY = NY;
			if (this->data == NULL) {
				this->data = new T[this->_NX*this->_NY];
			}
		}

		T operator()(const int i, const int j) {
			return this->data[i*this->_NY+j];
		}
		T& assign(const int i, const int j) {
			return this->data[i*this->_NY+j];
		}

	private:
		T* data;
		int _NX;
		int _NY;
};

template<class T>
class Array3D {
	public:
		Array3D() {
			this->data = NULL;
		}
		~Array3D() {
			if (this->data != NULL) {
				delete [] this->data;
			}
		}

		void init(const int NX, const int NY, const int NZ) {
			this->_NX = NX;
			this->_NY = NY;
			this->_NZ = NZ;
			if (this->data == NULL) {
				this->data = new T[this->_NX*this->_NY*this->_NZ];
			}
		}

		T operator()(const int i, const int j, const int k) {
			return this->data[i*this->_NY*this->_NZ+j*this->_NZ+k];
		}
		T& assign(const int i, const int j, const int k) {
			return this->data[i*this->_NY*this->_NZ+j*this->_NZ+k];
		}

	private:
		T* data;
		int _NX;
		int _NY;
		int _NZ;
};

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
