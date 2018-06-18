#ifndef IMAGE_H
#define IMAGE_H

#include <cstdlib>
#include <fstream>
#include <iostream>

template <typename T> class Image {
public:
  Image() : data_(NULL) { size_[0] = size_[1] = size_[2] = 0; }

  Image(size_t nx, size_t ny, size_t nz) {
    size_[0] = nx;
    size_[1] = ny;
    size_[2] = nz;
    // data_ = new T[nx * ny * nz];
    data_ = (T *)malloc(nx * ny * nz * sizeof(T));
    // std::cout << "Image() " << this << std::endl;
  }

  template <typename U> Image<T> &operator=(Image<U> &other) {
    reinitialize(other.size(0), other.size(1), other.size(2));
    for (size_t i = 0; i < length(); ++i)
      data_[i] = static_cast<T>(other[i]);
    return *this;
  }

  template <typename U> Image<T> &operator=(U *src) {
    for (size_t i = 0; i < length(); ++i)
      data_[i] = static_cast<U>(src[i]);
    return *this;
  }

  ~Image() {
    if (data_)
      free(data_);
  }

  inline T &operator[](size_t i) { return data_[i]; }
  inline const T &operator[](size_t i) const { return data_[i]; }

  inline T &operator()(size_t i, size_t j, size_t k) {
    return data_[i + size_[0] * j + size_[0] * size_[1] * k];
  }

  inline const T &operator()(size_t i, size_t j, size_t k) const {
    return data_[i + size_[0] * j + size_[0] * size_[1] * k];
  }

  inline size_t size(size_t i) const { return size_[i]; }

  inline size_t length() const { return size_[0] * size_[1] * size_[2]; }

  void load(const std::string &filename);
  void dump(const std::string &filename) const;

private:
  Image<T> &reinitialize(size_t nx, size_t ny, size_t nz) {
    if (data_ != NULL) {
      delete[] data_;
    }
    size_[0] = nx;
    size_[1] = ny;
    size_[2] = nz;
    data_ = new T[nx * ny * nz];

    return *this;
  }

private:
  size_t size_[3];
  T *data_;
};

template <typename T> void Image<T>::load(const std::string &filename) {
  std::ifstream input;
  input.open(filename.c_str(), std::ios::binary);
  if (!input.is_open()) {
    std::cout << "Error: Cannot open file" << std::endl;
  }
  std::cout << "Reading " << filename << std::endl;
  input.read(reinterpret_cast<char *>(data_),
             size_[0] * size_[1] * size_[2] * sizeof(T));
  input.close();
}

template <typename T> void Image<T>::dump(const std::string &filename) const {

  std::ofstream output;
  output.open(filename.c_str(), std::ios::binary);

  std::cout << "Writing " << filename << std::endl;
  output.write(reinterpret_cast<char *>(data_),
               size_[0] * size_[1] * size_[2] * sizeof(T));
  output.close();
}

#endif /* end of include guard: IMAGE_H */
