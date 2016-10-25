
//! file vector2d.hpp
/** This is a 2 dimentional vector (aka matrix) that is slightly faster than
 * vector<vector<T>> and slightly more convenient than keeping track of the
 * index offsets using vector<T>
 **/

#ifndef VECTOR2D_HPP
#define VECTOR2D_HPP

#include <cassert>
#include <vector>

namespace std{

  template<typename T>
  class _vector2d : public vector<T>
  {
  protected:
    //! translate 2d coordinates into 1d coordinates
    virtual size_t linearize(const pair<size_t,size_t>& coords) const = 0;

  public:
    virtual ~_vector2d(){}
    //! NOTE: coordinates are always (col, row) as in (x, y)
    /** also: why in the hell does C++ force operator[] to take exactly one argument??? **/
    T& operator[](const pair<size_t, size_t>& coords)
    {
      return vector<T>::operator[](linearize(coords));
    }
    const T& operator[](const pair<size_t, size_t>& coords) const
    {
      return vector<T>::operator[](linearize(coords));
    }
    const T& at(const pair<size_t, size_t>& coords) const
    {
      return vector<T>::at(linearize(coords));
    }
    void resize(const size_t cols, const size_t rows)
    {
      vector<T>::resize(linearize(cols - 1, rows - 1) + 1);
    }
    virtual pair<size_t, size_t> size() const = 0;
  };// class vector2d



  template<typename T>
  class vector2d : public _vector2d<T>
  {
  protected:
    size_t columns;

    size_t linearize(const pair<size_t,size_t>& coords) const
    {
      return coords.second * columns + coords.first;
    }
  public:
    vector2d(): vector<T>() {}

    vector2d(const size_t cols, const size_t rows, const T& element = T()):
      columns(cols)
    {
      // reserve size such that we can access the last index (cols-1,rows-1)
      vector<T>::resize(linearize({cols - 1, rows - 1}) + 1, element);
    }

    void resize(const size_t cols, const size_t rows)
    {
      _vector2d<T>::resize(cols, rows);
      columns = cols;
    }
    pair<size_t,size_t> size() const
    {
      return make_pair(columns, vector<T>::size() / columns);
    }
  };


  //! a symmetric version of the 2d vector that only uses half the space
  template<typename T>
  class symmetric_vector2d : public _vector2d<T>
  {
  protected:
    //! translate 2d coordinates into 1d coordinates
    size_t linearize(const pair<size_t,size_t>& coords) const
    {
      if(coords.first <= coords.second)
        return (coords.second * (coords.second + 1)) / 2  + coords.first;
      else
        return (coords.first * (coords.first + 1)) / 2  + coords.second;
    }
    
  public:
    symmetric_vector2d(): vector<T>() {}

    symmetric_vector2d(const size_t rows_and_cols, const T& element = T())
    {
      vector<T>::resize(linearize({rows_and_cols, rows_and_cols}), element);
    }

    pair<size_t,size_t> size() const
    {
      const size_t root = (size_t)(sqrt(vector<T>::size() * 2 + 0.25) - .5);
      return make_pair(root, root);
    }
  };// class vector2d



}// namespace


#endif

