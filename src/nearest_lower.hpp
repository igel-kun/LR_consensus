
#pragma once

#include <vector>
#include <stack>

enum Direction : char {DIR_LEFT = 1, DIR_RIGHT = -1};

// preprocessing a vector such as to allow querying the nearest element on the left (right) index with lower (higher) score
// per default: find the nearest element on the left with less score than the query element
// NOTE: Container should support random access by operator[](unsigned)
template<class Container, Direction dir = DIR_LEFT, class Compare = std::less<typename Container::value_type> >
class nearest_lower{
  std::vector<ssize_t> target_index; // vector storing the actual information

  void preprocess(const Container& container, const Compare& cmp)
  {
    target_index.resize(container.size());
    // stack of indices of container
    std::stack<ssize_t> prev_indices;
    // go from left to right, keeping the stack of "visible" values
    // (a value is "visible" from us if there is no other value of at least this score between them and us)
    const ssize_t lower_bound = (dir == DIR_LEFT ? 0 : container.size() - 1);
    const ssize_t upper_bound = (dir == DIR_LEFT ? container.size() : - 1);
    for(ssize_t i = lower_bound; i != upper_bound; i += dir){
      const auto& current = container[i];
      // pop everything that is being shadowed by us from the stack
      while(!prev_indices.empty() && cmp(current, container[prev_indices.top()]))
        prev_indices.pop();
      // the next item on the stack is the nearest with higher/lower score
      target_index[i] = prev_indices.empty() ? (ssize_t)(-1) : prev_indices.top();
      // add us to the stack
      prev_indices.push(i);
    }
  }

public:

  nearest_lower(const Container& _container, const Compare& _cmp = Compare())
  {
    preprocess(_container, _cmp);
  }

  // return the index in _container that holds the nearest element with higher/lower score
  ssize_t query(const unsigned index) const { return target_index[index]; }
};
