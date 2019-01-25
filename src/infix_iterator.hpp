//
//   Copyright 2019 University of Sheffield
//
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License.

#ifndef INFIX_ITERATOR_HPP
#define INFIX_ITERATOR_HPP

#include <iterator>
#include <ostream>
#include <string>

template <class T, class charT = char, class traits = std::char_traits<charT>>
class infix_ostream_iterator
    : public std::iterator<std::output_iterator_tag, void, void, void, void> {
  std::basic_ostream<charT, traits> *os;
  std::basic_string<charT> delimiter;
  std::basic_string<charT> real_delim;

public:
  typedef charT char_type;
  typedef traits traits_type;
  typedef std::basic_ostream<charT, traits> ostream_type;

  infix_ostream_iterator(ostream_type &s) : os(&s)
  {
  }

  infix_ostream_iterator(ostream_type &s, charT const *d) : os(&s), real_delim(d)
  {
  }

  infix_ostream_iterator<T, charT, traits> &operator=(T const &item)
  {
    *os << delimiter << item;
    delimiter = real_delim;
    return *this;
  }

  infix_ostream_iterator<T, charT, traits> &operator*()
  {
    return *this;
  }

  infix_ostream_iterator<T, charT, traits> &operator++()
  {
    return *this;
  }

  infix_ostream_iterator<T, charT, traits> &operator++(int)
  {
    return *this;
  }
};
#endif // INFIX_ITERATOR_HPP
