#pragma once

#include <string>
#include <sstream>
#include <cxxabi.h>
#include <cstdlib>

using std::string;

namespace meta {

// is_same
template <typename T1, typename T2>
struct is_same {
    const static bool value = false;
};

template <typename T>
struct is_same<T,T> {
    const static bool value = true;
};

// lexical_cast
template <typename T1, typename T2, bool Same>
struct lexical_cast_helper {
    static T1 cast(const T2& src) {
        T1 ret;
        std::stringstream ss;
        ss << src;
        ss >> ret;
        return ret;
    }
};

template <typename T>
struct lexical_cast_helper<T,T, true> {
    static T cast(const T& src) {
        return src;
    }
};

template <typename T>
struct lexical_cast_helper<std::string,T, false> {
    static std::string cast(const T& src) {
        std::ostringstream ss;
        ss << src;
        return ss.str();
    }
};

template <typename T>
struct lexical_cast_helper<T,std::string, false> {
    static T cast(const std::string& src) {
        T ret;
        std::istringstream ss(src);
        ss >> ret;
        return ret;
    }
};

template <typename T1, typename T2>
T1 lexical_cast(const T2& src) {
    return lexical_cast_helper<T1, T2, is_same<T1,T2>::value>::cast(src);
};


// typename
static inline std::string demangle(const std::string &name)
{
  int status=0;
  char *p=abi::__cxa_demangle(name.c_str(), 0, 0, &status);
  std::string ret(p);
  free(p);
  return ret;
}

template <class T>
std::string readable_typename()
{
  return demangle(typeid(T).name());
}

template<>
std::string readable_typename<std::string>()
{
    return "string";

}

} // namespace meta
