# scope_guard

[![Build Status](https://github.com/ricab/scope_guard/actions/workflows/test.yml/badge.svg)](https://github.com/ricab/scope_guard/actions/workflows/test.yml)
[![Build status](https://img.shields.io/appveyor/ci/ricab/scope-guard.svg)](https://ci.appveyor.com/project/ricab/scope-guard/branch/main)
[![codecov](https://img.shields.io/codecov/c/github/ricab/scope_guard.svg)](https://codecov.io/gh/ricab/scope_guard)

[![GitHub release](https://img.shields.io/github/v/release/ricab/scope_guard?include_prereleases)](https://github.com/ricab/scope_guard/releases)
[![Vcpkg](https://img.shields.io/vcpkg/v/scope-guard)](https://vcpkg.link/ports/scope-guard)
[![semver](https://img.shields.io/badge/semver-2.0.0-blue.svg)](https://semver.org/spec/v2.0.0.html)
[![GitHub license](https://img.shields.io/github/license/ricab/scope_guard.svg)](https://github.com/ricab/scope_guard/blob/main/LICENSE)

[![GitHub stars](https://img.shields.io/github/stars/ricab/scope_guard.svg?style=social&label=Stars)](https://github.com/ricab/scope_guard)

A public, general, simple, and fast C++11 scope guard that
defends against implicitly ignored returns and optionally enforces `noexcept`
at compile time (in C++17), all in a SFINAE-friendly maner.

##### TLDR

Get it [here](scope_guard.hpp).
Usage is simple:

```c++
#include "scope_guard.hpp"
...
{
  ...
  auto guard = sg::make_scope_guard(my_callback);
  ...
} // my_callback is invoked at this point
```

## Introduction

A scope guard is an object that employs RAII to execute a
provided callback when leaving scope, be it through a _fall-through_, a return,
or an exception. That callback can be a function, a function pointer, a functor,
a lambda, a bind result, a std::function, a reference to any of these, or any
other callable, as long as it respects a few [preconditions](docs/precond.md)
&ndash; most of which are enforced during compilation, the rest being hopefully
intuitive.

All necessary code is provided in a [single header](scope_guard.hpp) (the
remaining files are only for testing and documentation).

#### Acknowledgments

The concept of "scope guard" was [proposed](http://drdobbs.com/184403758)
by Andrei Alexandrescu and Petru Marginean and it is well known in the
C++ community. It has been proposed for standardization (see
[N4189](http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2014/n4189.pdf))
but is still not part of the standard library, as of March 2018.

#### Why

While there are several implementations available, I did not find any with all
the characteristics I aimed for here - safe, tested, documented, public domain, thin wrapping, general, standalone, simple interface... (see feature list below).

## Features

#### Main features
- [x] &ge; C++11
- [x] Reduced interface
- [x] Thin callback wrapping: no added `std::function` or virtual table
penalties
- [x] General: accepts any callable that respects a few
[preconditions](docs/precond.md)
- [x] No implicitly ignored return (details [here](docs/precond.md#void-return))
- [x] Option to enforce `noexcept` in C++17
(details [here](docs/interface.md#compilation-option-sg_require_noexcept_in_cpp17))
- [x] Modern exception specifications (`noexcept` with conditions when
necessary)
- [x] SFINAE friendliness (see [here](docs/design.md#sfinae-friendliness))

#### Other characteristics
- [x] No dependencies to use (besides &ge;C++11 compiler and standard library)
- [x] No macros to make guard &ndash; just write explicit lambda or bind or
what have you
- [x] Extensively tested, with both
[compile time tests](compile_time_tests.cpp) and
[runtime-tests](catch_tests.cpp)
- [x] Carefully documented (adhering to
[RFC2119](https://tools.ietf.org/html/rfc2119))
- [x] Standalone header that can be directly dumped into any project
- [x] Unlicense'd
- [x] `snake_case` style

#### Notes

<sup>_The key words "MUST", "MUST NOT", "REQUIRED", "SHALL", "SHALL
NOT", "SHOULD", "SHOULD NOT", "RECOMMENDED",  "MAY", and "OPTIONAL" in this
document are to be interpreted as described in RFC 2119._</sup>

#### Issues

Bug reports and suggestions are welcome. If you find that something is incorrect
or could be improved, feel free to open an issue.

## Setup

Setup consists merely of making the [header file](scope_guard.hpp) available to
the compiler. That can be achieved by any of the following options:

- placing it directly in the client project's include path
- placing it in a central include path that is known to the compiler
- placing it in an arbitrary path and configuring the compiler to include that
path

The preprocessor definition `SG_REQUIRE_NOEXCEPT_IN_CPP17` MAY be provided
to the compiler. The effect of this option is explained
[here](docs/interface.md#compilation-option-sg_require_noexcept_in_cpp17).

## Further documentation

#### Client interface

The client interface is documented in detail [here](docs/interface.md).

#### Preconditions in detail

Callback preconditions are explained [here](docs/precond.md).

#### Design choices and concepts

Design choices and concepts are discussed [here](docs/design.md).

#### Tests

Instructions on how to run the tests are [here](docs/tests.md).
