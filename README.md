Unofficial mirror of NOVAS: Naval Observatory Vector Astrometry Subroutines

### Status

[![GitHub release](https://img.shields.io/github/release/jacobwilliams/NOVAS.svg)](https://github.com/jacobwilliams/NOVAS/releases/latest)
[![Build Status](https://github.com/jacobwilliams/NOVAS/actions/workflows/CI.yml/badge.svg)](https://github.com/jacobwilliams/NOVAS/actions)
[![codecov](https://codecov.io/gh/jacobwilliams/NOVAS/branch/master/graph/badge.svg)](https://codecov.io/gh/jacobwilliams/NOVAS)
[![last-commit](https://img.shields.io/github/last-commit/jacobwilliams/NOVAS)](https://github.com/jacobwilliams/NOVAS/commits/master)

### Compiling

A [Fortran Package Manager](https://github.com/fortran-lang/fpm) manifest file is included, so that the library and test cases can be compiled with FPM. For example:

```
fpm build --profile release
fpm test --profile release
```

To use `NOVAS` within your fpm project, add the following to your `fpm.toml` file:
```toml
[dependencies]
NOVAS = { git="https://github.com/jacobwilliams/NOVAS.git" }
```

Or, to use a specific version:

```toml
[dependencies]
NOVAS = { git="https://github.com/jacobwilliams/NOVAS.git", tag = "3.1.0" }
```

### Documentation

The latest API documentation can be found [here](http://jacobwilliams.github.io/NOVAS/). This was generated from the source code using [FORD](https://github.com/Fortran-FOSS-Programmers/ford).

### See also

 * Original source: https://www.usno.navy.mil/USNO/astronomical-applications/software-products/novas
 * [astro-fortran](https://github.com/jacobwilliams/astro-fortran) -- Modern Fortran implementations of standard models used in fundamental astronomy