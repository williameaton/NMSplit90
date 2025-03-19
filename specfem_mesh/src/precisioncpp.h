#ifndef PRECISION_H
#define PRECISION_H

// Define precision toggle (comment/uncomment the desired type)
//#define USE_DOUBLE_PRECISION  // Comment this to switch to float

#ifdef USE_DOUBLE_PRECISION
    using CPPCUSTOM_REAL = double;
#else
    using CPPCUSTOM_REAL = float;
#endif

#endif  // PRECISION_H
