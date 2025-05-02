#ifndef GSMINRES_C_API_H
#define GSMINRES_C_API_H

#include <stddef.h>
#include <complex.h>

#ifdef __cplusplus
//extern "C" {
#endif

// Opaque solver handle (内部のC++クラスを隠蔽)
typedef void* gsminres_handle;

//
gsminres_handle gsminres_create(size_t, size_t m);

//
void gsminres_destroy(gsminres_handle handle);

//
void gsminres_initialize(gsminres_handle handle,
                         double complex* x,
                         const double complex* b,
                         double complex* w,
                         const double complex* sigma,
                         double threshold);

//
int gsminres_update(gsminrea_handle handle,
                    double complex* x,
                    size_t n, size_t m);


#endif // GSMINRES_C_API_H
