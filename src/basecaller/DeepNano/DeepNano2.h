#include <cstdarg>
#include <cstdint>
#include <cstdlib>
#include <new>

struct Caller;

using SgemmJitKernelT = void(*)(void *arg1, float *arg2, float *arg3, float *arg4);

extern "C" {

char *call_raw_signal(Caller *ptr, const float *raw_data, size_t size);

Caller *create_caller(const char *net_type,
                      const char *path,
                      uintptr_t beam_size,
                      float beam_cut_threshold);

extern uint32_t mkl_cblas_jit_create_sgemm(void **JITTER,
                                           uint32_t layout,
                                           uint32_t transa,
                                           uint32_t transb,
                                           uintptr_t m,
                                           uintptr_t n,
                                           uintptr_t k,
                                           float alpha,
                                           uintptr_t lda,
                                           uintptr_t ldb,
                                           float beta,
                                           uintptr_t ldc);

extern SgemmJitKernelT mkl_jit_get_sgemm_ptr(const void *JITTER);

extern void vmsExp(int n, const float *a, float *y, long long mode);

} // extern "C"
