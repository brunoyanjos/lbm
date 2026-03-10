#include "domain_tags.cuh"
#include "../../geometries/active_geometry.cuh"
#include "../../core/cuda_utils.cuh"

DomainTags domain_tags_allocate()
{
    DomainTags T{};
    T.N = static_cast<size_t>(Geometry::NX) * static_cast<size_t>(Geometry::NY);
    T.bytes_valid = T.N * sizeof(uint32_t);
    T.bytes_node = T.N * sizeof(uint8_t);

    T.h_valid = static_cast<uint32_t *>(std::malloc(T.bytes_valid));
    T.h_node = static_cast<uint8_t *>(std::malloc(T.bytes_node));
    if (!T.h_valid || !T.h_node)
    {
        std::free(T.h_valid);
        std::free(T.h_node);
        throw std::bad_alloc();
    }

    CUDA_CHECK(cudaMalloc(&T.d_valid, T.bytes_valid));
    CUDA_CHECK(cudaMalloc(&T.d_node, T.bytes_node));

    CUDA_CHECK(cudaMemset(T.d_valid, 0, T.bytes_valid));
    CUDA_CHECK(cudaMemset(T.d_node, 0, T.bytes_node));

    return T;
}

void domain_tags_free(DomainTags &T)
{
    std::free(T.h_valid);
    std::free(T.h_node);
    T.h_valid = nullptr;
    T.h_node = nullptr;

    if (T.d_valid)
        CUDA_CHECK(cudaFree(T.d_valid));
    if (T.d_node)
        CUDA_CHECK(cudaFree(T.d_node));
    T.d_valid = nullptr;
    T.d_node = nullptr;

    T.N = 0;
    T.bytes_valid = 0;
    T.bytes_node = 0;
}
