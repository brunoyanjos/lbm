#include "lbm_tuner.cuh"
#include <cstddef> // size_t
#include <cstdlib> // std::abs

__host__ dim3 find_optimal_block(size_t max_lattices)
{
    constexpr size_t MAX_THREADS = 1024;

    if (max_lattices == 0)
    {
        return dim3(1, 1, 1);
    }

    int best_x = 1, best_y = 1;
    size_t best_area = 1;
    int best_balance = 0;

    for (int bx = 32; bx >= 1; bx >>= 1)
    {
        for (int by = 32; by >= 1; by >>= 1)
        {
            const size_t area = static_cast<size_t>(bx) * static_cast<size_t>(by);

            if (area > max_lattices)
                continue;
            if (area > MAX_THREADS)
                continue;

            const int balance = std::abs(bx - by);

            // Critério:
            // 1) maior área
            // 2) mais "quadrado" (menor |bx-by|)
            // 3) se ainda empatar, prefira maior bx (tende a melhorar coalescing em x, se x for contíguo)
            if (area > best_area ||
                (area == best_area && balance < best_balance) ||
                (area == best_area && balance == best_balance && bx > best_x))
            {
                best_x = bx;
                best_y = by;
                best_area = area;
                best_balance = balance;
            }
        }
    }

    return dim3(best_x, best_y, 1);
}
