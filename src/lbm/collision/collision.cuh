#pragma once

#include "../../core/types.cuh"
#include "../../geometries/active_geometry.cuh"
#include "../moment/moment_values.cuh"
#include "../layout/collision_moment.cuh"
#include "stored_equilibrium_value.cuh"
#include "../meta/for_each_id.cuh"

template <class MomentList>
struct CollisionOperator
{
    MomentValues<MomentList> &M;
    real_t omega;
    real_t one_minus_omega;

    template <auto Id>
    __device__ __forceinline__ void operator()()
    {
        real_t &m = M.template get<Id>();
        const real_t meq = stored_equilibrium_value<Id>(M);
        m = one_minus_omega * m + omega * meq;
    }
};

template <class MomentList>
__device__ __forceinline__ void moment_space_collision(MomentValues<MomentList> &M)
{
    const real_t omega = Geometry::OMEGA;
    const real_t one_minus_omega = r::one - omega;

    for_each_id<CollisionMomentList>(CollisionOperator<MomentList>{M, omega, one_minus_omega});
}