#include "geometry_outputs.cuh"

// OBS: não inclui active_geometry.cuh
// Aqui a seleção é só por macros LBM_GEOM_* definidas no compile.sh.

#if defined(LBM_GEOM_ANNUL)
#include "geometries/annul/output_annul_profile_bin.cuh"
#include "geometries/annul/output_annul_pressure_bin.cuh"
#endif

#if defined(LBM_GEOM_SQUARE_CAVITY)
#include "geometries/square_cavity/output_centerline_bin.cuh"
#endif

#if defined(LBM_GEOM_CHANNEL)
#include "geometries/poiseuille/output_channel_profile_bin.cuh"
#endif

#if defined(LBM_GEOM_COUETTE)
#include "geometries/couette/output_couette_profile_bin.cuh"
#endif

#if defined(LBM_GEOM_JET)
#include "geometries/jet/output_jet_centerline_bin.cuh"
#include "geometries/jet/output_jet_sections_bin.cuh"
#endif

namespace io
{
    void geometry_step_outputs(const LBMState &, const CudaConfig &, int, const std::string &, const DomainTags &)
    {
        // Por padrão vazio (pra não pesar IO). Depois você pode ligar por geometria se quiser.
    }

    void geometry_final_outputs(const LBMState &state,
                                const CudaConfig &cfg,
                                int t_end,
                                const std::string &out_dir,
                                const DomainTags &tags)
    {
        (void)cfg;
        (void)tags;

#if defined(LBM_GEOM_ANNUL)
        write_annul_profile(state, t_end, out_dir, tags);
        write_annul_pressure(state, t_end, out_dir, tags);
#endif

#if defined(LBM_GEOM_SQUARE_CAVITY)
        write_centerline_profiles(state, t_end, out_dir);
#endif

#if defined(LBM_GEOM_CHANNEL)
        write_channel_profile(state, t_end, out_dir, tags);
#endif

#if defined(LBM_GEOM_COUETTE)
        write_couette_profile(state, t_end, out_dir, tags);
#endif

#if defined(LBM_GEOM_JET)
        write_jet_centerline(state, t_end, out_dir, tags);
        write_jet_sections(state, t_end, out_dir, tags);
#endif

        // Outras geometrias entram aqui
    }
}