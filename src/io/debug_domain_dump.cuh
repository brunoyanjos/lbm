#pragma once

#include <string>
#include "../lbm/domain/domain_tags.cuh"

namespace io
{
    // imprime mapa 0/1/D no stdout
    void debug_print_nodes(const DomainTags &tags);

    // imprime masks valid/in/out por célula (stdout)
    void debug_print_masks(const DomainTags &tags);

    // opcional: escreve um txt no out_dir (melhor que poluir stdout)
    void debug_write_masks_txt(const DomainTags &tags, const std::string &out_dir, const char *name);
}