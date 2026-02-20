#pragma once
#include "domain_tags.cuh"

// Cavidade quadrada:
// - domínio é [0..NX-1]x[0..NY-1]
// - bordas são CONTORNO (não sólido)
// - interior é bulk
// valid_mask: bit i = 1 se (x+cx(i), y+cy(i)) está dentro do domínio.
// bit 0 sempre 1 (direção rest) para todos os nós.
void build_cavity_square_tags(DomainTags &T);
