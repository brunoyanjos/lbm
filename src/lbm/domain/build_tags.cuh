#pragma once
#include "app/run_context.cuh"
#include "domain_tags.cuh"

#include <vector>

void build_tags(DomainTags &T);
[[nodiscard]] __host__ std::vector<DomainTags> allocate_partition_tags(const app::RunContext &ctx);
__host__ void free_partition_tags(std::vector<DomainTags> &tags, const app::RunContext &ctx);
