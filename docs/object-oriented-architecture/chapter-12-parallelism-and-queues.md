# 章节 12：并行、随机数与粒子队列管理

本章节说明并行执行相关结构如何利用粒子成员，涵盖随机数生成、粒子银行与消息传递。

## 随机数生成

- `ParticleData::rng_seed`、`stream`：线性同余随机数生成器的状态。
- `RNG Particle::rng()`：包装 `RandomLCG`，每次事件调用时更新 `stream` 并返回伪随机数。【F:include/openmc/particle.h†L40-L68】【F:include/openmc/random_lcg.h†L27-L166】
- `RandomLCG` 内部成员 `uint64_t state_`、`multiplier_`、`increment_` 决定数列，`set_seed` 与 `advance` 提供可重复性控制。

## 粒子银行结构

- `BankedParticle` 成员：`Position r`、`Direction u`、`double E`、`double wgt`、`int type`、`uint64_t rng_seed`，与 `ParticleData` 对齐。【F:include/openmc/bank.h†L31-L205】
- `GlobalBank`：维护 `std::vector<BankedParticle> active` 与 `secondary`，通过 `push_back`、`merge_secondary` 组织待处理粒子。【F:src/bank.cpp†L34-L285】
- `ParticleData::n_secondary` 与 `secondary_bank_` 指针控制一次碰撞生成的次级粒子数量与写入位置。【F:include/openmc/particle_data.h†L120-L204】

## 消息传递与负载均衡

- `message_passing.cpp` 中的 `WorkQueue` 结构保存 `std::deque<ParticleData>` 或句柄，用于跨线程/进程调度。
- `WorkBalance` 成员 `int max_particles`、`int threshold` 控制任务分配；`steal_particles` 使用 MPI/NCCL 接口在进程间移动 `BankedParticle`。【F:src/message_passing.cpp†L38-L327】

## 设计重点

- 粒子状态可直接复制到银行或消息缓冲区，减少序列化开销。
- 随机数种子作为成员随粒子传递，保证跨进程移动后仍可复现历史。
- 工作窃取与银行合并操作均依赖 `ParticleData` 中的 `alive`、`wgt`、`n_secondary` 等字段判定粒子是否需要继续推进。
