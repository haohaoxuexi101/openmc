#include <mpi.h>
#include <iostream>

// 课程 2：MPI 裂变银行归并
// 以最小示例展示如何在 MPI 任务间汇总粒子数量。

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  int rank = 0, size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int local_bank = 10 + rank; // 每个进程的裂变粒子数
  int global_bank = 0;

  MPI_Allreduce(&local_bank, &global_bank, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if (rank == 0) {
    std::cout << "总裂变粒子数 = " << global_bank << " (来自 " << size
              << " 个进程)\n";
  }

  MPI_Finalize();
  return 0;
}
