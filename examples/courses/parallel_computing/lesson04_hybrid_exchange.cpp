#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <vector>

// 课程 4：混合并行数据交换
// 演示 MPI + OpenMP 的协作：线程收集局部统计，MPI 汇总。

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::vector<double> thread_sum(omp_get_max_threads(), 0.0);

  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    thread_sum[tid] = 1.0 + rank + tid * 0.1;
  }

  double local_total = 0.0;
  for (double v : thread_sum) local_total += v;

  double global_total = 0.0;
  MPI_Allreduce(&local_total, &global_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if (rank == 0) {
    std::cout << "混合并行总和 = " << global_total << "\n";
  }

  MPI_Finalize();
  return 0;
}
