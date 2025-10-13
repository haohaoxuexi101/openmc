#include <mpi.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>

// 课程 6：MPI 非阻塞裂变银行流水线
// ---------------------------------
// 本课展示如何利用 MPI 的非阻塞通信 (Isend/Irecv) 构建裂变银行交换流水线。
// 每个进程维护两个缓冲区：一个发送给下游进程，一个从上游进程接收。
// 当粒子数量较大时，这种流水线可以与 OpenMC 主循环重叠通信与计算。

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  int rank = 0, size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const int prev = (rank - 1 + size) % size;
  const int next = (rank + 1) % size;

  std::vector<int> send_bank(8, rank * 10);
  std::vector<int> recv_bank(8, -1);

  MPI_Request send_req = MPI_REQUEST_NULL;
  MPI_Request recv_req = MPI_REQUEST_NULL;

  MPI_Irecv(recv_bank.data(), recv_bank.size(), MPI_INT, prev, 99, MPI_COMM_WORLD,
            &recv_req);
  MPI_Isend(send_bank.data(), send_bank.size(), MPI_INT, next, 99, MPI_COMM_WORLD,
            &send_req);

  // 模拟与通信重叠的计算工作
  std::this_thread::sleep_for(std::chrono::milliseconds(20));

  MPI_Wait(&send_req, MPI_STATUS_IGNORE);
  MPI_Wait(&recv_req, MPI_STATUS_IGNORE);

  std::cout << "Rank " << rank << " 收到来自 " << prev << " 的裂变银行样本: ";
  for (int value : recv_bank) {
    std::cout << value << ' ';
  }
  std::cout << "\n";

  MPI_Finalize();
  return 0;
}

