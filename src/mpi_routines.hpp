//
//   Copyright 2019 University of Sheffield
//
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License.

#include <exception>
#include <sstream>
#include <vector>

#include <mpi.h>

template <typename T>
std::vector<std::vector<T>> p2p_vecscatter(const std::vector<std::vector<T>> &src_data,
                                           MPI_Comm comm)
{
  size_t objsize = sizeof(T);
  int comm_size, rank;
  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &rank);

  if (src_data.size() != static_cast<size_t>(comm_size))
  {
    std::ostringstream errss;
    errss << "Vector length must be equal to communicator size (recieved " << src_data.size()
          << ", expected " << comm_size;
    throw std::length_error(errss.str());
  }

  // First, fire off ISSENDS to each other rank
  std::vector<MPI_Request> send_reqs(comm_size);
  for (int irank = 0; irank < comm_size; irank++)
  {
    if (irank == rank)
    {
      send_reqs[rank] = MPI_REQUEST_NULL; // This is likely unnecessary
      continue;
    }
    MPI_Issend(src_data[irank].data(), src_data[irank].size() * objsize, MPI_BYTE, irank, 0, comm,
               &send_reqs[irank]);
  }

  // Loop over comm size because we should receive as many messages as we send. If not, everything
  // is terrible and we deserve to deadlock
  std::vector<std::vector<T>> recvd_data(comm_size);
  for (int idx = 0; idx < comm_size; idx++)
  {
    if (idx == rank)
    {
      recvd_data[rank] = src_data[rank]; // Copy from src to recv for own rank
      continue;                          // Skip one because we did in sending
    }
    // Now check for a message and get the length
    MPI_Status status;
    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &status);
    int length;
    MPI_Get_count(&status, MPI_BYTE, &length);

    // Allocate needed elements in appropriate vector
    recvd_data[status.MPI_SOURCE].resize(length / objsize);

    // Now do the actual receive
    MPI_Recv(recvd_data[status.MPI_SOURCE].data(), length, MPI_BYTE, status.MPI_SOURCE,
             status.MPI_TAG, comm, &status);
  }

  // Wait for all ops to complete, avoids potential data corruption from caller changing src_data
  MPI_Barrier(comm);

  return recvd_data;
}
