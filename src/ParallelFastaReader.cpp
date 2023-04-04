// Created by Saliya Ekanayake on 1/6/19.

#include <iostream>
#include <fstream>
#include <memory>
#include <cassert>
#include <limits>
#include <unordered_set>
#include "../include/ParallelFastaReader.hpp"
#include "../include/DistributedFastaData.hpp"

void ParallelFastaReader::read_fasta(const char *file, uint64_t overlap, int rank,
                                     int world_size, char *&buff,
                                     uint64_t &l_start, uint64_t &l_end) {

  MPI_File f;

  int err = MPI_File_open(MPI_COMM_WORLD, file,
                          MPI_MODE_RDONLY, MPI_INFO_NULL, &f);
  if (err) {
    if (rank == 0) fprintf(stderr, "Couldn't open file %s\n", file);
    MPI_Finalize();
    exit(2);
  }

  /* The idea:
   * if sequence names and contents are roughly equal
   * then dividing the file row_size by the number of ranks
   * should land close to l_seq_count sequences. We'll fix
   * any mismatch later
   */

  /* Thanks, Rob Lathem for making my life a bit easy here
   * https://stackoverflow.com/a/12942718 */
  MPI_Offset g_start;
  uint64_t l_chunk_size;

  /* read in relevant chunk of file into "chunk",
   * which starts at location in the file g_start
   * and has row_size l_chunk_size
   */

  MPI_Offset g_end;
  MPI_Offset file_size = 0;

  std::cout << "Fsize " << file_size << std::endl; 
  /* figure out who reads what */
  MPI_File_get_size(f, &file_size);
  file_size--;  /* get rid of text file eof */
  std::cout << "Fsize " << file_size << std::endl; 
  l_chunk_size = static_cast<uint64_t >(file_size / world_size);

  g_start = rank    * l_chunk_size;
  g_end   = g_start + l_chunk_size - 1;

  if (rank == world_size - 1)
    g_end = file_size - 1;

  /* add overlap to the end of everyone's chunk except last proc */
  if (rank != world_size - 1)
    g_end += overlap;

  l_chunk_size = static_cast<uint64_t >(g_end - g_start + 1); // >> 2 l_chunk_size  8793948

  std::cout << "ChunkSize " << l_chunk_size << " rank " << rank << std::endl; 

  buff = new char[l_chunk_size+1];

  /* everyone reads in their part */

  // TODO - Saliya: fix if l_chunk_size > int max. Read multiple times
  // assert(l_chunk_size <= std::numeric_limits<int>::max());

  // MPI_File_read_at_all(f, g_start, buff, static_cast<int>(l_chunk_size),
  //                      MPI_CHAR, MPI_STATUS_IGNORE);
  for (size_t t = 0; t < l_chunk_size; t+= std::numeric_limits<int>::max()) {
    auto readlen = std::min((size_t)l_chunk_size - t, (size_t)std::numeric_limits<int>::max());
    MPI_File_read_at_all(f, g_start + t, &buff[t], static_cast<int>(readlen),
                       MPI_CHAR, MPI_STATUS_IGNORE);
  }
                    
  buff[l_chunk_size] = '\0';

  /*
   * everyone calculate what their start and end *really* are by going
   * from the first newline after start to the first newline after the
   * overlap region starts (eg, after end - overlap + 1)
   */

  l_start = 0, l_end = l_chunk_size - 1;

  if (rank != 0)
  {
    while (buff[l_start] != '>') l_start++;
  }

  if (rank != world_size - 1)
  {
    l_end -= overlap;
    while (buff[l_end] != '>') l_end++;
    // minus 2 because we don't need '>' as well as the '\n' before that
    l_end -= 2;
  }

  /* GGGG: terminate program if error detected (https://stackoverflow.com/questions/10818740/gracefully-exit-with-mpi) */
  // int error = 0, therank;
  // if(l_end <= l_start)
  // {
  //   error = 1;
  //   therank = rank;
  // }

  if(l_end <= l_start)
  {
    std::cout<< "Error: Program terminated because FASTA is too small for " << world_size << " procs. Run with fewer procs than the number of input reads." << std::endl;

    MPI_Abort(MPI_COMM_WORLD, 1);
    /* No further code will execute */
    MPI_Finalize();
    exit(1);
  }
}