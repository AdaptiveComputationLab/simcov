#include "utils.hpp"

using namespace upcxx_utils;

using std::min;
using std::string;
using std::string_view;
using std::to_string;

int pin_thread(pid_t pid, int cid) {
  cpu_set_t cpu_set;
  CPU_ZERO(&cpu_set);
  CPU_SET(cid, &cpu_set);
  if (sched_setaffinity(pid, sizeof(cpu_set), &cpu_set) == -1) {
    if (errno == 3) WARN("%s, pid: %d", strerror(errno), pid);
    return -1;
  }
  return 0;
}

void dump_single_file(const string &fname, const string &out_str) {
  auto sz = out_str.length();
  upcxx::atomic_domain<size_t> ad({upcxx::atomic_op::fetch_add, upcxx::atomic_op::load});
  upcxx::global_ptr<size_t> fpos = nullptr;
  // always give rank 0 the first chunk so it can write header info if needed
  if (!upcxx::rank_me()) fpos = upcxx::new_<size_t>(sz);
  fpos = upcxx::broadcast(fpos, 0).wait();
  size_t my_fpos = 0;
  if (upcxx::rank_me()) my_fpos = ad.fetch_add(fpos, sz, std::memory_order_relaxed).wait();
  // wait until all ranks have updated the global counter
  upcxx::barrier();
  int fileno = -1;
  size_t fsize = 0;
  if (!upcxx::rank_me()) {
    fsize = ad.load(fpos, std::memory_order_relaxed).wait();
    // rank 0 creates the file and truncates it to the correct length
    fileno = open(fname.c_str(), O_WRONLY | O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
    if (fileno == -1) WARN("Error trying to create file ", fname, ": ", strerror(errno), "\n");
    if (ftruncate(fileno, fsize) == -1)
      WARN("Could not truncate ", fname, " to ", fsize, " bytes\n");
  }
  upcxx::barrier();
  ad.destroy();
  // wait until rank 0 has finished setting up the file
  if (rank_me()) fileno = open(fname.c_str(), O_WRONLY, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
  if (fileno == -1) WARN("Error trying to open file ", fname, ": ", strerror(errno), "\n");
  auto bytes_written = pwrite(fileno, out_str.c_str(), sz, my_fpos);
  close(fileno);
  if (bytes_written != sz)
    DIE("Could not write all ", sz, " bytes; only wrote ", bytes_written, "\n");
  upcxx::barrier();
  auto tot_bytes_written = upcxx::reduce_one(bytes_written, upcxx::op_fast_add, 0).wait();
  SLOG_VERBOSE("Successfully wrote ", get_size_str(tot_bytes_written), " to ", fname, "\n");
}

std::shared_ptr<Random> _rnd_gen;
