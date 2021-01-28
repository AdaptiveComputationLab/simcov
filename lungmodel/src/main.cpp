// LungModel
// Author: Akil Andrews
#include <fcntl.h>
#include <unistd.h>
#include "lung.hpp"

int main(int argc, char** argv) {
  // Create a full or parital model of the human lung
  bool isFullModel = false;
  if (argc == 2) {
    std::string sel(argv[1]);
    if (sel == "true") {
      isFullModel = true;
    }
  }
  Lung lung(isFullModel);
  // Write alveolus epithileal cells
  std::string fname = "alveolus.dat";
  auto fileno =
      open(fname.c_str(), O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
  if (fileno == -1) {
    std::printf("Cannot create file %s\n", fname.c_str());
  }
  int bytes_written = 0;
  for (const int& a : lung.getAlveoliEpiCellIds()) {
    bytes_written += write(fileno, reinterpret_cast<const char*>(&a), sizeof(int));
  }
  if (bytes_written != sizeof(int) * (lung.getAlveoliEpiCellIds()).size()) {
    std::printf("Could not write all %ld bytes; only wrote %d\n",
                (lung.getAlveoliEpiCellIds()).size(), bytes_written);
  }
  close(fileno);
  // Write bronchiole epithileal cells
  fname = "bronchiole.dat";
  fileno = open(fname.c_str(), O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
  if (fileno == -1) {
    std::printf("Cannot create file %s\n", fname.c_str());
  }
  bytes_written = 0;
  for (const int& a : lung.getAirwayEpiCellIds()) {
    bytes_written += write(fileno, reinterpret_cast<const char*>(&a), sizeof(int));
  }
  if (bytes_written != sizeof(int) * (lung.getAirwayEpiCellIds()).size()) {
    std::printf("Could not write all %ld bytes; only wrote %d\n",
                (lung.getAirwayEpiCellIds()).size(), bytes_written);
  }
  close(fileno);
  return 0;
}
