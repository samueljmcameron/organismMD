add_executable(cellsplit main.cpp cellsplit.cpp input.cpp utility.cpp domain.cpp comm_brick.cpp atom.cpp read_atoms.cpp group.cpp fix.cpp fix_organism_birthdeath.cpp integrate.cpp dump.cpp compute.cpp)

target_link_libraries(cellsplit PUBLIC MPI::MPI_CXX Eigen3::Eigen)

