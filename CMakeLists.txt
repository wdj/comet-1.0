#==============================================================================
#
# CMakeLists.txt
#
#==============================================================================

cmake_minimum_required(VERSION 2.8)

project("genomics_gpu")

#--------------------

SET(INCLUDE_DIRS
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/driver
  )

SET(LIB_SOURCES
  src/env.cc
  src/assertions.cc
  src/decomp_mgr.cc
  src/mirrored_buf.cc
  src/linalg.cc
  src/linalg_accel.cc
  src/metrics.cc
  src/checksum.cc
  src/vectors.cc
  src/vector_sums.cc
  src/comm_xfer_utils.cc
  src/compute_metrics_2way_block_nums.cc
  src/compute_metrics_2way_block_combine.cc
  src/compute_metrics_3way_block_nongpu.cc
  src/compute_metrics_3way_block_gpu.cc
  src/compute_metrics_3way_block.cc
  src/compute_metrics_2way.cc
  src/compute_metrics_3way.cc
  src/compute_metrics.cc
  )

SET(LIB_INCLUDES
  src/env.hh
  src/assertions.hh
  src/types.hh
  src/decomp_mgr.hh
  src/mirrored_buf.hh
  src/linalg.hh
  src/linalg_accel.hh
  src/metrics.hh
  src/metrics_2way_indexing.hh
  src/metrics_2way_accessors.hh
  src/metrics_3way_indexing.hh
  src/metrics_3way_accessors.hh
  src/checksum.hh
  src/vectors.hh
  src/vector_sums.hh
  src/comm_xfer_utils.hh
  src/compute_metrics_2way_block_nums.hh
  src/compute_metrics_2way_block_combine.hh
  src/compute_metrics_3way_block_nongpu.hh
  src/compute_metrics_3way_block_gpu.hh
  src/compute_metrics_3way_block.hh
  src/compute_metrics_2way.hh
  src/compute_metrics_3way.hh
  src/compute_metrics.cc
  )

#--------------------

if(TESTING) 
  set(C_AND_CXX_FLAGS "${C_AND_CXX_FLAGS} -isystem ${GTEST_DIR}/include -pthread -DTESTING")
endif(TESTING) 

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${C_AND_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${C_AND_CXX_FLAGS}")

if(NOT NOMPI)
  find_package(MPI REQUIRED)
endif()

if(NOT NOCUDA)
  find_package(CUDA REQUIRED)
  #SET(CUDA_SOURCES, src/linalg_cuda.cu)
  CUDA_INCLUDE_DIRECTORIES(${INCLUDE_DIRS})
  # see https://stackoverflow.com/questions/36021233/cmake-cuda-compile-cpp-files-in-cuda-mode-x-cu
  set_source_files_properties( src/linalg_accel.cc PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ )
  CUDA_ADD_LIBRARY(genomics_metrics STATIC ${LIB_SOURCES} ${CUDA_SOURCES})
else()
  ADD_LIBRARY(genomics_metrics STATIC ${LIB_SOURCES} ${CUDA_SOURCES})
endif()

INCLUDE_DIRECTORIES(${INCLUDE_DIRS})

ADD_EXECUTABLE(genomics_metric driver/genomics_metric.cc driver/driver.cc driver/test_problems.cc driver/input_output.cc)
TARGET_LINK_LIBRARIES(genomics_metric genomics_metrics)
TARGET_LINK_LIBRARIES(genomics_metric ${CMAKE_EXE_LINKER_FLAGS})
if(TESTING) 
  TARGET_LINK_LIBRARIES(genomics_metric "-L${GTEST_DIR}/lib -lgtest")
endif(TESTING) 

SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)

install(TARGETS genomics_metric DESTINATION bin)
install(TARGETS genomics_metrics DESTINATION lib)
install(FILES ${LIB_INCLUDES} DESTINATION include)

MESSAGE("Compiler ID: ${CMAKE_C_COMPILER_ID}")

#--------------------

if(TESTING) 
  enable_testing()

  foreach(test_target driver_test metrics_test)
    ADD_EXECUTABLE(${test_target} testing/${test_target}.cc driver/driver.cc driver/test_problems.cc driver/input_output.cc)
    TARGET_LINK_LIBRARIES(${test_target} genomics_metrics)
    TARGET_LINK_LIBRARIES(${test_target} ${CMAKE_EXE_LINKER_FLAGS})
    TARGET_LINK_LIBRARIES(${test_target} "-L${GTEST_DIR}/lib -lgtest")

    #add_test(${test_target} env CRAY_CUDA_PROXY=1 OMP_NUM_THREADS=16 aprun -n64 ${CMAKE_CURRENT_BINARY_DIR}/${test_target})
    add_test(${test_target} bash -c "cd ${CMAKE_INSTALL_PREFIX} ; ${TEST_COMMAND} ${CMAKE_CURRENT_BINARY_DIR}/${test_target}")
  endforeach(test_target)
endif(TESTING) 

#==============================================================================
