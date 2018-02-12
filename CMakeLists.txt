cmake_minimum_required (VERSION 2.8.5)
project (mscbdecomp Fortran C)

find_package(lapack)
SET(BLAS_LIBRARIES ${BLAS_LIBRARIES} CACHE INTERNAL "BLAS LIBS" FORCE)
SET(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} CACHE INTERNAL "LAPACK LIBS" FORCE)

add_executable ( mscbdecomp mscbdecomp.F90 )
target_link_libraries(mscbdecomp ${LAPACK_LIBRARIES})


