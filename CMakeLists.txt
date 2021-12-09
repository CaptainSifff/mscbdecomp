cmake_minimum_required (VERSION 2.8.12)
project (mscbdecomp Fortran C)

find_package(lapack)
SET(BLAS_LIBRARIES ${BLAS_LIBRARIES} CACHE INTERNAL "BLAS LIBS" FORCE)
SET(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} CACHE INTERNAL "LAPACK LIBS" FORCE)

add_library( libmscbdecomp Node.F90 Exponentials.F90 MvG_decomp.F90 vertex.F90 colorvertex.F90 GeneralSingleColExp.F90 HomogeneousSingleColExp.F90 SingleColExpBase.F90 TraceLessSingleColExp.F90 ZeroDiagSingleColExp.F90 graphdata.F90 vertex.F90)

add_executable ( mscbdecomp mscbdecomp.F90)
target_link_libraries(mscbdecomp libmscbdecomp)
target_link_libraries(mscbdecomp ${LAPACK_LIBRARIES})
target_link_libraries(mscbdecomp "-llapack -lblas")

add_subdirectory(testsuite)
enable_testing()
