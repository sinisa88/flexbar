# Add FindTBB directory to CMake's module path
# set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${FLEXBAR_SOURCE_DIR}/FindTBB/CMake/FindTBB/")

#find_package(TBB REQUIRED)
#if(NOT TBB_FOUND)
#       MESSAGE (STATUS "TBB not found. Project will not be built.")
#else(NOT TBB_FOUND)

#LINK_DIRECTORIES(${FLEXBAR_SOURCE_DIR}/lib)
#set(FLEXBAR_64BIT_ARCHITECTURE 1 CACHE INTERNAL "Architecture-bits")


INCLUDE_DIRECTORIES(${FLEXBAR_SOURCE_DIR}/include)

if(CMAKE_SIZEOF_VOID_P MATCHES "8")
	message(STATUS "Architecture: 64 bit")
	if (WIN32)
		LINK_DIRECTORIES(${FLEXBAR_SOURCE_DIR}/lib/win64)
	elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
		LINK_DIRECTORIES(${FLEXBAR_SOURCE_DIR}/lib/macosx)
	else() ## assuming linux, or adjust to your system's tbb lib
		LINK_DIRECTORIES(${FLEXBAR_SOURCE_DIR}/lib/linux64)
	endif()
else()
	message(STATUS "Architecture: 32 bit")
	if (WIN32)
		LINK_DIRECTORIES(${FLEXBAR_SOURCE_DIR}/lib/win32)
	else()
		message(FATAL_ERROR "Platform not supported!")
	endif()
endif()


ADD_EXECUTABLE(flexbar Flexbar.cpp)
TARGET_LINK_LIBRARIES(flexbar tbb)

#ADD_EXECUTABLE(PairedreadFinder PairedreadFinder.cpp)
#TARGET_LINK_LIBRARIES(PairedreadFinder tbb)

