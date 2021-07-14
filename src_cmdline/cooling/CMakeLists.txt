#
# Cooling executable.
#
add_executable(cooling cooling.cpp)
target_link_libraries(cooling ${Boost_LIBRARIES})
install(TARGETS cooling DESTINATION bin)
