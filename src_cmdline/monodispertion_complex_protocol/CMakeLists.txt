#
# Protocol executable.
#
add_executable(monodispertion_complex_protocol monodispertion_complex_protocol.cpp)
target_link_libraries(monodispertion_complex_protocol ${Boost_LIBRARIES})
install(TARGETS monodispertion_complex_protocol DESTINATION bin)
