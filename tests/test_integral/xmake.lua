target("test_integral")
    set_kind("binary")
    add_files("./*.cpp")
    add_deps("msl")