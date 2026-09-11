target("test_difference")
    set_kind("binary")
    add_files("./*.cpp")
    add_deps("msl")