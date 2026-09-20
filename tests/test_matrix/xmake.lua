target("test_matrix")
    set_kind("binary")
    add_files("./*.cpp")
    add_deps("msl")