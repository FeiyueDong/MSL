target("test_polynomial")
    set_kind("binary")
    add_files("./*.cpp")
    add_deps("msl")