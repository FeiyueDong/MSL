target("test_signal")
    set_kind("binary")
    add_files("./*.cpp")
    add_deps("msl")