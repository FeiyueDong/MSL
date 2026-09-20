target("test_interp")
    set_kind("binary")
    add_files("./*.cpp")
    add_deps("msl")