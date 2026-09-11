target("test_equation")
    set_kind("binary")
    add_files("./*.cpp")
    add_deps("msl")
