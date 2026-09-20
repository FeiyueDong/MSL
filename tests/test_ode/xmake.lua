target("test_ode")
    set_kind("binary")
    add_files("./*.cpp")
    add_deps("msl")
