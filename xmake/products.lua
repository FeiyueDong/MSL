local msl_root = path.join(os.scriptdir(), "..")

add_requires("eigen")

target("msl")
    set_kind("headeronly")
    add_packages("eigen", {public = true})
    add_includedirs(path.join(msl_root, "include"),
                    path.join(msl_root, "include/msl"),
                    {public = true})
    add_headerfiles(path.join(msl_root, "include/msl/**.hpp"))
