set_project("MSL")
set_version("1.0.0")
set_xmakever("3.0.3")
set_warnings("all")
set_allowedplats("windows", "linux", "macosx", "mingw")

add_rules("mode.debug", "mode.release")

set_languages("c++20")

if is_plat("windows") then
    set_toolchains("msvc")
elseif is_plat("mingw") then
    local msys2_root = os.getenv("MSYS2_ROOT")
    if msys2_root and #msys2_root > 0 then
        set_config("sdk", msys2_root)
    end
    set_toolchains("gcc")
end

add_requires("eigen")

if is_plat("linux", "macosx", "mingw") then
    add_cxflags("-fPIC")
end

target("msl")
    set_kind("headeronly")
    add_packages("eigen", {public = true})
    add_includedirs("include", {public = true})
    add_headerfiles("include/msl/**.hpp")

includes("tests")
