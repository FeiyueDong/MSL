set_project("MSL")
set_version("1.0.0")
set_xmakever("3.0.3")
set_warnings("all")
set_allowedplats("windows", "linux", "macosx", "mingw")

add_rules("mode.debug", "mode.release")

if is_plat("windows") then
    set_toolchains("msvc")
elseif is_plat("mingw") then
    -- MSYS2/GCC installations on PATH are auto-detected. Set the environment
    -- variable MSYS2_ROOT only when the SDK lives outside PATH.
    local msys2_root = os.getenv("MSYS2_ROOT")
    if msys2_root and #msys2_root > 0 then
        set_config("sdk", msys2_root)
    end
    set_toolchains("gcc")
end

set_languages("c++20")

if is_plat("mingw") then
    set_targetdir("$(projectdir)/build/mingw",{ bindir = "bin", libdir = "lib" })
elseif is_plat("windows") then
    set_targetdir("$(projectdir)/build/windows",{ bindir = "bin", libdir = "lib" })
elseif is_plat("linux") then
    set_targetdir("$(projectdir)/build/linux",{ bindir = "bin", libdir = "lib" })
elseif is_plat("macosx") then
    set_targetdir("$(projectdir)/build/macosx",{ bindir = "bin", libdir = "lib" })
end

if is_plat("linux", "macosx", "mingw") then
    add_cxflags("-fPIC")
end

-- Header-only library target. Tests depend on "msl" and inherit its public
-- include directories.
--
-- Eigen is the single external dependency. MSYS2/GCC, Linux and macOS system
-- installations are found through the compiler default include paths; an
-- alternative location can be supplied with the EIGEN_INCLUDE_DIR environment
-- variable.
target("msl")
    set_kind("headeronly")
    add_includedirs("src/msl", {public = true})
    add_headerfiles("src/msl/**.hpp")

    local eigen_include = os.getenv("EIGEN_INCLUDE_DIR")
    if eigen_include and #eigen_include > 0 then
        add_includedirs(eigen_include, {public = true})
    end
    if is_plat("windows") then
        add_includedirs("vcpkg_installed/x64-windows/x64-windows/include",
            {public = true})
    elseif is_plat("linux", "macosx") then
        add_includedirs("/usr/include/eigen3", "/usr/local/include/eigen3",
            "/opt/homebrew/include/eigen3", {public = true})
    end

includes("src/tests")
