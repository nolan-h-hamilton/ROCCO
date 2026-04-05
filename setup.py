import glob
import os
import subprocess
import sys
import sysconfig
import textwrap

from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
VENDORED_HTSLIB_DIR = os.path.join(ROOT_DIR, "vendor", "htslib")
VERSION_PY_PATH = os.path.join(ROOT_DIR, "rocco", "_version.py")
PREFIX_INCLUDE_DIR = os.path.join(sys.prefix, "include")
PREFIX_LIB_DIR = os.path.join(sys.prefix, "lib")
HTSLIB_CONFIG_MK_PATH = os.path.join(VENDORED_HTSLIB_DIR, "config.mk")
HTSLIB_CONFIG_H_PATH = os.path.join(VENDORED_HTSLIB_DIR, "config.h")
HTSCODECS_CONFIGURE_AC_PATH = os.path.join(
    VENDORED_HTSLIB_DIR,
    "htscodecs",
    "configure.ac",
)
HTSCODECS_VERSION_H_PATH = os.path.join(
    VENDORED_HTSLIB_DIR,
    "htscodecs",
    "htscodecs",
    "version.h",
)
BASE_COMPILE_ARGS = [
    "-O3",
    "-fno-trapping-math",
    "-fno-math-errno",
    "-mtune=generic",
]


class get_numpy_include:
    def __str__(self) -> str:
        import numpy

        return numpy.get_include()


def has_vendored_htslib() -> bool:
    return os.path.exists(os.path.join(VENDORED_HTSLIB_DIR, "Makefile"))


def get_rocco_version() -> str:
    with open(VERSION_PY_PATH, "r", encoding="utf-8") as handle:
        for line in handle:
            line_ = line.strip()
            if line_.startswith("__version__"):
                _, value = line_.split("=", 1)
                return value.strip().strip('"').strip("'")
    raise RuntimeError("Could not determine ROCCO version.")


def get_htslib_include_dirs() -> list[str]:
    include_dirs = []
    if has_vendored_htslib():
        include_dirs.extend(
            [
                VENDORED_HTSLIB_DIR,
                os.path.join(VENDORED_HTSLIB_DIR, "htslib"),
            ]
        )
    include_dirs.extend(
        [
            PREFIX_INCLUDE_DIR,
            os.path.join(PREFIX_INCLUDE_DIR, "htslib"),
        ]
    )
    return include_dirs


def get_library_dirs() -> list[str]:
    return [PREFIX_LIB_DIR]


def get_bundled_htslib_archive() -> str:
    return os.path.join(VENDORED_HTSLIB_DIR, "libhts.a")


def find_static_library(library_name: str) -> str | None:
    candidate_dirs = [
        PREFIX_LIB_DIR,
        "/usr/local/lib",
        "/opt/homebrew/lib",
        "/usr/lib",
        "/usr/lib64",
        "/lib",
        "/lib64",
    ]
    library_pattern = f"lib{library_name}.a"
    for candidate_dir in candidate_dirs:
        candidate_path = os.path.join(candidate_dir, library_pattern)
        if os.path.exists(candidate_path):
            return candidate_path
    glob_patterns = [
        f"/usr/lib/*/{library_pattern}",
        f"/lib/*/{library_pattern}",
        f"/usr/local/lib/*/{library_pattern}",
    ]
    for pattern in glob_patterns:
        matches = sorted(glob.glob(pattern))
        if matches:
            return matches[0]
    return None


def get_bundled_dependency_archives() -> list[str]:
    dependency_archives = []
    if sys.platform == "darwin":
        static_zlib = find_static_library("z")
        if static_zlib is not None:
            dependency_archives.append(static_zlib)
    return dependency_archives


def get_bundled_htslib_libraries() -> list[str]:
    libraries = []
    if sys.platform != "darwin" or find_static_library("z") is None:
        libraries.append("z")
    if sys.platform.startswith("linux"):
        libraries.extend(["m", "pthread"])
    return libraries


def get_bundled_htslib_extra_objects() -> list[str]:
    return [get_bundled_htslib_archive()] + get_bundled_dependency_archives()


def write_text_if_changed(path: str, contents: str) -> None:
    existing_contents = None
    if os.path.exists(path):
        with open(path, "r", encoding="utf-8") as handle:
            existing_contents = handle.read()
    if existing_contents == contents:
        return
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(contents)


def join_compiler_flags(*flag_groups: list[str]) -> str:
    return " ".join(
        flag
        for flag_group in flag_groups
        for flag in flag_group
        if flag
    )


def get_vendored_htslib_cppflags() -> str:
    return join_compiler_flags(
        [f"-I{PREFIX_INCLUDE_DIR}"],
        [os.environ.get("CPPFLAGS", "").strip()],
    )


def get_vendored_htslib_cflags() -> str:
    base_cflags = [
        "-g",
        "-Wall",
        "-O3",
        "-fvisibility=hidden",
        "-fPIC",
    ]
    if sys.platform == "darwin":
        deployment_target = os.environ.get("MACOSX_DEPLOYMENT_TARGET")
        if not deployment_target:
            deployment_target = sysconfig.get_config_var(
                "MACOSX_DEPLOYMENT_TARGET"
            )
        if deployment_target:
            base_cflags.append(f"-mmacosx-version-min={deployment_target}")
    extra_cflags = os.environ.get("CFLAGS", "").strip()
    return join_compiler_flags([extra_cflags], base_cflags)


def get_vendored_htslib_ldflags() -> str:
    return join_compiler_flags([os.environ.get("LDFLAGS", "").strip()])


def get_vendored_htslib_config_mk() -> str:
    return textwrap.dedent(
        f"""\
        CC = {os.environ.get("CC", "cc")}
        RANLIB = {os.environ.get("RANLIB", "ranlib")}

        CPPFLAGS = {get_vendored_htslib_cppflags()}
        CFLAGS = {get_vendored_htslib_cflags()}
        LDFLAGS = {get_vendored_htslib_ldflags()}
        VERSION_SCRIPT_LDFLAGS =
        LIBS = -lz -lm

        NONCONFIGURE_OBJS =
        plugin_OBJS =
        noplugin_LDFLAGS =
        noplugin_LIBS =

        REF_CACHE_PROGRAMS =
        HTS_CFLAGS_AVX2 =
        HTS_CFLAGS_AVX512 =
        HTS_CFLAGS_SSE4 =
        """
    )


def get_vendored_htslib_config_h() -> str:
    return textwrap.dedent(
        """\
        /* rocco vendored htslib config */
        #ifndef _XOPEN_SOURCE
        #define _XOPEN_SOURCE 600
        #endif
        #define HAVE_DRAND48 1
        #if defined __x86_64__
        #define HAVE_X86INTRIN_H 1
        #endif
        #if defined __x86_64__ || defined __arm__ || defined __aarch64__
        #define HAVE_ATTRIBUTE_CONSTRUCTOR 1
        #endif
        #if defined __linux__
        #define HAVE_GETAUXVAL 1
        #elif defined __FreeBSD__
        #define HAVE_ELF_AUX_INFO 1
        #elif defined __OpenBSD__
        #define HAVE_OPENBSD 1
        #endif
        """
    )


def get_vendored_htscodecs_version() -> str:
    if not os.path.exists(HTSCODECS_CONFIGURE_AC_PATH):
        raise FileNotFoundError("Vendored htscodecs configure.ac is missing.")

    with open(HTSCODECS_CONFIGURE_AC_PATH, "r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped.startswith("AC_INIT("):
                continue
            fields = [field.strip() for field in stripped.split(",")]
            if len(fields) < 2:
                break
            version_field = fields[1].rstrip(")")
            version_field = version_field.strip().strip("[]")
            version_field = version_field.strip().strip('"').strip("'")
            if version_field:
                return version_field
            break
    raise RuntimeError("Could not determine vendored htscodecs version.")


def get_vendored_htscodecs_version_h() -> str:
    return textwrap.dedent(
        f"""\
        #define HTSCODECS_VERSION_TEXT "{get_vendored_htscodecs_version()}"
        """
    )


def prepare_vendored_htslib_build() -> None:
    write_text_if_changed(
        HTSLIB_CONFIG_MK_PATH,
        get_vendored_htslib_config_mk(),
    )
    write_text_if_changed(
        HTSLIB_CONFIG_H_PATH,
        get_vendored_htslib_config_h(),
    )
    write_text_if_changed(
        HTSCODECS_VERSION_H_PATH,
        get_vendored_htscodecs_version_h(),
    )


def build_vendored_htslib() -> None:
    if not has_vendored_htslib():
        raise FileNotFoundError("Vendored HTSlib source tree is missing.")
    prepare_vendored_htslib_build()
    subprocess.check_call(
        ["make", "-C", VENDORED_HTSLIB_DIR, "clean"],
        cwd=ROOT_DIR,
    )
    prepare_vendored_htslib_build()
    subprocess.check_call(
        ["make", "-C", VENDORED_HTSLIB_DIR, "lib-static"],
        cwd=ROOT_DIR,
    )
    if not os.path.exists(get_bundled_htslib_archive()):
        raise FileNotFoundError("Failed to build vendored libhts.a.")


class build_rocco_ext(build_ext):
    def run(self):
        if has_vendored_htslib():
            build_vendored_htslib()
        super().run()


with open("README.md", "r", encoding="utf-8") as readme_file:
    long_description = readme_file.read()


extensions = [
    Extension(
        "rocco._baseline",
        sources=[
            "rocco/_baseline.c",
            "rocco/native/baseline_backend.c",
        ],
        include_dirs=[get_numpy_include(), "rocco"],
        extra_compile_args=BASE_COMPILE_ARGS,
    ),
    Extension(
        "rocco._chain_dp",
        sources=["rocco/_chain_dp.c"],
        include_dirs=[get_numpy_include()],
        extra_compile_args=BASE_COMPILE_ARGS,
    ),
    Extension(
        "rocco._wls",
        sources=[
            "rocco/_wls.c",
            "rocco/native/wls_backend.c",
        ],
        include_dirs=[get_numpy_include(), "rocco"],
        extra_compile_args=BASE_COMPILE_ARGS,
    ),
    Extension(
        "rocco._hts_counts",
        sources=[
            "rocco/_hts_counts.c",
            "rocco/native/ccounts_backend.c",
        ],
        include_dirs=[get_numpy_include(), "rocco"] + get_htslib_include_dirs(),
        libraries=get_bundled_htslib_libraries(),
        library_dirs=get_library_dirs(),
        extra_objects=get_bundled_htslib_extra_objects(),
        extra_compile_args=BASE_COMPILE_ARGS,
    ),
]


setup(
    name="rocco",
    version=get_rocco_version(),
    author="Nolan Holt Hamilton",
    author_email="nolan.hamilton@unc.edu",
    description="Multisample Consensus Peak Calling via Convex Optimization",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/nolan-h-hamilton/rocco",
    packages=find_packages(),
    license="MIT",
    license_files=[],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Programming Language :: Python :: 3.14",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    keywords=[
        "genomics",
        "functional genomics",
        "epigenomics",
        "epigenetics",
        "peak calling",
        "chromatin accessibility",
        "chromatin",
        "consensus peak",
        "DNase",
        "ATAC",
        "ChIP-seq",
    ],
    python_requires=">=3.10, <4",
    setup_requires=["numpy"],
    ext_modules=extensions,
    cmdclass={"build_ext": build_rocco_ext},
    install_requires=[
        "numpy",
        "scipy",
        "pandas",
        "pysam",
        "pyBigWig",
    ],
    extras_require={"pytest": ["pytest>=6.0.1"]},
    entry_points={"console_scripts": ["rocco = rocco.rocco:main"]},
    include_package_data=True,
    zip_safe=False,
)
