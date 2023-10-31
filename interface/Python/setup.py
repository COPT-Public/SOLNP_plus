from setuptools import setup
import platform

if platform.system() == "Linux":
    solnp_lib = "../../build/libsolnp.so"
elif platform.system() == "Windows":
    # solnp_lib = "../build/Debug/solnp.dll"
    solnp_lib = "../../build/libsolnp.dll"
else:
    # raise Exception("%s Platform not supported yet" % platform.system())
    print(platform.system())
    solnp_lib = "../../build/libsolnp.dylib"

setup(
    name="pysolnp",
    version="1.0",
    author="LLT",
    py_modules=["pysolnp", "SOLNP_CONST"],
    data_files=[solnp_lib],
)
