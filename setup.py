from setuptools import setup, find_packages
from pathoscope._version import VERSION

setup(
    name="PathoScope",
    version=VERSION,
    packages=find_packages(),
    install_requires=[],
    zip_safe=False,
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "pathoscope = pathoscope.__main__:main"
        ]
    },
    namespace_packages=["pathoscope"],
    scripts=['scripts/pathoscope'],
    )
