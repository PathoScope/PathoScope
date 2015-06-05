from setuptools import setup, find_packages

setup(
    name="PathoScope",
    version="2.0.6",
    packages=find_packages(),
    install_requires=[],
    zip_safe=False,
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "pathoscope = pathoscope.__main__:main"
        ]
    },
    namespace_packages=["pathoscope"]
    )
