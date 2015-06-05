import subprocess

has_setuptools = True
try:
    from setuptools import setup
except ImportError:
    has_setuptools = False

def main():
    if not has_setuptools:
        print("Couldn't import setuptools. Bootstrapping a user-site setuptools ")
        # Install setuptools using ez_setup.py bootstrapper
        subprocess.check_call("python ez_setup.py --user", shell=True)
    print("Installing PathoScope")
    subprocess.check_call("python setup.py develop --user", shell=True)
    import pathoscope

if __name__ == "__main__":
    main()
