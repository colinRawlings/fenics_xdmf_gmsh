# !python3

from setuptools import setup

if __name__ == "__main__":
    setup(name="fenics-utils",
          version="0.1.0",
          packages=[
              "fenics_utils",
          ],
          package_data={"": [
              "assets/*.geo",
          ]})
