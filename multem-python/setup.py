#
# Copyright (C) 2019 James Parkhurst
#
# This code is distributed under the GPLv3 license.
#
from skbuild import setup
from setuptools import Command
from setuptools import find_packages


class GcovCommand(Command):

    description = "run gcovr to output C/C++ coverage statistics"

    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        import os.path
        from gcovr.__main__ import main

        main(
            [
                "--root",
                os.path.dirname(os.path.abspath(__file__)),
                "--exclude",
                ".*CMakeCCompilerId.c$",
            ]
        )


def main():
    """
    Setup the package

    """
    tests_require = ["pytest", "pytest-cov", "mock"]

    setup(
        package_dir={"": "src"},
        packages=find_packages(where="src"),
        install_requires=["numpy"],
        setup_requires=["pytest-runner"],
        tests_require=tests_require,
        test_suite="test",
        cmdclass={"gcov": GcovCommand},
        extras_require={
            "build_sphinx": ["sphinx", "sphinx_rtd_theme"],
            "test": tests_require,
            "gcov": ["gcovr"],
        },
    )


if __name__ == "__main__":
    main()
