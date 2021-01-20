#
# Copyright (C) 2019 James Parkhurst
#
# This code is distributed under the GPLv3 license.
#
from skbuild import setup


def main():
    """
    Setup the package

    """
    tests_require = ["pytest", "pytest-cov", "mock"]

    setup(
        package_dir={"multem": "src"},
        packages=["multem"],
        install_requires=["numpy", "scipy"],
        setup_requires=["pytest-runner"],
        tests_require=tests_require,
        test_suite="test",
        extras_require={
            "build_sphinx": ["sphinx", "sphinx_rtd_theme"],
            "test": tests_require,
        },
    )


if __name__ == "__main__":
    main()
