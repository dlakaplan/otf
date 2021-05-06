from setuptools import setup, find_packages

setup(
    name="otf",
    version="0.1.0",
    description="Fitting and analysis of on-the-fly (OTF) pulsar mapping",
    author="David Kaplan, Joe Swiggum",
    author_email="kaplan@uwm.edu",
    url="",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "extract_pointing=otf.scripts.extract_pointing:main",
            "fit_otf=otf.scripts.fit_otf:main",
        ],
    },
    python_requires=">=3.7",
    include_package_data=True,
    zip_safe=False,
)
