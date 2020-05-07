import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="continuous-peak-fit",
    version="0.0.1.dev1",
    author="Simon Hunt & Danielle Fenech",
    author_email="simon.hunt@manchester.ac.uk",
    description="Fits azimuth/time dependency of peaks with Fourier Series descriptions",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ExperimentalMineralPhysics/Continuous-Peak-Fit",
    packages=setuptools.find_packages(),
    scripts=['bin/CPF_XRD_FitPattern', 'bin/CPF_generate_inputs'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=["matplotlib",
                      "pycairo",
                      "lmfit",
                      "h5py",
                      "pillow",
                      "pyFAI"
                      ]
)