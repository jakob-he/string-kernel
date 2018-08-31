from setuptools import setup, find_packages


setup(name='strkernel',
      version='0.2',
      description='Collection of string kernels',
      long_description='This packages includes a variety of string kernel methods. Each method transforms the list of input strings into a format that can be used by machine learning algorithm (e.g. SVM).',
      classifiers=[
          'Development Status :: 3 - Alpha',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3.6',
          'Topic :: Text Processing',
      ],
      keywords='string kernel SVM machine learning',
      url='https://github.com/jakobhertzberg/string-kernel',
      author='Meng Zhang, Mitra Darvish, Jakob Hertzberg',
      author_email='jakob.hertzberg@gmail.com, mdarvish@posteo.de',
      test_suite='nose.collector',
      tests_require=['nose'],
      packages=find_packages(),
      install_requires=[
        'scipy',
        'numpy',
        'Biopython'
      ],
      include_package_data=True,
      zip_safe=False)
