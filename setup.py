from setuptools import setup, find_packages


setup(name='strkernel',
      version='0.1',
      description='Collection of string kernels',
      long_description='This packages includes a variety of string kernel methods. Each method transform the list of input strings into a format that can be used by machine learning algorithm (e.g. SVM).',
      classifiers=[
          'Development Status :: 3 - Alpha',
          'License :: MIT License',
          'Programming Language :: Python :: 3.6',
          'Topic :: Text Processing :: Machine Learning',
      ],
      keywords='string kernel SVM machine learning',
      url='https://github.com/jakobhertzberg/string-kernel',
      author='Meng Zhang, Mitra Darvish, Jakob Hertzberg',
      author_email='jakob.hertzberg@gmail.com, mdarvish@posteo.de',
      license='MIT',
      test_suite='nose.collector',
      tests_require=['nose'],
      packages=find_packages(),
      install_requires=[
      'numpy',
      'scipy',
        'Biopython',
        'nbsphinx'
      ],
      include_package_data=True,
      zip_safe=False)
