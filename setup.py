from setuptools import setup

setup(name='randommut',
      version='0.1.1',
      description='A package to randomze somatic mutaitons in the genome',
      url='http://vm-ubuntu-dmp:3000/davidmasp/randommut',
      author='David Mas',
      author_email='david.mas@irbbarcelona.org',
      license='MIT',
      packages=['randommut'],
      install_requires=[
          'biopython','tqdm','numpy','pandas'
      ],
      zip_safe=False,
      entry_points = {
        'console_scripts': [
            'randommut = randommut.__main__:main'
        ]
       })
