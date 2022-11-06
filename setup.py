import platform
import sysconfig
import os
import shutil
from setuptools import setup, find_packages

# macOS-12.6-arm64-arm-64bit
# macOS-10.15.7-x86_64-i386-64bit
# Linux-4.18.0-193.28.1.el8_2.x86_64-x86_64-with-glibc2.10
# Windows-10-10.0.19041-SP0
python_minor_version = platform.python_version().split('.')[1]

sysos = platform.platform(aliased=True).lower()
if 'win' in sysos:
    file = 'win64/ispace.cp3'+python_minor_version+'-win_amd64.pyd'
elif 'macos' in sysos and 'arm64' in sysos:    
    file = 'macarm/ispace.cpython-3'+python_minor_version+'-darwin.so'
elif 'macos' in sysos and 'x86_64' in sysos:
    file = 'macintel/ispace.cpython-3'+python_minor_version+'-darwin.so'
elif 'linux' in sysos:
    file = 'linux64/ispace.cpython-3'+python_minor_version+'-x86_64-linux-gnu.so'
else: 
    raise OSError('Unknown OS:'+sysos)

# Check if iSpace module files exists
proj_dir = os.path.dirname(os.path.abspath(__file__))
file_abspath = proj_dir+'/iSpace/libs/'+file
if not os.path.isfile(file_abspath):
    raise FileNotFoundError('No iSpace module found in project folder:'+ file_abspath)

# copy the iSpace module file to the site-package folder.
site_path = sysconfig.get_paths()["purelib"]
shutil.copy(file_abspath, site_path)
print('iSpace module is copied to: ', site_path+os.sep+file)


setup(
    name='iSpace',
    version='1.0.0',
    description='iSpace',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/xnchu/iSpace',
    author='Xiangning Chu',
    author_email='xiangning.chu@colorado.edu',
    license='MIT',
    classifiers=[ 'Development Status :: 5 - Production/Stable', 
                 'Intended Audience :: Science/Research',
                 'Topic :: Scientific/Engineering',
                 'License :: OSI Approved :: MIT License',
                 'Programming Language :: Python :: 3',
                 ],
    keywords='iSpace geopack Tsyganenko magnetic field data tools',
    project_urls={'Information': 'https://github.com/xnchu/iSpace',
                  },
    packages=find_packages(),
    install_requires=['numpy'],
    python_requires='>=3.7',
    include_package_data=True,
)
