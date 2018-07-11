# pyREES
Python framework for Robust Easy Environment for Simulation (pyREES)


## Installing OpenMesh Dependency for pyREES

One need to setup OpenMesh for Python in order for PyREES to work.

First one need to ensure that boost is correctly installed on the system. We
need to make sure we have python 3.5 for Boost. We are using MacPort for this.
One can test if one has the correct boost setup by writing

```bash
port installed boost
```

If one does not have "boost @1.59.0_2+no_single+python35+universal" then one
can install this by writing:

```bash
sudo port install boost @1.66.0_3+no_single+python36+universal
```

Next one must download OpenMesh source code from

https://www.openmesh.org/media/Releases/6.3/OpenMesh-6.3.zip

And uzip this in any local folder, it does not really matter where, for instance

```bash
/Users/name/Documents/
```

Now we have to build and install OpenMesh

```bash
cd /Users/name/Documents/OpenMesh-6.3
mkdir build
cmake .. -DOPENMESH_PYTHON_VERSION=3
make
sudo make install
sudo make doc
```

If successfully installed then one may use the OpenMesh python package by writing

```python
import sys

sys.path.append('/usr/local/lib/python/')
```

One may test the installation by using the Python console, like this:

```python
import sys
sys.path.append('/usr/local/lib/python/')
from openmesh import *
mesh = TriMesh()
```
