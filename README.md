#  Interpret intermittent bugs in Bugzilla
This project attempts to interpret the exact happening time of intermittent bugs in Bugzilla, it leverages intermittent failure rates and textual comments, and make use of tools including change detector, TF-IDF and Word2Vec to analyse whether or not there is difference of language used before and after intermittent bugs. All the tools used in this project are based on python framework.

## Getting started
There are two python change detectors "SEED" and "ANGLE" in drift_detection converted from Java, which are mainly codes used in this project, and they could be configured and made compatible with Scikit-Multiflow package on your local machines. The other Jupyter files contains the documented experiments I have conducted in the project, which are easily runnable as reference of similar works.

### Prerequisites
In order to run the experiments, a lot of python packages are required to install in your local machine firstly, including Numpy, Pandas, Scipy, Scikit-Multiflow, Sklearn, Matplotlib and Gensim. Most of them are handy installed if you have Jupyter.

### Installing
Here is the steps showing how to make "SEED" and "ANGLE" compatible with Scikit-Multiflow on your local machine. Firstly, make sure Scikit-Multiflow installed successfully,  and put code files (SEED.py, ANGLE.py) into the correct file path \src\skmultiflow\drift_detection. Then the corresponding new class names "SEED", "ANGLE" should be added into constructor file, please make sure with the correct format as their originals. Or you could simply replace drift_detection with mine in this repository. By that, all the setup are finished, please reinstall the package as long as you have modified the codes.

