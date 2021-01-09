import pickle
import os
import csv
import json
import subprocess
import numpy as np
import time
from functools import update_wrapper
import gzip
import sys

proto = pickle.HIGHEST_PROTOCOL

def zdump(obj, f_name):
    """
    save and compress an object with pickle.HIGHEST_PROTOCOL

    Parameters
    ----------
    obj: python object
    f_name: file name used to save compressed file 

    Returns
    -------
    A saved file

    """
    f = gzip.open(f_name,'wb', proto)
    pickle.dump(obj,f)
    f.close()

def zload(f_name):
	"""
    load compressed object with pickle.HIGHEST_PROTOCOL

    Parameters
    ----------
    obj: python object
    f_name: file name used to save compressed file 

    Returns
    -------
    the object in the file

    """
	f = gzip.open(f_name,'rb', proto)
	obj = pickle.load(f)
	f.close()
	return obj

def dump(obj, f_name):
    """
    save an object with pickle.HIGHEST_PROTOCOL

    Parameters
    ----------
    obj: python object
    f_name: file name used to save compressed file 

    Returns
    -------
    A saved file

    """
    f = open(f_name,'wb', proto)
    pickle.dump(obj,f)
    f.close()

def load(f_name):
	"""
    load file with pickle.HIGHEST_PROTOCOL

    Parameters
    ----------
    f_name: file name used to save compressed file 

    Returns
    -------
    the object in the file

    """
	f = open(f_name,'rb', proto)
	obj = pickle.load(f)
	f.close()
	return obj

def mkdir_n(dirName):
	"""
    Function to create directorry in case it does not exist

    Parameters
    ----------
    f_name: file name used to save compressed file 

    Returns
    -------
    the object in the file

    """
	if os.path.isdir(dirName) == False:
		os.mkdir(dirName)

def csv2dict(fname):
	"""
    read a csv file to a dict when keys and values are
    separate by a comma

    Parameters
    ----------
    fname: csv fname with two columns separated by a comma

    Returns
    -------
    a python dictionary

    """
	reader = csv.reader(open(fname, 'r'))
	d = {}
	for row in reader:
		k, v = row
		d[k] = v
	return d

def dict2csv(dict_, fname):
    """
    writes a csv file from a directory with comma separated
    values

    Parameters
    ----------
    dict_: a dictionary we aim to save
    fname: name of the output file

    Returns
    -------
    a csv file

    """
    w = csv.writer(open(fname, "w"))
    for key, val in dict_.items():
        w.writerow([key, val])

def e_vect(n, i):
    """
    get a vector of zeros with a one in the i-th position

    Parameters
    ----------
    n: vector length 
    i: position

    Returns
    -------
    an array with zeros and 1 in the i-th position

    """    
    zeros = [0 for n_i in range(n)]
    zeros[i] = 1
    return zeros

def dict2json(dict_, fname):
    """
    writes a json file from a directory 

    Parameters
    ----------
    dict_: a dictionary we aim to save
    fname: name of the output file

    Returns
    -------
    a csv file

    """
    with open(fname, 'w') as fp:
        json.dump(dict_, fp)


def json2dict(fname):
    """
    reads a json dictionary from a file

    Parameters
    ----------
    dict_: a dictionary we aim to save
    fname: name of the output file

    Returns
    -------
    a csv file

    """
    with open(fname, 'r') as f:
        dict = json.load(f)
    return dict


def shell(command, printOut=True):
    """
    Run shell commands in Linux, decide if printing or not the output in console

    Parameters
    ----------
    command: text command
    printOut: decide if print output or not

    Returns
    -------
    None

    """   
    if printOut == False:   
        proc = subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        proc.wait(timeout=None)
    else:
        proc = subprocess.Popen(command, shell=True)
        proc.wait(timeout=None)

def write_file(string, fname):
    """
    writes a text file from a string 

    Parameters
    ----------
    string: some text in string format
    fname: name of the output file

    Returns
    -------
    a text file

    """
    f = open(fname, 'w')
    f.write(string)
    f.close()

def read_vector(fname, delimiter):
    """
    read a numpy vector

    Parameters
    ----------
    fname: name of the file where the vector is stored
    delimiter

    Returns
    -------
    a numpy array

    """
    x = np.loadtxt(fname, delimiter=delimiter)

def decorator(d):
    """
    Make function d a decorator: d wraps a function fn.

    Parameters
    ----------
    d: function 
    
    """
    def _d(fn):
        return update_wrapper(d(fn), fn)
    update_wrapper(_d, d)
    return _d

def list2txt(list_, fname):
    """
    writes a text file from a list

    Parameters
    ----------
    list_: a python list
    fname: name of the output file

    Returns
    -------
    a text file

    """
    with open(fname, 'w') as filehandle:
        for listitem in list_:
            filehandle.write('%s\n' % listitem)

def txt2list(fname):
    """
    writes a text file from a list

    Parameters
    ----------
    list_: a python list
    fname: name of the output file

    Returns
    -------
    a text file

    """
    crimefile = open(fname, 'r')
    result = [line.split(',') for line in crimefile.readlines()]
    return result

@decorator
def timeit(f):
    """
    time a function by applying a decorator

    Parameters
    ----------
    f: function 

    Returns
    -------
    a decorated function 

    """
    def new_f(*args, **kwargs):
        bt = time.time()
        r = f(*args, **kwargs)
        et = time.time()
        print("time spent on {0}: {1:.2f}s".format(f.__name__, et - bt))
        return r
    return new_f

def get_integers(str_):
    return ''.join(filter(lambda i: i.isdigit(), str_))

class progressBar:
    def __init__(self, toolbar_width):
        self.toolbar_width = toolbar_width

    def set(self):
        # setup toolbar
        print("\nProgress:")
        sys.stdout.write("[%s]" % (" " * self.toolbar_width))
        sys.stdout.flush()
        sys.stdout.write("\b" * (self.toolbar_width+1)) # return to start of line, after '['

    def tic(self):
        sys.stdout.write("-")
        sys.stdout.flush()

    def finish(slef):
        sys.stdout.write("]\n")