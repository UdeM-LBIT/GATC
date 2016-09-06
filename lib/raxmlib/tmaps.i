// Typemaps that are independent of scripting language

// cleans up the char ** array after the function call
%typemap(freearg) char ** {
    free($1);
}

// cleans up the char ** array after the function call
%typemap(freearg) (int argc, char **argv) {
    free((char *) $2);
}


// Python typemaps
#ifdef SWIGPYTHON

// treat char ** as a special case
%typemap(in) char ** {
    /* Check if is a list */
    if (PyList_Check($input)) {
        int size = PyList_Size($input);
        int i = 0;
        $1 = (char **) malloc((size+1)*sizeof(char *));
        for (i = 0; i < size; i++) {
            PyObject *o = PyList_GetItem($input,i);
            if (PyString_Check(o))
                $1[i] = PyString_AsString(PyList_GetItem($input,i));
            else {
                PyErr_SetString(PyExc_TypeError,"list must contain strings");
                free($1);
                return NULL;
            }
        }
        $1[i] = 0;
    }
    else {
        PyErr_SetString(PyExc_TypeError,"not a list");
        return NULL;
    }
}

// treat (int argc, char **argv) as a special case
%typemap(in) (int argc, char **argv) {
    /* Check if is a list */
    if (PyList_Check($input)) {
        int i;
        $1 = PyList_Size($input);
        $2 = (char **) malloc(($1+1)*sizeof(char *));
        for (i = 0; i < $1; i++) {
            PyObject *o = PyList_GetItem($input,i);
            if (PyString_Check(o))
                $2[i] = PyString_AsString(PyList_GetItem($input,i));
            else {
                PyErr_SetString(PyExc_TypeError,"list must contain strings");
                free($2);
                return NULL;
            }
        }
        $2[i] = 0;
    }
    else {
        PyErr_SetString(PyExc_TypeError,"not a list");
        return NULL;
    }
}

// convert between python and C file handle
%typemap(in) FILE * {
    if (!PyFile_Check($input)) {
        PyErr_SetString(PyExc_TypeError, "$1_name must be a file type.");
        return NULL;
    }
    $1 = PyFile_AsFile($input);
}

#endif



// a few test functions
//%inline %{
//int print_args(char **argv) {
//    int i = 0;
//    while (argv[i]) {
//        printf("argv[%d] = %s\n", i,argv[i]);
//        i++;
//    }
//    return i;
//}

// Returns a char ** list
//char **get_args() {
//    static char *values[] = {"Dave","Mike","Susan","John","Michelle",0};
//    return &values[0];
//}
//%}
