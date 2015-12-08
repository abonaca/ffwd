#include <Python.h>
#include <numpy/arrayobject.h>
#include "ffwd.h"

/* Docstrings */
static char module_docstring[] = "This module provides an interface for speedups in the fast-forward potential recovery method.";
static char compare_docstring[] = "Compare two distributions of points in arbitrary number of dimensions.";
    
/* Available functions */
static PyObject *ffwd_compare(PyObject *self, PyObject *args);

/* Module specification */
static PyMethodDef module_methods[] = {
	{"compare", ffwd_compare, METH_VARARGS, compare_docstring},
	{NULL, NULL, 0, NULL}
};

/* Initialize the module */
PyMODINIT_FUNC initffwd(void)
{
	PyObject *m = Py_InitModule3("ffwd", module_methods, module_docstring);
	if (m == NULL)
		return;

	/* Load `numpy` functionality. */
	import_array();
}

double *pyvector_to_Carrayptrs(PyArrayObject *arrayin);

// static PyObject* func(PyObject* self, PyObject* args) {
//   PyObject *list2_obj;
//   PyObject *list3_obj;
//   if (!PyArg_ParseTuple(args, "OO", &list2_obj, &list3_obj))
//     return NULL;
// 
//   double **list2;
//   double ***list3;
// 
//   //Create C arrays from numpy objects:
//   int typenum = NPY_DOUBLE;
//   
//   PyArray_Descr *descr;
//   
//   descr = PyArray_DescrFromType(typenum);
//   
//   npy_intp dims[3];
//   
//   if (PyArray_AsCArray(&list2_obj, (void **)&list2, dims, 2, descr) < 0 || PyArray_AsCArray(&list3_obj, (void ***)&list3, dims, 3, descr) < 0) {
//     PyErr_SetString(PyExc_TypeError, "error converting to c array");
//     return NULL;
//   }
//   
//   printf("2D: %f, 3D: %f.\n", list2[3][1], list3[1][0][2]);
// }

static PyObject *ffwd_compare(PyObject *self, PyObject *args)
{
	int i, j, l, N, k, K;
	double *logdet, logp;
	
	PyObject *x_obs_obj, *x_mod_obj, *err_obs_obj, *sigma_inv_obj, *logdet_obj;
	PyObject *x_obs_array, *x_mod_array, *err_obs_array, *sigma_inv_array, *logdet_array;

	//printf("here\t");

	// Parse the input tuple
	if (!PyArg_ParseTuple(args, "OOOOOiii", &x_obs_obj, &x_mod_obj, &err_obs_obj, &sigma_inv_obj, &logdet_obj, &N, &k, &K))	// reads in input parameters
		return NULL;
	
// 	printf("here\n");
	// Interpret the input parameters as numpy arrays
	x_obs_array = PyArray_FROM_OTF(x_obs_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
	x_mod_array = PyArray_FROM_OTF(x_mod_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
	err_obs_array = PyArray_FROM_OTF(err_obs_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
	sigma_inv_array = PyArray_FROM_OTF(sigma_inv_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
	logdet_array = PyArray_FROM_OTF(logdet_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
	
 	//printf("here\t");
	//If that didn't work, throw an exception
	if (x_obs_array == NULL) {
		Py_XDECREF(x_obs_array);
		return NULL;
	}
	if (x_mod_array == NULL) {
		Py_XDECREF(x_mod_array);
		return NULL;
	}
	if (err_obs_array == NULL) {
		Py_XDECREF(err_obs_array);
		return NULL;
	}
	if (sigma_inv_array == NULL) {
		Py_XDECREF(sigma_inv_array);
		return NULL;
	}
	if (logdet_array == NULL) {
		Py_XDECREF(logdet_array);
		return NULL;
	}
	
	
	// Read data from python objects (or rather get pointers to the first element of the array)
	double *x_obs_ptr, *x_mod_ptr, *err_obs_ptr, *sigma_inv_ptr;
	
	x_obs_ptr = (double*)PyArray_DATA(x_obs_array);
	x_mod_ptr = (double*)PyArray_DATA(x_mod_array);
	err_obs_ptr = (double*)PyArray_DATA(err_obs_array);
	sigma_inv_ptr = (double*)PyArray_DATA(sigma_inv_array);
	logdet = (double*)PyArray_DATA(logdet_array);
	
	// Map 1d arrays (pointer) to 2d and 3d arrays (pointers of pointers)
	double **x_obs = (double**) malloc(k*sizeof(double *));
	double **x_mod = (double**) malloc(k*sizeof(double *));
	double **err_obs = (double**) malloc(k*sizeof(double *));
	double ***sigma_inv = (double***) malloc(N*sizeof(double **));
	
// 	for(i=0;i<N;i++){
// 		sigma_inv[i] = (double**) malloc(k*sizeof(double*));
// 		for(j=0;j<k;j++){
// 			sigma_inv[i][j] = (double*) malloc(k*sizeof(double));
// 		}
// 	}
	
	for(i=0;i<k;i++){
		x_obs[i] = (double*) malloc(N*sizeof(double));
		x_mod[i] = (double*) malloc(K*sizeof(double));
		err_obs[i] = (double*) malloc(N*sizeof(double));
		for(j=0;j<K;j++){
			x_mod[i][j] = x_mod_ptr[i*K+j];
			if(j<N){
				x_obs[i][j] = x_obs_ptr[i*N+j];
				err_obs[i][j] = err_obs_ptr[i*N+j];
				//printf("%d %d %lf %lf %lf\n", i, j, x_mod[i][j],  x_obs[i][j],  err_obs[i][j]);
			}
			//printf("%d %d %lf\n", i, j, x_mod[i][j]);
		}
	}
	
	//printf("here\t");
	
	for(i=0;i<N;i++){
		sigma_inv[i] = (double**) malloc(k*sizeof(double*));
		for(j=0;j<k;j++){
			sigma_inv[i][j] = (double*) malloc(k*sizeof(double));
			for(l=0;l<k;l++){
				sigma_inv[i][j][l] = sigma_inv_ptr[k*k*i+j*k+l];
				//printf("%d %d %d %lf\n", i, j, l, sigma_inv[i][j][l]);
			}
		}
	}

// 	// Call the external C function to compare two distributions
	logp = compare(x_obs, x_mod, err_obs, sigma_inv, logdet, N, k, K);
	//logp = 0;

 	//printf("%lf %d %d %d \n", logp, N, k, K);
// 
// 	// Check if error raised
// 	if(err!=0) {
// 		PyErr_SetString(PyExc_RuntimeError, "Error occured in the leapfrog integrator.");
// 		return NULL;
// 	}
// 	
	//printf("here1\t");
	// Store return array
	PyObject *out = Py_BuildValue("d", logp);
	
	// Clean up
	for(i=0; i<k; i++){
		free(x_mod[i]);
		free(x_obs[i]);
		free(err_obs[i]);
	}
	
	//printf("here2\t");
	
	for(i=0; i<N; i++){
		for(j=0; j<k; j++)
			free(sigma_inv[i][j]);
		free(sigma_inv[i]);
	}
	
	//printf("here3\t");
	
	//free(x_obs_ptr);
	//free(x_mod_ptr);
	//free(err_obs_ptr);
	//free(sigma_inv_ptr);
	
	//printf("here\t");
	
	free(x_obs);
	
	free(x_mod);
	
	free(err_obs);
	
	free(sigma_inv);
	
	//free(logdet);
	
	//printf("here\t");
	
	
	
	//Py_XDECREF(x_obs_array);
	//Py_XDECREF(x_mod_array);
	//Py_XDECREF(err_obs_array);
	//Py_XDECREF(sigma_inv_array);
	//Py_XDECREF(logdet_array);
	
	//printf("here\t");
	
	// Return likelihood
	return out;
}

double *pyvector_to_Carrayptrs(PyArrayObject *arrayin)
{
// 	int n=arrayin->dimensions[0];
	return (double *) arrayin->data;  /* pointer to arrayin data as double */
}
