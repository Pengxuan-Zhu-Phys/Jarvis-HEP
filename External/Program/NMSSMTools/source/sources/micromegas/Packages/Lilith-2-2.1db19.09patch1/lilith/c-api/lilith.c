/*
##########################################################################
#
#  This file is part of Lilith
#  v1 (2015) by Jeremy Bernon and Beranger Dumont 
#  v2 (2019) by Sabine Kraml, Tran Quang Loc, Dao Thi Nhung, Le Duc Ninh 
#
#  Web page: http://lpsc.in2p3.fr/projects-th/lilith/
#
#  In case of questions email sabine.kraml@lpsc.in2p3.fr 
#
#
#    Lilith is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Lilith is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Lilith.  If not, see <http://www.gnu.org/licenses/>.
#
##########################################################################
*/

#include <Python.h>
#include "lilith.h"

int userread = 0;

/**
Import Lilith, create a Lilith object and read experimental results
**/
PyObject* initialize_lilith(char* experimental_input)
{
      PyObject* args;

      // |--------------------------------|
      // | Adding Lilith path to sys.path |
      // |--------------------------------|

      char pathtolilith[200];
      // remove "lilith.c" from the path (8 characters + '\0')
      strncpy(pathtolilith, __FILE__, sizeof(__FILE__)-9);
      pathtolilith[sizeof(__FILE__)-9] = '\0';
      // add "../.." to the path
      strcat(pathtolilith, "../..");

      PyObject* sys_path = PySys_GetObject((char*)"path");
      if (sys_path == NULL) {
            printf("sys_path error\n");
            exit(1);
      }

      if(PyList_Append(sys_path, PyString_FromString(pathtolilith)) < 0) {
            printf("Error in appending sys_path and path_to_lh\n");
            exit(1);
      }


      // |------------------|
      // | Importing Lilith |
      // |------------------|

      // Defining lilith module
      PyObject* lilithModuleString = PyString_FromString((char*)"lilith");
      // Importing lilith module
      PyObject* lilithModule = PyImport_Import(lilithModuleString);

      if(lilithModule == NULL) {
            printf("Lilith was not found\n");
            printf("Please check the path to Lilith\n");
            exit(1);
      }


      // |-------------------------------------------------|
      // | Initialize Lilith and read experimental results |
      // |-------------------------------------------------|

      // Creating a Lilith object constructor
      PyObject* lilithconstructor = PyObject_GetAttrString(lilithModule,(char*)"Lilith");
      args = PyTuple_Pack(1, Py_False);
      // Creating lilithcalc: an object of the Lilith class
      PyObject* lilithcalc = PyObject_CallObject(lilithconstructor, args);

      // Getting the function that reads the experimental results
      PyObject* lilithsetexpinput = PyObject_GetAttrString(lilithcalc,(char*)"readexpinput");

      // Creating the argument for the previous function
      args = PyTuple_Pack(1, PyString_FromString((char*)experimental_input));

      // Calling the function
      // lastestLHC.list is used by default if no other list was provided
      if(strcmp(experimental_input, "") == 0) {
            PyObject_CallObject(lilithsetexpinput, NULL);
      } else {
            PyObject_CallObject(lilithsetexpinput, args);
      }

      // Checking if error has occured
      // If so, Lilith has no been initialized and cannot run: the code will exit
      if (PyErr_Occurred()) {
        printf("Error during the initialization of Lilith:\n\n");
        PyErr_PrintEx(0);
        printf("\nCode will now exit.\n");
        exit(1);
        }
      else{
        return lilithcalc;
        }
}


/**
Read the user-input XML string XMLinputstring
**/
PyObject* lilith_readuserinput(PyObject* lilithcalc, char* XMLinputstring)
{

      // Getting the function to read the user input
      PyObject* readuserinput = PyObject_GetAttrString(lilithcalc,(char*)"readuserinput");

      // Passing the user input file to this function
      PyObject* args;
      args = PyTuple_Pack(1, PyString_FromString(XMLinputstring));
      PyObject_CallObject(readuserinput, args);

      // Checking if error has occured
      if(PyErr_Occurred())
      {
        printf("Error during the reading of the user input:\n\n");
        PyErr_PrintEx(0);
        printf("\n-2LogL will be set to -1 when evaluated.\n\n");
        return lilithcalc;
      }
  
      else
      {
        userread = 1;
        return lilithcalc;
      }
}


/**
Read the user-input XML file located at XMLinputfile
**/
PyObject* lilith_readuserinput_fromfile(PyObject* lilithcalc, char* XMLinputfile)
{
  
      // Getting the function to read the user input
      PyObject* readuserinput = PyObject_GetAttrString(lilithcalc,(char*)"readuserinputfile");
      // Passing the user input file to this function
      PyObject* args;
      args = PyTuple_Pack(1, PyString_FromString(XMLinputfile));
      PyObject* lilithcalc_read = PyObject_CallObject(readuserinput, args);

      // Checking if error has occured
      if(PyErr_Occurred()){
        printf("Error during the reading of the user input\n\n");
        PyErr_PrintEx(0);
        printf("\nCode will now exit.\n");
        exit(1);
        }
  
      else{
        userread = 1;
        return lilithcalc_read;
        }
  
}


/**
Evaluate -2LogL
**/
float lilith_computelikelihood(PyObject* lilithcalc)
{
      if(userread==0){
          printf("\nError occured while reading the user input file\n");
          printf("-2LogL is set to -1.\n\n");
          return -1.;
      }
      // Getting the function that computes the likelihood
      PyObject* computelikelihood = PyObject_GetAttrString(lilithcalc,(char*)"computelikelihood");
      PyObject_CallObject(computelikelihood, NULL);
      // Extracting the likelihood
      PyObject* likelihood = PyObject_GetAttrString(lilithcalc,(char*)"l");
      float my_likelihood = PyFloat_AsDouble(likelihood);
  
      // Checking if error has occured
      // If so, the -2logL is set to -1 and the code keeps on running
      if(PyErr_Occurred()){
        printf("Error during the computation of the likelihood:\n\n");
        PyErr_PrintEx(0);
        printf("\n-2LogL is set to -1 for this point.\n\n");
        return -1.;
        }
  
      else{
        return my_likelihood;
        }
  
}


/**
Getting the experimental number of degrees of freedom
**/
float lilith_exp_ndf(PyObject* lilithcalc)
{
      // Computing the experimental ndf
      PyObject* ndf = PyObject_GetAttrString(lilithcalc,(char*)"exp_ndf");
      float my_ndf = PyFloat_AsDouble(ndf);
  
      // Checking if error has occured
      // If so, the ndf is set to -1 and the code keeps on running
      if(PyErr_Occurred()){
        printf("\nError during the computation of exp_ndf:\n\n");
        PyErr_PrintEx(0);
        printf("\nexp_ndf is set to -1 for this point.\n");
        return -1;
        }
  
      else{
        return my_ndf;
        }
  
}

/**
Writing the couplings in an output file
**/
void lilith_couplings_output(PyObject* lilithcalc, char* outputfilepath)
{

      if(userread==0){
          printf("\nError occured while reading the user input file\n");
          printf("Couplings output cannot be written.\n");
          return;
      }
      // Getting the output function
      PyObject* writecouplings = PyObject_GetAttrString(lilithcalc,(char*)"writecouplings");
      // Passing the output filepath to this function
      PyObject* args;
      args = PyTuple_Pack(1, PyString_FromString(outputfilepath));

      PyObject_CallObject(writecouplings, args);
  
      // Checking if error has occured
      if(PyErr_Occurred()){
        printf("\nError during the signal strength output writing:\n\n");
        PyErr_PrintEx(0);
        return;
        }
  
      else{
        return;
        }
  
}

/**
-2LogL output in XML (SLHA=0) or SLHA-like format (otherwise)
**/
void lilith_likelihood_output(PyObject* lilithcalc, char* outputfilepath, int slha)
{

      if(userread==0){
          printf("\nError occured while reading the user input file\n");
          if(slha==0){
          printf("XML likelihood output cannot be written.\n");
          }
          else
          {
          printf("SLHA likelihood output cannot be written.\n");
          }
          return;
      }
      // Getting the output function
      PyObject* writeresults = PyObject_GetAttrString(lilithcalc,(char*)"writeresults");
      // Passing the output filepath to this function
      PyObject* args;
      if(slha==0){
        args = PyTuple_Pack(2, PyString_FromString(outputfilepath), Py_False);
        }
      else{
        args = PyTuple_Pack(2, PyString_FromString(outputfilepath), Py_True);
        }
      PyObject_CallObject(writeresults, args);
  
      // Checking if error has occured
      if(PyErr_Occurred()){
        printf("\nError during the likelihood output writing:\n\n");
        PyErr_PrintEx(0);
        return;
        }
  
      else{
        return;
        }
  
}


/**
Signal strength output for individual Higgs(es) (tot=0) or
total signal strength (otherwise)
tot is not relevant is only one Higgs has been defined
**/
void lilith_mu_output(PyObject* lilithcalc, char* outputfilepath, int tot)
{

      if(userread==0){
          printf("\nError occured while reading the user input file\n");
          printf("Signal strengths output cannot be written.\n");
          return;
      }
  
      // Getting the output function
      PyObject* writesignalstrengths = PyObject_GetAttrString(lilithcalc,(char*)"writesignalstrengths");
      // Passing the output filepath to this function
      PyObject* args;
      if(tot==0){
        args = PyTuple_Pack(2, PyString_FromString(outputfilepath), Py_False);
        }
      else{
        args = PyTuple_Pack(2, PyString_FromString(outputfilepath), Py_True);
        }
      PyObject_CallObject(writesignalstrengths, args);
  
      // Checking if error has occured
      if(PyErr_Occurred()){
        printf("\nError during the signal strength output writing:\n\n");
        PyErr_PrintEx(0);
        return;
        }
  
      else{
        return;
        }
  
}

