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

PyObject* initialize_lilith(char* );
PyObject* lilith_readuserinput(PyObject* , char*);
PyObject* lilith_readuserinput_fromfile(PyObject* , char*);

float lilith_computelikelihood(PyObject*);
float lilith_exp_ndf(PyObject*);
void lilith_likelihood_output(PyObject*, char*, int);
void lilith_mu_output(PyObject*, char*, int);
void lilith_couplings_output(PyObject*, char*);

