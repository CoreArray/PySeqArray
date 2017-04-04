// ===========================================================
//
// GetData.cpp: Get data from the GDS file
//
// Copyright (C) 2017    Xiuwen Zheng
//
// This file is part of PySeqArray.
//
// PySeqArray is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License Version 3 as
// published by the Free Software Foundation.
//
// PySeqArray is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with PySeqArray.
// If not, see <http://www.gnu.org/licenses/>.

#include "Index.h"
#include "ReadByVariant.h"
// #include "ReadBySample.h"


using namespace PySeqArray;

extern "C"
{

// ===========================================================
// Get data from a working space
// ===========================================================

static bool is_logical(PdGDSObj Node)
{
	char classname[32];
	classname[0] = 0;
	GDS_Node_GetClassName(Node, classname, sizeof(classname));
	return (strcmp(classname, "dBit1") == 0);
}


// get data
static PyObject* VarGetData(CFileInfo &File, const char *name)
{
	static const char *ERR_DIM = "Invalid dimension of '%s'.";

	PyObject *rv_ans = NULL;
	TSelection &Sel = File.Selection();

	if (strcmp(name, "sample.id") == 0)
	{
		// ===========================================================
		// sample.id

		PdAbstractArray N = File.GetObj(name, TRUE);
		// check
		if ((GDS_Array_DimCnt(N) != 1) ||
				(GDS_Array_GetTotalCount(N) != File.SampleNum()))
			throw ErrSeqArray(ERR_DIM, name);
		// read
		C_BOOL *ss = Sel.pSample();
		rv_ans = GDS_Py_Array_Read(N, NULL, NULL, &ss, svCustom);

	} else if (strcmp(name, "position") == 0)
	{
		int n = File.VariantSelNum();
		rv_ans = numpy_new_int32(n);
		if (n > 0)
		{
			const int *base = &File.Position()[0];
			int *p = (int*)numpy_getptr(rv_ans);
			C_BOOL *s = Sel.pVariant();
			for (size_t m=File.VariantNum(); m > 0; m--)
			{
				if (*s++) *p++ = *base;
				base ++;
			}
		}

	} else if (strcmp(name, "chromosome") == 0)
	{
		int n = File.VariantSelNum();
		rv_ans = numpy_new_string(n);
		if (n > 0)
		{
			CChromIndex &Chrom = File.Chromosome();
			PyObject **p = (PyObject**)numpy_getptr(rv_ans);
			C_BOOL *s = Sel.pVariant();
			size_t m = File.VariantNum();
			string lastss;
			PyObject *last = NULL;
			for (size_t i=0; i < m; i++)
			{
				if (*s++)
				{
					const string &ss = Chrom[i];
					if (ss != lastss)
					{
						lastss = ss;
						last = NULL;
					}
					if (!last)
						last = PYSTR_SET2(&lastss[0], lastss.size());
					numpy_setval(rv_ans, p++, last);
				}
			}
		}
	
	} else if ( (strcmp(name, "variant.id")==0) ||
		(strcmp(name, "allele")==0) ||
		(strcmp(name, "annotation/id")==0) ||
		(strcmp(name, "annotation/qual")==0) ||
		(strcmp(name, "annotation/filter")==0) )
	{
		// ===========================================================
		// variant.id, allele, annotation/id, annotation/qual, annotation/filter

		PdAbstractArray N = File.GetObj(name, TRUE);
		// check
		if ((GDS_Array_DimCnt(N) != 1) ||
				(GDS_Array_GetTotalCount(N) != File.VariantNum()))
			throw ErrSeqArray(ERR_DIM, name);
		// read
		C_BOOL *ss = Sel.pVariant();
		rv_ans = GDS_Py_Array_Read(N, NULL, NULL, &ss, svCustom);

	} else if (strcmp(name, "genotype") == 0)
	{
		// ===========================================================
		// genotypic data

		int nSample  = File.SampleSelNum();
		int nVariant = File.VariantSelNum();

		if ((nSample > 0) && (nVariant > 0))
		{
			// initialize GDS genotype Node
			CApply_Variant_Geno NodeVar(File);
			// set
			rv_ans = numpy_new_uint8_dim3(nVariant, nSample, File.Ploidy());
			C_UInt8 *base = (C_UInt8*)numpy_getptr(rv_ans);
			ssize_t SIZE = (ssize_t)nSample * File.Ploidy();
			do {
				NodeVar.ReadGenoData(base);
				base += SIZE;
			} while (NodeVar.Next());
		} else
			rv_ans = numpy_new_uint8(0);

	} else if (strcmp(name, "@genotype") == 0)
	{
		static const char *VarName = "genotype/@data";
		PdAbstractArray N = File.GetObj(VarName, TRUE);
		// check
		if ((GDS_Array_DimCnt(N) != 1) ||
				(GDS_Array_GetTotalCount(N) != File.VariantNum()))
			throw ErrSeqArray(ERR_DIM, VarName);
		// read
		C_BOOL *ss = Sel.pVariant();
		rv_ans = GDS_Py_Array_Read(N, NULL, NULL, &ss, svInt32);

	} else if (strcmp(name, "$dosage") == 0)
	{
		// ===========================================================
		// dosage data

		ssize_t nSample  = File.SampleSelNum();
		ssize_t nVariant = File.VariantSelNum();

		if ((nSample > 0) && (nVariant > 0))
		{
			// initialize GDS genotype Node
			CApply_Variant_Dosage NodeVar(File);
			// set
			rv_ans = numpy_new_uint8_mat(nVariant, nSample);
			C_UInt8 *base = (C_UInt8*)numpy_getptr(rv_ans);
			do {
				NodeVar.ReadDosage(base);
				base += nSample;
			} while (NodeVar.Next());
		} else
			rv_ans = numpy_new_uint8(0);

	} else if (strcmp(name, "phase") == 0)
	{
		// ===========================================================
		// phase/

		PdAbstractArray N = File.GetObj("phase/data", TRUE);
		// check
		int ndim = GDS_Array_DimCnt(N);
		C_Int32 dim[4];
		GDS_Array_GetDim(N, dim, 3);
		if (ndim<2 || ndim>3 || dim[0]!= File.VariantNum() ||
				dim[1]!=File.SampleNum())
			throw ErrSeqArray(ERR_DIM, name);
		// read
		C_BOOL *ss[3] = { Sel.pVariant(), Sel.pSample(), NULL };
		if (ndim == 3)
			ss[2] = NeedArrayTRUEs(dim[2]);
		rv_ans = GDS_Py_Array_Read(N, NULL, NULL, ss, svCustom);

	} else if (strncmp(name, "annotation/info/@", 17) == 0)
	{
		if (File.GetObj(name, FALSE) != NULL)
		{
			CIndex &V = File.VarIndex(name);
			rv_ans = V.GetLen_Sel(Sel.pVariant());
		}

	} else if (strncmp(name, "annotation/info/", 16) == 0)
	{
		// ===========================================================
		// annotation/info

		GDS_PATH_PREFIX_CHECK(name);
		PdAbstractArray N = File.GetObj(name, TRUE);
		int ndim = GDS_Array_DimCnt(N);
		if ((ndim!=1) && (ndim!=2))
			throw ErrSeqArray(ERR_DIM, name);

		string name2 = GDS_PATH_PREFIX(name, '@');
		PdAbstractArray N_idx = File.GetObj(name2.c_str(), FALSE);
		if (N_idx == NULL)
		{
			// no index
			C_Int32 dim[4];
			GDS_Array_GetDim(N, dim, 2);
			C_BOOL *ss[2] = { Sel.pVariant(), NULL };
			if (ndim == 2)
				ss[1] = NeedArrayTRUEs(dim[1]);
			C_SVType SV = svCustom;  // is_logical
			rv_ans = GDS_Py_Array_Read(N, NULL, NULL, ss, SV);

		} else {
			// with index
			CIndex &V = File.VarIndex(name2);
			int var_start, var_count;
			vector<C_BOOL> var_sel;
			PyObject *Index = V.GetLen_Sel(Sel.pVariant(), var_start, var_count, var_sel);

			C_BOOL *ss[2] = { &var_sel[0], NULL };
			C_Int32 dimst[2]  = { var_start, 0 };
			C_Int32 dimcnt[2] = { var_count, 0 };
			if (ndim == 2)
			{
				GDS_Array_GetDim(N, dimcnt, 2);
				dimcnt[0] = var_count;
			}
			PyObject *Val = GDS_Py_Array_Read(N, dimst, dimcnt, ss, svCustom);

			rv_ans = Py_BuildValue("{s:N,s:N}", "index", Index, "data", Val);
		}

	} else if (strncmp(name, "annotation/format/@", 19) == 0)
	{
		string name2(name);
		name2.erase(18, 1).append("/@data");
		if (File.GetObj(name2.c_str(), FALSE) != NULL)
		{
			CIndex &V = File.VarIndex(name2.c_str());
			rv_ans = V.GetLen_Sel(Sel.pVariant());
		}

	} else if (strncmp(name, "annotation/format/", 18) == 0)
	{
		// ===========================================================
		// annotation/format

		GDS_PATH_PREFIX_CHECK(name);
		string name1 = string(name) + "/data";
		string name2 = string(name) + "/@data";
		PdAbstractArray N = File.GetObj(name1.c_str(), TRUE);

		// with index
		CIndex &V = File.VarIndex(name2);
		int var_start, var_count;
		vector<C_BOOL> var_sel;
		PyObject *Index = V.GetLen_Sel(Sel.pVariant(), var_start, var_count, var_sel);

		C_BOOL *ss[2] = { &var_sel[0], Sel.pSample() };
		C_Int32 dimst[2]  = { var_start, 0 };
		C_Int32 dimcnt[2];
		GDS_Array_GetDim(N, dimcnt, 2);
		dimcnt[0] = var_count;
		PyObject *Val = GDS_Py_Array_Read(N, dimst, dimcnt, ss, svCustom);

		rv_ans = Py_BuildValue("{s:N,s:N}", "index", Index, "data", Val);

	} else if (strncmp(name, "sample.annotation/", 18) == 0)
	{
		// ===========================================================
		// sample.annotation

		GDS_PATH_PREFIX_CHECK(name);
		PdAbstractArray N = File.GetObj(name, TRUE);
		// check
		int ndim = GDS_Array_DimCnt(N);
		if ((ndim!=1) && (ndim!=2))
			throw ErrSeqArray(ERR_DIM, name);
		C_Int32 dim[2];
		GDS_Array_GetDim(N, dim, 2);
		if (dim[0] != File.SampleNum())
			throw ErrSeqArray(ERR_DIM, name);

		C_BOOL *ss[2] = { Sel.pSample(), NULL };
		if (ndim == 2)
			ss[1] = NeedArrayTRUEs(dim[1]);
		rv_ans = GDS_Py_Array_Read(N, NULL, NULL, ss, svCustom);

	} else if (strcmp(name, "$chrom_pos") == 0)
	{
		// ===========================================================
		// chromosome-position

		PdAbstractArray N1 = File.GetObj("chromosome", TRUE);
		PdAbstractArray N2 = File.GetObj("position", TRUE);
		C_Int64 n1 = GDS_Array_GetTotalCount(N1);
		C_Int64 n2 = GDS_Array_GetTotalCount(N2);
		if ((n1 != n2) || (n1 != File.VariantNum()))
			throw ErrSeqArray("Invalid dimension of 'chromosome' and 'position'.");

		vector<string> chr;
		vector<C_Int32> pos;

		int n = File.VariantSelNum();
		chr.resize(n);
		pos.resize(n);
		C_BOOL *ss = Sel.pVariant();

		GDS_Array_ReadDataEx(N1, NULL, NULL, &ss, &chr[0], svStrUTF8);
		GDS_Array_ReadDataEx(N2, NULL, NULL, &ss, &pos[0], svInt32);

		char buf1[1024] = { 0 };
		char buf2[1024] = { 0 };
		char *p1 = buf1, *p2 = buf2;
		int dup = 0;
		rv_ans = numpy_new_string(n1);
		PyObject **p = (PyObject**)numpy_getptr(rv_ans);
		for (size_t i=0; i < (size_t)n1; i++,p++)
		{
			snprintf(p1, sizeof(buf1), "%s_%d", chr[i].c_str(), pos[i]);
			if (strcmp(p1, p2) == 0)
			{
				dup ++;
				snprintf(p1, sizeof(buf1), "%s_%d_%d", chr[i].c_str(),
					pos[i], dup);
				numpy_setval(rv_ans, p, PYSTR_SET(p1));
			} else {
				char *tmp;
				tmp = p1; p1 = p2; p2 = tmp;
				numpy_setval(rv_ans, p, PYSTR_SET(p2));
				dup = 0;
			}
		}

	} else if (strcmp(name, "$num_allele") == 0)
	{
		// ===========================================================
		// the number of distinct alleles

		ssize_t nVariant = File.VariantSelNum();
		rv_ans = numpy_new_int32(nVariant);
		int *p = (int*)numpy_getptr(rv_ans);

		CApply_Variant_NumAllele NodeVar(File);
		for (ssize_t i=0; i < nVariant; i++)
		{
			p[i] = NodeVar.GetNumAllele();
			NodeVar.Next();
		}

	} else {
		throw ErrSeqArray(
			"'%s' is not a standard variable name, and the standard format:\n"
			"    sample.id, variant.id, position, chromosome, allele, genotype\n"
			"    annotation/id, annotation/qual, annotation/filter\n"
			"    annotation/info/VARIABLE_NAME, annotation/format/VARIABLE_NAME\n"
			"    sample.annotation/VARIABLE_NAME", name);
	}

	return rv_ans;
}


/// Get data from a working space
COREARRAY_DLL_EXPORT PyObject* SEQ_GetData(PyObject *self, PyObject *args)
{
	int file_id;
	const char *name;
	if (!PyArg_ParseTuple(args, "is", &file_id, &name))
		return NULL;

	COREARRAY_TRY
		// File information
		CFileInfo &File = GetFileInfo(file_id);
		// Get data
		return VarGetData(File, name);
	COREARRAY_CATCH_NONE
}

} // extern "C"
