// ===========================================================
//
// ReadByVariant.h: Read data variant by variant
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


namespace PySeqArray
{

using namespace Vectorization;


// =====================================================================

/// Object for reading basic variables variant by variant
class COREARRAY_DLL_LOCAL CApply_Variant_Basic: public CApply_Variant
{
protected:
	C_SVType SVType;
public:
	/// constructor
	CApply_Variant_Basic(CFileInfo &File, const char *var_name);
	virtual void ReadData(PyObject *val);
	virtual PyObject *NeedArray();
};


/// Object for reading positions variant by variant
class COREARRAY_DLL_LOCAL CApply_Variant_Pos: public CApply_Variant
{
protected:
	int *PtrPos;
	PyObject *VarNode;  ///< R object
public:
	/// constructor
	CApply_Variant_Pos(CFileInfo &File);
	virtual void ReadData(PyObject *val);
	virtual PyObject *NeedArray(int &nProtected);
};


/// Object for reading chromosomes variant by variant
class COREARRAY_DLL_LOCAL CApply_Variant_Chrom: public CApply_Variant
{
protected:
	CChromIndex *ChromIndex;
	PyObject *VarNode;  ///< R object
public:
	/// constructor
	CApply_Variant_Chrom(CFileInfo &File);
	virtual void ReadData(PyObject *val);
	virtual PyObject *NeedArray(int &nProtected);
};


// =====================================================================

/// Object for reading genotypes variant by variant
class COREARRAY_DLL_LOCAL CApply_Variant_Geno: public CApply_Variant
{
protected:
	CGenoIndex *GenoIndex;  ///< indexing genotypes
	ssize_t SiteCount;  ///< the total number of entries at a site
	ssize_t CellCount;  ///< the selected number of entries at a site
	vector<C_BOOL> Selection;  ///< the buffer of selection
	VEC_AUTO_PTR ExtPtr;       ///< a pointer to the additional buffer
	PyObject *VarIntGeno;      ///< genotype R integer object

	inline int _ReadGenoData(int *Base);
	inline C_UInt8 _ReadGenoData(C_UInt8 *Base);

public:
	ssize_t SampNum;  ///< the number of selected samples
	int Ploidy;       ///< ploidy

	/// constructor
	CApply_Variant_Geno();
	CApply_Variant_Geno(CFileInfo &File);
	~CApply_Variant_Geno();

	void Init(CFileInfo &File);

	virtual PyObject *NeedArray();
	virtual void ReadData(PyObject *val);

	/// read genotypes in 32-bit integer
	void ReadGenoData(int *Base);
	/// read genotypes in unsigned 8-bit intetger
	void ReadGenoData(C_UInt8 *Base);
};


// =====================================================================

/// Object for reading genotypes (dosages) variant by variant
class COREARRAY_DLL_LOCAL CApply_Variant_Dosage: public CApply_Variant_Geno
{
public:
	/// constructor
	CApply_Variant_Dosage(CFileInfo &File);

	virtual PyObject *NeedArray();
	virtual void ReadData(PyObject *val);

	/// read dosages in 32-bit integer
	void ReadDosage(int *Base);
	/// read dosages in unsigned 8-bit intetger
	void ReadDosage(C_UInt8 *Base);
};


// =====================================================================

/// Object for reading phasing information variant by variant
class COREARRAY_DLL_LOCAL CApply_Variant_Phase: public CApply_Variant
{
protected:
	ssize_t SiteCount;  ///< the total number of entries at a site
	ssize_t CellCount;  ///< the selected number of entries at a site
	bool UseRaw;  ///< whether use RAW type
	vector<C_BOOL> Selection;  ///< the buffer of selection
	PyObject *VarPhase;  ///< genotype R object

public:
	ssize_t SampNum;  ///< the number of selected samples
	int Ploidy;       ///< ploidy

	/// constructor
	CApply_Variant_Phase();
	CApply_Variant_Phase(CFileInfo &File, bool use_raw);

	void Init(CFileInfo &File, bool use_raw);

	virtual void ReadData(PyObject *val);
	virtual PyObject *NeedArray(int &nProtected);
};


// =====================================================================

/// Object for reading info variables variant by variant
class COREARRAY_DLL_LOCAL CApply_Variant_Info: public CApply_Variant
{
protected:
	CIndex *VarIndex;  ///< indexing the format variable
	C_SVType SVType;        ///< data type for GDS reading
	C_Int32 BaseNum;        ///< if 2-dim, the size of the first dimension
	map<int, PyObject*> VarList;  ///< a list of PyObject variables

public:
	/// constructor
	CApply_Variant_Info(CFileInfo &File, const char *var_name);

	virtual void ReadData(PyObject *val);
	virtual PyObject *NeedArray(int &nProtected);
};


// =====================================================================

/// Object for reading format variables variant by variant
class COREARRAY_DLL_LOCAL CApply_Variant_Format: public CApply_Variant
{
protected:
	CIndex *VarIndex;  ///< indexing the format variable
	ssize_t _TotalSampNum;  ///< the total number of samples

	C_SVType SVType;        ///< data type for GDS reading
	C_BOOL *SelPtr[2];      ///< pointers to selection
	map<int, PyObject*> VarList;  ///< a list of PyObject variables

public:
	ssize_t SampNum;  ///< the number of selected samples

	/// constructor
	CApply_Variant_Format();
	CApply_Variant_Format(CFileInfo &File, const char *var_name);

	void Init(CFileInfo &File, const char *var_name);

	virtual void ReadData(PyObject *val);
	virtual PyObject *NeedArray(int &nProtected);
};


// =====================================================================

/// Object for calculating the number of distinct alleles variant by variant
class COREARRAY_DLL_LOCAL CApply_Variant_NumAllele: public CApply_Variant
{
private:
	string strbuf;
public:
	/// constructor
	CApply_Variant_NumAllele(CFileInfo &File);

	virtual PyObject *NeedArray();
	virtual void ReadData(PyObject *val);
	int GetNumAllele();
};

}


extern "C"
{

/// Apply functions over margins on a working space
COREARRAY_DLL_EXPORT PyObject *SEQ_Apply_Variant(PyObject *gdsfile, PyObject *var_name,
	PyObject *FUN, PyObject *as_is, PyObject *var_index, PyObject *param, PyObject *rho);

} // extern "C"
