/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__dCaAP
#define _nrn_initial _nrn_initial__dCaAP
#define nrn_cur _nrn_cur__dCaAP
#define _nrn_current _nrn_current__dCaAP
#define nrn_jacob _nrn_jacob__dCaAP
#define nrn_state _nrn_state__dCaAP
#define _net_receive _net_receive__dCaAP 
#define dCaAP dCaAP__dCaAP 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define w _p[0]
#define vth _p[1]
#define refract_period _p[2]
#define tauA _p[3]
#define tauB _p[4]
#define D _p[5]
#define sigma_diff _p[6]
#define K _p[7]
#define vrest _p[8]
#define t_dCaAP _p[9]
#define dcaap_count _p[10]
#define i _p[11]
#define A _p[12]
#define B _p[13]
#define DA _p[14]
#define DB _p[15]
#define v _p[16]
#define _g _p[17]
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 /* declaration of user functions */
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_dCaAP", _hoc_setdata,
 0, 0
};
 /* declare global and static user variables */
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "vth_dCaAP", "mV",
 "i_dCaAP", "nA",
 0,0
};
 static double A0 = 0;
 static double B0 = 0;
 static double delta_t = 0.01;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void _ba1() ;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[0]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"dCaAP",
 "w_dCaAP",
 "vth_dCaAP",
 "refract_period_dCaAP",
 "tauA_dCaAP",
 "tauB_dCaAP",
 "D_dCaAP",
 "sigma_diff_dCaAP",
 "K_dCaAP",
 "vrest_dCaAP",
 "t_dCaAP_dCaAP",
 "dcaap_count_dCaAP",
 0,
 "i_dCaAP",
 0,
 "A_dCaAP",
 "B_dCaAP",
 0,
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 18, _prop);
 	/*initialize range parameters*/
 	w = 0;
 	vth = 0;
 	refract_period = 0;
 	tauA = 0;
 	tauB = 0;
 	D = 0;
 	sigma_diff = 0;
 	K = 0;
 	vrest = 0;
 	t_dCaAP = 0;
 	dcaap_count = 0;
 	_prop->param = _p;
 	_prop->param_size = 18;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 1, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _dCaAP_reg() {
	int _vectorized = 1;
  _initlists();
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 18, 1);
  hoc_register_dparam_semantics(_mechtype, 0, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_reg_ba(_mechtype, _ba1, 11);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 dCaAP D:/Dendritic action potentials and computation in/GidonEtAl2020_fig3andS9/NEURON_Code/_mod/dCaAP.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[2], _dlist1[2];
 static int dCaAP(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   DA = A * ( 1.0 - A ) / tauA ;
   DB = B * ( 1.0 - B ) / tauB ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 DA = DA  / (1. - dt*( ( (( 1.0 )*( ( 1.0 - A ) ) + ( A )*( ( ( - 1.0 ) ) )) ) / tauA )) ;
 DB = DB  / (1. - dt*( ( (( 1.0 )*( ( 1.0 - B ) ) + ( B )*( ( ( - 1.0 ) ) )) ) / tauB )) ;
  return 0;
}
 /*END CVODE*/
 static int dCaAP (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
    A = A + (1. - exp(dt*(( (( 1.0 )*( ( 1.0 - A ) ) + ( A )*( ( ( - 1.0 ) ) )) ) / tauA)))*( - A) ;
    B = B + (1. - exp(dt*(( (( 1.0 )*( ( 1.0 - B ) ) + ( B )*( ( ( - 1.0 ) ) )) ) / tauB)))*( - B) ;
   }
  return 0;
}
 /* BEFORE BREAKPOINT */
 static void _ba1(Node*_nd, double* _pp, Datum* _ppd, Datum* _thread, _NrnThread* _nt)  {
   double* _p; Datum* _ppvar; _p = _pp; _ppvar = _ppd;
  v = NODEV(_nd);
 if ( t >= t_dCaAP + refract_period  && v > vth  && w > 0.0 ) {
     t_dCaAP = t ;
     A = 0.001 ;
     B = 0.0 ;
     K = exp ( - ( v - vth ) / ( vth - vrest ) / D ) ;
     if ( K > 1.0 ) {
       K = 1.0 ;
       }
     dcaap_count = dcaap_count + 1.0 ;
     }
   if ( B  == 0.0  && t >= t_dCaAP + sigma_diff  && w > 0.0 ) {
     B = 0.001 ;
     }
   }
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 (_p, _ppvar, _thread, _nt);
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
 _ode_matsol_instance1(_threadargs_);
 }}

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  A = A0;
  B = B0;
 {
   A = 0.0 ;
   B = 0.0 ;
   K = 0.0 ;
   t_dCaAP = - refract_period ;
   dcaap_count = 0.0 ;
   }
 
}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
 initmodel(_p, _ppvar, _thread, _nt);
}
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   i = - ( A - B ) * w * K ;
   }
 _current += i;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
 	}
 _g = (_g - _rhs)/.001;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}
 
}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
 {   dCaAP(_p, _ppvar, _thread, _nt);
  }}}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(A) - _p;  _dlist1[0] = &(DA) - _p;
 _slist1[1] = &(B) - _p;  _dlist1[1] = &(DB) - _p;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "dCaAP.mod";
static const char* nmodl_file_text = 
  "\n"
  ": Novel dendritic action potentials shape the computational \n"
  ": 	properties of human layer 2/3 cortical neurons( Gidon et al., 2019)\n"
  ": \n"
  ": written by \n"
  ": Athanasia Papoutsi and Albert Gidon\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX dCaAP\n"
  "	RANGE w, vth\n"
  "	RANGE K, tauA,tauB, A, B, D\n"
  "	RANGE t_dCaAP\n"
  "	RANGE sigma_diff,refract_period\n"
  "	RANGE vrest : to normalize the K\n"
  "	RANGE dcaap_count\n"
  "	NONSPECIFIC_CURRENT i\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "    (mV) = (millivolt)\n"
  "    (nA) = (nanoamp)\n"
  "    (uS) = (microsiemens)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "    w = 0 : by default the mechanism does not work.\n"
  "    vth    (mV)   : spike threshold \n"
  "	refract_period\n"
  "	tauA\n"
  "	tauB\n"
  "	D\n"
  "	sigma_diff\n"
  "	K\n"
  "	vrest\n"
  "	t_dCaAP \n"
  "	dcaap_count\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "    v (mV)\n"
  "   	i (nA)\n"
  "}\n"
  "\n"
  "STATE{\n"
  "	A\n"
  "	B\n"
  "}\n"
  "INITIAL {\n"
  "	A=0\n"
  "	B=0\n"
  ": the first spike is a complication, because we only describe\n"
  ": the inactivation mechanism phenomemologically, we just skip \n"
  ": the first spike	\n"
  "	K=0 \n"
  ":just make sure we are not in the refractory period	\n"
  "	t_dCaAP = -refract_period \n"
  "	dcaap_count = 0\n"
  "}\n"
  "\n"
  "BREAKPOINT {	\n"
  "	SOLVE dCaAP METHOD cnexp\n"
  "	i =  -(A - B) * w * K\n"
  "}\n"
  "\n"
  "DERIVATIVE dCaAP{\n"
  "	: derivative of this \n"
  "	: 		sigmoid(y) = 1 / (1 + exp(-y))\n"
  "	: given y = (t-b)/tau\n"
  "	: is \n"
  "	: 		(1/tau)*sigmoid(y) * (1 - sigmoid(y) )\n"
  "	A' =  A * (1 - A) / tauA\n"
  "	B' =  B * (1 - B) / tauB\n"
  "}\n"
  "\n"
  "\n"
  "BEFORE BREAKPOINT {\n"
  "	: 1. activate dcaap after refractory period\n"
  "	: 2. activate dcaap only if threshold was crossed\n"
  "	: 3. ... and only in compartments where w is set to > 0\n"
  "	if (t >= t_dCaAP + refract_period && v > vth && w > 0) {\n"
  "		t_dCaAP = t\n"
  "		A=0.001\n"
  "		B=0	\n"
  "		: K is set only in the begining of the dcaap.\n"
  "		K = exp( -(v - vth) / (vth - vrest) / D ) \n"
  "		if(K > 1){\n"
  "			K = 1\n"
  "		}\n"
  "		dcaap_count = dcaap_count + 1\n"
  "		:printf(\"Fired %.0f dCaAP at %f ms with K = %f\\n\",dcaap_count, t_dCaAP,K)\n"
  "	}\n"
  "	: B is equal to 0 only when spike is initially triggered.\n"
  "	if(B == 0 && t >= t_dCaAP + sigma_diff && w > 0){\n"
  "		B=0.001\n"
  "	}		\n"
  "}\n"
  ;
#endif
