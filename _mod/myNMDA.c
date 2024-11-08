/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
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
 
#define nrn_init _nrn_init__myNMDA
#define _nrn_initial _nrn_initial__myNMDA
#define nrn_cur _nrn_cur__myNMDA
#define _nrn_current _nrn_current__myNMDA
#define nrn_jacob _nrn_jacob__myNMDA
#define nrn_state _nrn_state__myNMDA
#define _net_receive _net_receive__myNMDA 
#define EPSP EPSP__myNMDA 
#define play play__myNMDA 
#define remove_play remove_play__myNMDA 
#define state state__myNMDA 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define e _p[0]
#define tau_r _p[1]
#define tau_d _p[2]
#define dist _p[3]
#define factor _p[4]
#define gmax _p[5]
#define tsyn _p[6]
#define mg _p[7]
#define mggate _p[8]
#define gamma _p[9]
#define Vtrg _p[10]
#define tau_H _p[11]
#define space_synutils _p[12]
#define use_play_synutils _p[13]
#define i _p[14]
#define g _p[15]
#define A _p[16]
#define B _p[17]
#define lSF _p[18]
#define DA _p[19]
#define DB _p[20]
#define DlSF _p[21]
#define _g _p[22]
#define _tsav _p[23]
#define _nd_area  *_ppvar[0]._pval
 
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
 /* external NEURON variables */
 /* declaration of user functions */
 static double _hoc_EPSP();
 static double _hoc_get_factor();
 static double _hoc_play();
 static double _hoc_remove_play();
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

 extern Prop* nrn_point_prop_;
 static int _pointtype;
 static void* _hoc_create_pnt(_ho) Object* _ho; { void* create_point_process();
 return create_point_process(_pointtype, _ho);
}
 static void _hoc_destroy_pnt();
 static double _hoc_loc_pnt(_vptr) void* _vptr; {double loc_point_process();
 return loc_point_process(_pointtype, _vptr);
}
 static double _hoc_has_loc(_vptr) void* _vptr; {double has_loc_point();
 return has_loc_point(_vptr);
}
 static double _hoc_get_loc_pnt(_vptr)void* _vptr; {
 double get_loc_point_process(); return (get_loc_point_process(_vptr));
}
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata(void* _vptr) { Prop* _prop;
 _prop = ((Point_process*)_vptr)->_prop;
   _setdata(_prop);
 }
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 0,0
};
 static Member_func _member_func[] = {
 "loc", _hoc_loc_pnt,
 "has_loc", _hoc_has_loc,
 "get_loc", _hoc_get_loc_pnt,
 "EPSP", _hoc_EPSP,
 "get_factor", _hoc_get_factor,
 "play", _hoc_play,
 "remove_play", _hoc_remove_play,
 0, 0
};
#define get_factor get_factor_myNMDA
 extern double get_factor( double , double );
 /* declare global and static user variables */
#define mg_dep mg_dep_myNMDA
 double mg_dep = 3.57;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "e", "mV",
 "tau_r", "ms",
 "tau_d", "ms",
 "mg", "mM",
 "tau_H", "ms",
 "i", "nA",
 "g", "uS",
 0,0
};
 static double A0 = 0;
 static double B0 = 0;
 static double delta_t = 0.01;
 static double lSF0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "mg_dep_myNMDA", &mg_dep_myNMDA,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 static void _hoc_destroy_pnt(_vptr) void* _vptr; {
   destroy_point_process(_vptr);
}
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[3]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"myNMDA",
 "e",
 "tau_r",
 "tau_d",
 "dist",
 "factor",
 "gmax",
 "tsyn",
 "mg",
 "mggate",
 "gamma",
 "Vtrg",
 "tau_H",
 "space_synutils",
 "use_play_synutils",
 0,
 "i",
 "g",
 0,
 "A",
 "B",
 "lSF",
 0,
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
  if (nrn_point_prop_) {
	_prop->_alloc_seq = nrn_point_prop_->_alloc_seq;
	_p = nrn_point_prop_->param;
	_ppvar = nrn_point_prop_->dparam;
 }else{
 	_p = nrn_prop_data_alloc(_mechtype, 24, _prop);
 	/*initialize range parameters*/
 	e = 0;
 	tau_r = 0;
 	tau_d = 0;
 	dist = 0;
 	factor = 0;
 	gmax = 0;
 	tsyn = 0;
 	mg = 0;
 	mggate = 0;
 	gamma = 0.062;
 	Vtrg = -60;
 	tau_H = 1e+300;
 	space_synutils = 0;
 	use_play_synutils = 0;
  }
 	_prop->param = _p;
 	_prop->param_size = 24;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 
#define _tqitem &(_ppvar[2]._pvoid)
 static void _net_receive(Point_process*, double*, double);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _myNMDA_reg() {
	int _vectorized = 0;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 0,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 24, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "netsend");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_size[_mechtype] = 1;
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 myNMDA D:/Dendritic action potentials and computation in/GidonEtAl2020_fig3andS9/NEURON_Code/_mod/myNMDA.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int EPSP(double);
static int play();
static int remove_play();
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[3], _dlist1[3];
 static int state(_threadargsproto_);
 
double get_factor (  double _ltau_fast , double _ltau_slow ) {
   double _lget_factor;
 double _ltp ;
 _ltp = ( _ltau_fast * _ltau_slow ) / ( _ltau_slow - _ltau_fast ) * log ( _ltau_slow / _ltau_fast ) ;
   _lget_factor = - exp ( - _ltp / _ltau_fast ) + exp ( - _ltp / _ltau_slow ) ;
   _lget_factor = 1.0 / _lget_factor ;
   
return _lget_factor;
 }
 
static double _hoc_get_factor(void* _vptr) {
 double _r;
    _hoc_setdata(_vptr);
 _r =  get_factor (  *getarg(1) , *getarg(2) );
 return(_r);
}
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   DA = - A / tau_r ;
   DB = - B / tau_d ;
   DlSF = ( Vtrg - v ) / tau_H ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 DA = DA  / (1. - dt*( ( - 1.0 ) / tau_r )) ;
 DB = DB  / (1. - dt*( ( - 1.0 ) / tau_d )) ;
 DlSF = DlSF  / (1. - dt*( 0.0 )) ;
  return 0;
}
 /*END CVODE*/
 static int state () {_reset=0;
 {
    A = A + (1. - exp(dt*(( - 1.0 ) / tau_r)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau_r ) - A) ;
    B = B + (1. - exp(dt*(( - 1.0 ) / tau_d)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau_d ) - B) ;
    lSF = lSF - dt*(- ( ( ( Vtrg - v ) ) / tau_H ) ) ;
   }
  return 0;
}
 
static int  EPSP (  double _lw ) {
   A = A + _lw * factor ;
   B = B + _lw * factor ;
   tsyn = t ;
    return 0; }
 
static double _hoc_EPSP(void* _vptr) {
 double _r;
    _hoc_setdata(_vptr);
 _r = 1.;
 EPSP (  *getarg(1) );
 return(_r);
}
 
/*VERBATIM*/

extern double* vector_vec(); 
extern int vector_capacity(); 
extern void* vector_arg(); 
 
static int  play (  ) {
   use_play_synutils = 1.0 ;
   
/*VERBATIM*/

	void** vv; 
	vv = (void**)(&space_synutils); 
	*vv = (void*)0; 
	if (ifarg(1)) { 
	 *vv = vector_arg(1); 
	} 
  return 0; }
 
static double _hoc_play(void* _vptr) {
 double _r;
    _hoc_setdata(_vptr);
 _r = 1.;
 play (  );
 return(_r);
}
 
static int  remove_play (  ) {
   use_play_synutils = 0.0 ;
    return 0; }
 
static double _hoc_remove_play(void* _vptr) {
 double _r;
    _hoc_setdata(_vptr);
 _r = 1.;
 remove_play (  );
 return(_r);
}
 
static void _net_receive (_pnt, _args, _lflag) Point_process* _pnt; double* _args; double _lflag; 
{    _p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
  if (_tsav > t){ extern char* hoc_object_name(); hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t;   if (_lflag == 1. ) {*(_tqitem) = 0;}
 {
   EPSP ( _threadargscomma_ 1.0 ) ;
   } }
 
static int _ode_count(int _type){ return 3;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 ();
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 3; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 ();
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
 _ode_matsol_instance1(_threadargs_);
 }}

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  A = A0;
  B = B0;
  lSF = lSF0;
 {
   A = 0.0 ;
   B = 0.0 ;
   g = 0.0 ;
   lSF = 1.0 ;
   tsyn = - 1e80 ;
   if ( tau_r / tau_d > 0.999 ) {
     tau_r = tau_d * 0.999 ;
     }
   if (  ! tau_r  ||  ! tau_d ) {
     printf ( "User must set tau_r and tau_d (zero by default)" ) ;
     }
   factor = get_factor ( _threadargscomma_ tau_r , tau_d ) ;
   }
 {
   double _letime ;
 
/*VERBATIM*/
	void* vv1; int size; double *px; int ii;
 	vv1 = *((void**)(&space_synutils)); 
	if (vv1 && use_play_synutils) { 
	 	size  = vector_capacity(vv1);
		px = vector_vec(vv1); 
		for(ii=0;ii<size;ii++){
			_letime = px[ii];
 net_send ( _tqitem, (double*)0, _ppvar[1]._pvoid, t +  _letime , - 1.0 ) ;
   
/*VERBATIM*/
		} 
	}
 }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _tsav = -1e20;
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
 initmodel();
}}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   mggate = 1.0 / ( 1.0 + exp ( gamma * - v ) * ( mg / mg_dep ) ) ;
   lSF = lSF * ( lSF > 0.0 ) ;
   g = gmax * ( B - A ) * mggate * lSF ;
   i = g * ( v - e ) ;
   }
 _current += i;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
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
 _g = _nrn_current(_v + .001);
 	{ _rhs = _nrn_current(_v);
 	}
 _g = (_g - _rhs)/.001;
 _g *=  1.e2/(_nd_area);
 _rhs *= 1.e2/(_nd_area);
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
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
 
}}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
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
 { error =  state();
 if(error){fprintf(stderr,"at line 81 in file ISyn.inc:\n	: from Jahr & Stevens\n"); nrn_complain(_p); abort_run(error);}
 }}}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(A) - _p;  _dlist1[0] = &(DA) - _p;
 _slist1[1] = &(B) - _p;  _dlist1[1] = &(DB) - _p;
 _slist1[2] = &(lSF) - _p;  _dlist1[2] = &(DlSF) - _p;
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "myNMDA.mod";
static const char* nmodl_file_text = 
  "\n"
  "COMMENT\n"
  "\n"
  "\n"
  "Use only what is in ISyn and IUtils\n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON {\n"
  "	POINT_PROCESS myNMDA\n"
  "}\n"
  "\n"
  ":::INCLUDE \"ISyn.inc\"\n"
  "\n"
  "NEURON {\n"
  "	:NOTE: since this is an interface, there is no POINT_PROCESS name\n"
  "	RANGE e, i, gmax, g\n"
  "	RANGE tau_r\n"
  "	RANGE tau_d\n"
  "	RANGE dist :distance from the soma\n"
  "	RANGE ndist :normalized distance from the soma\n"
  "	RANGE factor\n"
  "	RANGE tsyn,mg,mggate,gamma\n"
  "	GLOBAL mg_dep\n"
  "	\n"
  "	RANGE lSF,gSF,Vtrg	:local and global synaptic homeostatic plasticity factor between 0-1\n"
  "	RANGE tau_H	: time constant for the homeostasis, should be orders of magnitude larger than the membrane time constant.\n"
  "\n"
  "	NONSPECIFIC_CURRENT i\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(nA) = (nanoamp)\n"
  "	(mV) = (millivolt)\n"
  "	(uS) = (microsiemens)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	e = 0			(mV)	\n"
  "	tau_r = 0	 	(ms) 	\n"
  "	tau_d = 0		(ms) \n"
  "	dist\n"
  "	factor\n"
  "	gmax\n"
  "	tsyn\n"
  "	mg			(mM)\n"
  "	mggate\n"
  "	gamma = 0.062\n"
  "	Vtrg = -60\n"
  "	:the default is that there is not homeostasis with such a high time constant\n"
  "	tau_H = 1e300	 	(ms) :assuming membrane potential within the range of 20ms tau_H is three orders of magnitude higher. 	\n"
  "	mg_dep = 3.57\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	v (mV)\n"
  "	i (nA)\n"
  "	g (uS)\n"
  "}\n"
  "\n"
  "STATE {\n"
  "	A\n"
  "	B\n"
  "	lSF\n"
  "}\n"
  "\n"
  "FUNCTION get_factor(tau_fast,tau_slow){\n"
  "	LOCAL tp\n"
  "	tp = (tau_fast*tau_slow)/(tau_slow - tau_fast) * log(tau_slow/tau_fast)\n"
  "	get_factor = -exp(-tp/tau_fast) + exp(-tp/tau_slow)\n"
  "	get_factor = 1/get_factor\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	:start fresh\n"
  "	A = 0\n"
  "	B = 0\n"
  "	g = 0\n"
  "	lSF = 1\n"
  "	tsyn = -1e80 			:last time the synapse was active is very long ego\n"
  "	if(tau_r/tau_d>0.999){\n"
  "		tau_r = tau_d*0.999\n"
  "	}\n"
  "	if(!tau_r || !tau_d){\n"
  "		printf(\"User must set tau_r and tau_d (zero by default)\")\n"
  "	}\n"
  "	factor = get_factor(tau_r,tau_d)\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE state METHOD cnexp\n"
  "	:see Jahr & Stevens, J. Neurosci 10: 1830-1837, 1990;\n"
  "	:Jahr & Stevens, J. Neurosci 10: 3178-3182, 1990\n"
  "	: from Jahr & Stevens\n"
  "	mggate = 1 / (1 + exp( gamma *  -v) * (mg / mg_dep))\n"
  "	lSF = lSF*(lSF>0)\n"
  "	g = gmax*(B - A)*mggate*lSF\n"
  "	i = g*(v - e)\n"
  "}\n"
  "\n"
  "\n"
  "DERIVATIVE state {\n"
  "	A' = -A/tau_r\n"
  "	B' = -B/tau_d\n"
  "	lSF' = (Vtrg-v)/tau_H :the global implementation requires a pointer to the somatic voltage.\n"
  "}\n"
  "\n"
  "\n"
  "PROCEDURE EPSP(w){\n"
  "    \n"
  "	A = A + w*factor\n"
  "	B = B + w*factor\n"
  "	tsyn = t\n"
  "}\n"
  ":::end INCLUDE ISyn.inc\n"
  ":::INCLUDE \"IUtils.inc\"\n"
  "NEURON{\n"
  "	: NOTE that there is no POINT_PROCESS name since this is an interface\n"
  "	RANGE space_synutils\n"
  "	RANGE use_play_synutils\n"
  "}\n"
  "\n"
  "VERBATIM \n"
  "extern double* vector_vec(); \n"
  "extern int vector_capacity(); \n"
  "extern void* vector_arg(); \n"
  "ENDVERBATIM      \n"
  " \n"
  "\n"
  "PARAMETER{\n"
  "	space_synutils	\n"
  "	use_play_synutils\n"
  "}\n"
  "\n"
  "\n"
  "PROCEDURE play() { \n"
  "	use_play_synutils = 1\n"
  "	VERBATIM \n"
  "	void** vv; \n"
  "	vv = (void**)(&space_synutils); \n"
  "	*vv = (void*)0; \n"
  "	if (ifarg(1)) { \n"
  "	 *vv = vector_arg(1); \n"
  "	} \n"
  "	ENDVERBATIM \n"
  "} \n"
  "\n"
  "PROCEDURE remove_play(){\n"
  "	use_play_synutils = 0\n"
  "}\n"
  "\n"
  "INITIAL{LOCAL etime\n"
  "VERBATIM\n"
  "	void* vv1; int size; double *px; int ii;\n"
  " 	vv1 = *((void**)(&space_synutils)); \n"
  "	if (vv1 && use_play_synutils) { \n"
  "	 	size  = vector_capacity(vv1);\n"
  "		px = vector_vec(vv1); \n"
  "		for(ii=0;ii<size;ii++){\n"
  "			_letime = px[ii];\n"
  "			ENDVERBATIM \n"
  "				net_send (etime,-1) \n"
  "				:printf(\"%1f etime=\",etime)\n"
  "			VERBATIM\n"
  "		} \n"
  "	}\n"
  "ENDVERBATIM\n"
  "}\n"
  ":::end INCLUDE IUtils.inc\n"
  "\n"
  "NET_RECEIVE(w) {\n"
  "	EPSP(1): finally trigger an epsp\n"
  "}\n"
  ;
#endif
