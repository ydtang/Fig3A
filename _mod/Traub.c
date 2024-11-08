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
 
#define nrn_init _nrn_init__traub
#define _nrn_initial _nrn_initial__traub
#define nrn_cur _nrn_cur__traub
#define _nrn_current _nrn_current__traub
#define nrn_jacob _nrn_jacob__traub
#define nrn_state _nrn_state__traub
#define _net_receive _net_receive__traub 
#define _f_rates _f_rates__traub 
#define rates rates__traub 
#define states states__traub 
 
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
#define gNabar _p[0]
#define gKbar _p[1]
#define gLbar _p[2]
#define eL _p[3]
#define eK _p[4]
#define eNa _p[5]
#define v_shft _p[6]
#define i _p[7]
#define iL _p[8]
#define iNa _p[9]
#define iK _p[10]
#define m _p[11]
#define h _p[12]
#define n _p[13]
#define a _p[14]
#define b _p[15]
#define Dm _p[16]
#define Dh _p[17]
#define Dn _p[18]
#define Da _p[19]
#define Db _p[20]
#define cm _p[21]
#define gNa _p[22]
#define gK _p[23]
#define minf _p[24]
#define hinf _p[25]
#define ninf _p[26]
#define mtau _p[27]
#define htau _p[28]
#define ntau _p[29]
#define v _p[30]
#define _g _p[31]
 
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
 static void _hoc_rates(void);
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
 "setdata_traub", _hoc_setdata,
 "rates_traub", _hoc_rates,
 0, 0
};
 
static void _check_rates(double*, Datum*, Datum*, _NrnThread*); 
static void _check_table_thread(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, int _type) {
   _check_rates(_p, _ppvar, _thread, _nt);
 }
 /* declare global and static user variables */
 static int _thread1data_inuse = 0;
static double _thread1data[1];
#define _gth 0
#define totG_traub _thread1data[0]
#define totG _thread[_gth]._pval[0]
#define usetable usetable_traub
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_traub", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "gNabar_traub", "S/cm2",
 "gKbar_traub", "S/cm2",
 "gLbar_traub", "S/cm2",
 "eL_traub", "mV",
 "eK_traub", "mV",
 "eNa_traub", "mV",
 "i_traub", "mA/cm2",
 "iK_traub", "mA/cm2",
 0,0
};
 static double a0 = 0;
 static double b0 = 0;
 static double delta_t = 0.01;
 static double h0 = 0;
 static double m0 = 0;
 static double n0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "totG_traub", &totG_traub,
 "usetable_traub", &usetable_traub,
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
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[0]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"traub",
 "gNabar_traub",
 "gKbar_traub",
 "gLbar_traub",
 "eL_traub",
 "eK_traub",
 "eNa_traub",
 "v_shft_traub",
 0,
 "i_traub",
 "iL_traub",
 "iNa_traub",
 "iK_traub",
 0,
 "m_traub",
 "h_traub",
 "n_traub",
 "a_traub",
 "b_traub",
 0,
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 32, _prop);
 	/*initialize range parameters*/
 	gNabar = 0.03;
 	gKbar = 0.015;
 	gLbar = 0.00014;
 	eL = -62;
 	eK = -80;
 	eNa = 90;
 	v_shft = 49.2;
 	_prop->param = _p;
 	_prop->param_size = 32;
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
 static void _thread_mem_init(Datum*);
 static void _thread_cleanup(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _Traub_reg() {
	int _vectorized = 1;
  _initlists();
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 2);
  _extcall_thread = (Datum*)ecalloc(1, sizeof(Datum));
  _thread_mem_init(_extcall_thread);
  _thread1data_inuse = 0;
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 1, _thread_mem_init);
     _nrn_thread_reg(_mechtype, 0, _thread_cleanup);
     _nrn_thread_table_reg(_mechtype, _check_table_thread);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 32, 1);
  hoc_register_dparam_semantics(_mechtype, 0, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 traub D:/Dendritic action potentials and computation in/GidonEtAl2020_fig3andS9/NEURON_Code/_mod/Traub.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double *_t_mtau;
 static double *_t_ntau;
 static double *_t_htau;
 static double *_t_minf;
 static double *_t_ninf;
 static double *_t_hinf;
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int _f_rates(_threadargsprotocomma_ double);
static int rates(_threadargsprotocomma_ double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static void _n_rates(_threadargsprotocomma_ double _lv);
 static int _slist1[3], _dlist1[3];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   rates ( _threadargscomma_ v ) ;
   Dm = ( minf - m ) / mtau ;
   Dh = ( hinf - h ) / htau ;
   Dn = ( ninf - n ) / ntau ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 rates ( _threadargscomma_ v ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / mtau )) ;
 Dh = Dh  / (1. - dt*( ( ( ( - 1.0 ) ) ) / htau )) ;
 Dn = Dn  / (1. - dt*( ( ( ( - 1.0 ) ) ) / ntau )) ;
  return 0;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
   rates ( _threadargscomma_ v ) ;
    m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / mtau)))*(- ( ( ( minf ) ) / mtau ) / ( ( ( ( - 1.0 ) ) ) / mtau ) - m) ;
    h = h + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / htau)))*(- ( ( ( hinf ) ) / htau ) / ( ( ( ( - 1.0 ) ) ) / htau ) - h) ;
    n = n + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / ntau)))*(- ( ( ( ninf ) ) / ntau ) / ( ( ( ( - 1.0 ) ) ) / ntau ) - n) ;
   }
  return 0;
}
 static double _mfac_rates, _tmin_rates;
  static void _check_rates(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  if (!usetable) {return;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_rates =  - 100.0 ;
   _tmax =  70.0 ;
   _dx = (_tmax - _tmin_rates)/1000.; _mfac_rates = 1./_dx;
   for (_i=0, _x=_tmin_rates; _i < 1001; _x += _dx, _i++) {
    _f_rates(_p, _ppvar, _thread, _nt, _x);
    _t_mtau[_i] = mtau;
    _t_ntau[_i] = ntau;
    _t_htau[_i] = htau;
    _t_minf[_i] = minf;
    _t_ninf[_i] = ninf;
    _t_hinf[_i] = hinf;
   }
  }
 }

 static int rates(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _lv) { 
#if 0
_check_rates(_p, _ppvar, _thread, _nt);
#endif
 _n_rates(_p, _ppvar, _thread, _nt, _lv);
 return 0;
 }

 static void _n_rates(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_rates(_p, _ppvar, _thread, _nt, _lv); return; 
}
 _xi = _mfac_rates * (_lv - _tmin_rates);
 if (isnan(_xi)) {
  mtau = _xi;
  ntau = _xi;
  htau = _xi;
  minf = _xi;
  ninf = _xi;
  hinf = _xi;
  return;
 }
 if (_xi <= 0.) {
 mtau = _t_mtau[0];
 ntau = _t_ntau[0];
 htau = _t_htau[0];
 minf = _t_minf[0];
 ninf = _t_ninf[0];
 hinf = _t_hinf[0];
 return; }
 if (_xi >= 1000.) {
 mtau = _t_mtau[1000];
 ntau = _t_ntau[1000];
 htau = _t_htau[1000];
 minf = _t_minf[1000];
 ninf = _t_ninf[1000];
 hinf = _t_hinf[1000];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 mtau = _t_mtau[_i] + _theta*(_t_mtau[_i+1] - _t_mtau[_i]);
 ntau = _t_ntau[_i] + _theta*(_t_ntau[_i+1] - _t_ntau[_i]);
 htau = _t_htau[_i] + _theta*(_t_htau[_i+1] - _t_htau[_i]);
 minf = _t_minf[_i] + _theta*(_t_minf[_i+1] - _t_minf[_i]);
 ninf = _t_ninf[_i] + _theta*(_t_ninf[_i+1] - _t_ninf[_i]);
 hinf = _t_hinf[_i] + _theta*(_t_hinf[_i+1] - _t_hinf[_i]);
 }

 
static int  _f_rates ( _threadargsprotocomma_ double _lv ) {
   double _lalpha , _lbeta , _lsum , _lvt , _lQ ;
 _lvt = _lv + v_shft ;
   _lQ = pow( 3.0 , ( ( 35.0 - 32.0 ) / 10.0 ) ) ;
   if ( _lvt  == 13.1 ) {
     _lalpha = 0.32 * 4.0 ;
     }
   else {
     _lalpha = 0.32 * ( 13.1 - _lvt ) / ( exp ( ( 13.1 - _lvt ) / 4.0 ) - 1.0 ) ;
     }
   if ( _lvt  == 40.1 ) {
     _lbeta = 0.28 * 5.0 ;
     }
   else {
     _lbeta = 0.28 * ( _lvt - 40.1 ) / ( exp ( ( _lvt - 40.1 ) / 5.0 ) - 1.0 ) ;
     }
   _lsum = _lalpha + _lbeta ;
   mtau = 1.0 / _lsum ;
   mtau = mtau / _lQ ;
   minf = _lalpha / _lsum ;
   _lalpha = 0.128 * exp ( ( 17.0 - _lvt ) / 18.0 ) ;
   _lbeta = 4.0 / ( 1.0 + exp ( ( 40.0 - _lvt ) / 5.0 ) ) ;
   _lsum = _lalpha + _lbeta ;
   htau = 1.0 / _lsum ;
   htau = htau / _lQ ;
   hinf = _lalpha / _lsum ;
   if ( _lvt  == 35.1 ) {
     _lalpha = 0.016 * 5.0 ;
     }
   else {
     _lalpha = 0.016 * ( 35.1 - _lvt ) / ( exp ( ( 35.1 - _lvt ) / 5.0 ) - 1.0 ) ;
     }
   _lbeta = 0.25 * exp ( ( 20.0 - _lvt ) / 40.0 ) ;
   _lsum = _lalpha + _lbeta ;
   ntau = 1.0 / _lsum ;
   ntau = ntau / _lQ ;
   ninf = _lalpha / _lsum ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 
#if 1
 _check_rates(_p, _ppvar, _thread, _nt);
#endif
 _r = 1.;
 rates ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 3;}
 
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
	for (_i=0; _i < 3; ++_i) {
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
 
static void _thread_mem_init(Datum* _thread) {
  if (_thread1data_inuse) {_thread[_gth]._pval = (double*)ecalloc(1, sizeof(double));
 }else{
 _thread[_gth]._pval = _thread1data; _thread1data_inuse = 1;
 }
 }
 
static void _thread_cleanup(Datum* _thread) {
  if (_thread[_gth]._pval == _thread1data) {
   _thread1data_inuse = 0;
  }else{
   free((void*)_thread[_gth]._pval);
  }
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  a = a0;
  b = b0;
  h = h0;
  m = m0;
  n = n0;
 {
   rates ( _threadargscomma_ v ) ;
   m = minf ;
   h = hinf ;
   n = ninf ;
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

#if 0
 _check_rates(_p, _ppvar, _thread, _nt);
#endif
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
   gNa = gNabar * h * m * m ;
   iNa = gNa * ( v - eNa ) ;
   gK = gKbar * n ;
   iK = gK * ( v - eK ) ;
   iL = gLbar * ( v - eL ) ;
   i = iL + iK + iNa ;
   totG = gNa + gK + gLbar ;
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
 {   states(_p, _ppvar, _thread, _nt);
  }}}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(m) - _p;  _dlist1[0] = &(Dm) - _p;
 _slist1[1] = &(h) - _p;  _dlist1[1] = &(Dh) - _p;
 _slist1[2] = &(n) - _p;  _dlist1[2] = &(Dn) - _p;
   _t_mtau = makevector(1001*sizeof(double));
   _t_ntau = makevector(1001*sizeof(double));
   _t_htau = makevector(1001*sizeof(double));
   _t_minf = makevector(1001*sizeof(double));
   _t_ninf = makevector(1001*sizeof(double));
   _t_hinf = makevector(1001*sizeof(double));
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "Traub.mod";
static const char* nmodl_file_text = 
  "\n"
  "COMMENT\n"
  "	All the channels are taken from same good old classic articles.\n"
  "	The arrengment was done after:\n"
  "	Kang, S., Kitano, K., and Fukai, T. (2004). \n"
  "		Self-organized two-state membrane potential \n"
  "		transitions in a network of realistically modeled \n"
  "		cortical neurons. Neural Netw 17, 307-312.\n"
  "	\n"
  "	Whenever available I used the same parameters they used,\n"
  "	except in n gate:\n"
  "		n' = phi*(ninf-n)/ntau\n"
  "	Kang used phi = 12\n"
  "	I used phi = 1\n"
  "	\n"
  "	Written by Albert Gidon & Leora Menhaim (2004).\n"
  "ENDCOMMENT\n"
  "\n"
  "UNITS {\n"
  " 	(mA) = (milliamp)\n"
  " 	(mV) = (millivolt)\n"
  "	(S) = (siemens)		\n"
  "}\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX traub\n"
  "	NONSPECIFIC_CURRENT i\n"
  "	RANGE iL,iNa,iK\n"
  "	RANGE eL, eNa, eK\n"
  "	RANGE gLbar, gNabar, gKbar\n"
  "	RANGE v_shft\n"
  " }\n"
  "	\n"
  "	\n"
  "PARAMETER {\n"
  "        gNabar = .03 (S/cm2)	:Traub et. al. 1991\n"
  "        gKbar = .015 (S/cm2) 	:Traub et. al. 1991\n"
  "        gLbar = 0.00014 (S/cm2) :Siu Kang - by email.\n"
  "        eL = -62.0 (mV) :Siu Kang - by email.\n"
  "        eK = -80 (mV)	:Siu Kang - by email.\n"
  "        eNa = 90 (mV)	:Leora\n"
  "        totG = 0\n"
  "		v_shft = 49.2 : shift to apply to all curves\n"
  "}\n"
  " \n"
  "STATE {\n"
  "        m h n a b\n"
  "}\n"
  " \n"
  "ASSIGNED {\n"
  "        v (mV)\n"
  "        i (mA/cm2)\n"
  "        cm (uF)\n"
  "        iL iNa iK(mA/cm2)\n"
  "        gNa gK (S/cm2)\n"
  "	    minf hinf ninf \n"
  "		mtau (ms) htau (ms) ntau (ms) \n"
  "}\n"
  "\n"
  "\n"
  "BREAKPOINT {\n"
  "        SOLVE states METHOD cnexp \n"
  "        :-------------------------\n"
  "        :Traub et. al. 1991\n"
  "        gNa = gNabar*h*m*m\n"
  "		iNa = gNa*(v - eNa)\n"
  "		\n"
  "        gK = gKbar*n : - Traub et. al. 1991\n"
  "		iK = gK*(v - eK)\n"
  "        :-------------------------\n"
  "		iL = gLbar*(v - eL) \n"
  "		i = iL + iK + iNa\n"
  "		:to calculate the input resistance get the sum of\n"
  "		:	all the conductance.\n"
  "		totG = gNa + gK + gLbar\n"
  "			\n"
  "}\n"
  " \n"
  "\n"
  "INITIAL {\n"
  "	rates(v)\n"
  "	m = minf\n"
  "	h = hinf\n"
  "	n = ninf\n"
  "}\n"
  "\n"
  "? states\n"
  "DERIVATIVE states {  \n"
  "	rates(v)\n"
  "	:Traub Spiking channels\n"
  "	m' = (minf-m)/mtau\n"
  "	h' = (hinf-h)/htau\n"
  "	:n' = 2*(ninf-n)/ntau :phi=12 from Kang et. al. 2004\n"
  "	n' = (ninf-n)/ntau :phi=12 from Kang et. al. 2004\n"
  "}\n"
  "\n"
  "? rates\n"
  "DEFINE Q10 3\n"
  "PROCEDURE rates(v(mV)) {  \n"
  "	:Computes rate and other constants at current v.\n"
  "	:Call once from HOC to initialize inf at resting v.\n"
  "	LOCAL  alpha, beta, sum, vt, Q\n"
  "	TABLE 	mtau,ntau,htau,minf,ninf,hinf\n"
  "	FROM -100 TO 70 WITH 1000\n"
  "	: see Resources/The unreliable Q10.htm for details\n"
  "	: remember that not only Q10 is temprature dependent \n"
  "	: and just astimated here, but also the calculation of\n"
  "	: Q is itself acurate only in about 10% in this range of\n"
  "	: temperatures. the transformation formulation is:\n"
  "	: Q = Q10^(( new(degC) - from_original_experiment(degC) )/ 10)\n"
  "	\n"
  "		:--------------------------------------------------------\n"
  "		\n"
  "		: This part was taken **directly** from:\n"
  "		: Traub, R. D., Wong, R. K., Miles, R., and Michelson, H. (1991). \n"
  "		:	A model of a CA3 hippocampal pyramidal neuron incorporating \n"
  "		:	voltage-clamp data on intrinsic conductances. \n"
  "		:	J Neurophysiol 66, 635-650.\n"
  "		:	Experiments were done in >=32degC for m,h\n"
  "		: Traub et al uses their -60mV as 0mV thus here is the shift\n"
  "		vt = v + v_shft :49.2\n"
  "		Q = Q10^((35 - 32)/ 10)\n"
  "		:\"m\" sodium activation system\n"
  "		if(vt == 13.1){alpha = 0.32*4}\n"
  "		else{alpha = 0.32*(13.1 - vt)/(exp((13.1 - vt)/4) - 1)}\n"
  "		if(vt == 40.1){beta = 0.28*5}\n"
  "		else{beta = 0.28*(vt - 40.1)/(exp((vt - 40.1)/5)-1)}\n"
  "        sum = alpha + beta\n"
  "		mtau = 1/sum\n"
  "		mtau = mtau/Q\n"
  "        minf = alpha/sum\n"
  "\n"
  "       :\"h\" sodium inactivation system\n"
  "		alpha = 0.128*exp((17 - vt)/18)\n"
  "		beta = 4/(1 + exp((40 - vt)/5))\n"
  "        sum = alpha + beta\n"
  "		htau = 1/sum\n"
  "		htau = htau/Q\n"
  "        hinf = alpha/sum\n"
  "\n"
  "    	:\"n\" potassium activation system\n"
  "    	if(vt == 35.1){ alpha = 0.016*5 }\n"
  "		else{alpha =0.016*(35.1 - vt)/(exp((35.1 - vt)/5) - 1)}\n"
  "		beta = 0.25*exp((20 - vt)/40)\n"
  "		sum = alpha + beta\n"
  "        ntau = 1/sum\n"
  "        ntau = ntau/Q\n"
  "        ninf = alpha/sum\n"
  "}\n"
  "\n"
  "\n"
  ;
#endif
