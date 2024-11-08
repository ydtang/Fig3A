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
 
#define nrn_init _nrn_init__namr
#define _nrn_initial _nrn_initial__namr
#define nrn_cur _nrn_cur__namr
#define _nrn_current _nrn_current__namr
#define nrn_jacob _nrn_jacob__namr
#define nrn_state _nrn_state__namr
#define _net_receive _net_receive__namr 
#define rates rates__namr 
#define states states__namr 
 
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
#define gnabar _p[0]
#define vhalf _p[1]
#define vhalfn _p[2]
#define vhalfl _p[3]
#define rmax _p[4]
#define vvh _p[5]
#define b _p[6]
#define Ir _p[7]
#define n _p[8]
#define l _p[9]
#define r _p[10]
#define fr _p[11]
#define ena _p[12]
#define Dn _p[13]
#define Dl _p[14]
#define Dr _p[15]
#define Dfr _p[16]
#define ina _p[17]
#define _g _p[18]
#define _ion_ena	*_ppvar[0]._pval
#define _ion_ina	*_ppvar[1]._pval
#define _ion_dinadv	*_ppvar[2]._pval
 
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
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_alpl(void);
 static void _hoc_alpr(void);
 static void _hoc_alpv(void);
 static void _hoc_alpn(void);
 static void _hoc_alp(void);
 static void _hoc_betl(void);
 static void _hoc_betr(void);
 static void _hoc_betn(void);
 static void _hoc_rates(void);
 static void _hoc_states(void);
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
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_namr", _hoc_setdata,
 "alpl_namr", _hoc_alpl,
 "alpr_namr", _hoc_alpr,
 "alpv_namr", _hoc_alpv,
 "alpn_namr", _hoc_alpn,
 "alp_namr", _hoc_alp,
 "betl_namr", _hoc_betl,
 "betr_namr", _hoc_betr,
 "betn_namr", _hoc_betn,
 "rates_namr", _hoc_rates,
 "states_namr", _hoc_states,
 0, 0
};
#define alpl alpl_namr
#define alpr alpr_namr
#define alpv alpv_namr
#define alpn alpn_namr
#define alp alp_namr
#define betl betl_namr
#define betr betr_namr
#define betn betn_namr
 extern double alpl( double , double );
 extern double alpr( double );
 extern double alpv( double , double );
 extern double alpn( double , double );
 extern double alp( double , double );
 extern double betl( double , double );
 extern double betr( double );
 extern double betn( double , double );
 /* declare global and static user variables */
#define a0r a0r_namr
 double a0r = 0.0003;
#define a0n a0n_namr
 double a0n = 2;
#define a0l a0l_namr
 double a0l = 0.1;
#define b0r b0r_namr
 double b0r = 0.0003;
#define b0n b0n_namr
 double b0n = 2;
#define b0l b0l_namr
 double b0l = 0.1;
#define gml gml_namr
 double gml = 0.65;
#define gmr gmr_namr
 double gmr = 0.2;
#define gmn gmn_namr
 double gmn = 0.9;
#define inf inf_namr
 double inf = 0;
#define lmax lmax_namr
 double lmax = 1;
#define linf linf_namr
 double linf = 0;
#define nmax nmax_namr
 double nmax = 0.01;
#define ninf ninf_namr
 double ninf = 0;
#define rinf rinf_namr
 double rinf = 0;
#define taur taur_namr
 double taur = 0;
#define taun taun_namr
 double taun = 0;
#define taul taul_namr
 double taul = 0;
#define vvs vvs_namr
 double vvs = 2;
#define vhalfr vhalfr_namr
 double vhalfr = -60;
#define zetal zetal_namr
 double zetal = 4;
#define zetar zetar_namr
 double zetar = 12;
#define zetan zetan_namr
 double zetan = -4;
#define zeta zeta_namr
 double zeta = 20;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "zeta_namr", "1",
 "vhalfr_namr", "mV",
 "a0l_namr", "/ms",
 "b0l_namr", "/ms",
 "a0n_namr", "/ms",
 "b0n_namr", "/ms",
 "a0r_namr", "ms",
 "b0r_namr", "ms",
 "zetan_namr", "1",
 "zetar_namr", "1",
 "zetal_namr", "1",
 "gmn_namr", "1",
 "gmr_namr", "1",
 "gml_namr", "1",
 "lmax_namr", "1",
 "nmax_namr", "1",
 "vvs_namr", "1",
 "gnabar_namr", "mho/cm2",
 "vhalf_namr", "mV",
 "vhalfn_namr", "mV",
 "vhalfl_namr", "mV",
 "rmax_namr", "1",
 "vvh_namr", "mv",
 "Ir_namr", "mA/cm2",
 0,0
};
 static double delta_t = 0.01;
 static double fr0 = 0;
 static double l0 = 0;
 static double n0 = 0;
 static double r0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "zeta_namr", &zeta_namr,
 "vhalfr_namr", &vhalfr_namr,
 "a0l_namr", &a0l_namr,
 "b0l_namr", &b0l_namr,
 "a0n_namr", &a0n_namr,
 "b0n_namr", &b0n_namr,
 "a0r_namr", &a0r_namr,
 "b0r_namr", &b0r_namr,
 "zetan_namr", &zetan_namr,
 "zetar_namr", &zetar_namr,
 "zetal_namr", &zetal_namr,
 "gmn_namr", &gmn_namr,
 "gmr_namr", &gmr_namr,
 "gml_namr", &gml_namr,
 "lmax_namr", &lmax_namr,
 "nmax_namr", &nmax_namr,
 "vvs_namr", &vvs_namr,
 "inf_namr", &inf_namr,
 "ninf_namr", &ninf_namr,
 "linf_namr", &linf_namr,
 "rinf_namr", &rinf_namr,
 "taul_namr", &taul_namr,
 "taun_namr", &taun_namr,
 "taur_namr", &taur_namr,
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
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"namr",
 "gnabar_namr",
 "vhalf_namr",
 "vhalfn_namr",
 "vhalfl_namr",
 "rmax_namr",
 "vvh_namr",
 "b_namr",
 0,
 "Ir_namr",
 0,
 "n_namr",
 "l_namr",
 "r_namr",
 "fr_namr",
 0,
 0};
 static Symbol* _na_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 19, _prop);
 	/*initialize range parameters*/
 	gnabar = 0.01;
 	vhalf = -50;
 	vhalfn = -30;
 	vhalfl = -57;
 	rmax = 3;
 	vvh = -58;
 	b = 0;
 	_prop->param = _p;
 	_prop->param_size = 19;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_na_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ena */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ina */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dinadv */
 
}
 static void _initlists();
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _namir_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("na", -10000.);
 	_na_sym = hoc_lookup("na_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 19, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "na_ion");
 	hoc_register_cvode(_mechtype, _ode_count, 0, 0, 0);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 namr D:/Dendritic action potentials and computation in/GidonEtAl2020_fig3andS9/NEURON_Code/_mod/namir.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double _zfacn , _zfacl , _zfacr ;
static int _reset;
static char *modelname = "Na channel";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates(double, double, double, double, double, double);
static int states();
 
double alp (  double _lv , double _lvf ) {
   double _lalp;
 _lalp = exp ( 1.e-3 * zeta * ( _lv - _lvf ) * 9.648e4 / ( 8.315 * ( 273.16 + celsius ) ) ) ;
   
return _lalp;
 }
 
static void _hoc_alp(void) {
  double _r;
   _r =  alp (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double alpn (  double _lv , double _lvn ) {
   double _lalpn;
 _lalpn = exp ( 1.e-3 * zetan * ( _lv - _lvn ) * 9.648e4 / ( 8.315 * ( 273.16 + celsius ) ) ) ;
   
return _lalpn;
 }
 
static void _hoc_alpn(void) {
  double _r;
   _r =  alpn (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double betn (  double _lv , double _lvn ) {
   double _lbetn;
 _lbetn = exp ( 1.e-3 * zetan * gmn * ( _lv - _lvn ) * 9.648e4 / ( 8.315 * ( 273.16 + celsius ) ) ) ;
   
return _lbetn;
 }
 
static void _hoc_betn(void) {
  double _r;
   _r =  betn (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double alpv (  double _lv , double _lvh ) {
   double _lalpv;
 _lalpv = ( 1.0 + b * exp ( ( _lv - _lvh ) / vvs ) ) / ( 1.0 + exp ( ( _lv - _lvh ) / vvs ) ) ;
   
return _lalpv;
 }
 
static void _hoc_alpv(void) {
  double _r;
   _r =  alpv (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double alpr (  double _lv ) {
   double _lalpr;
 _lalpr = exp ( 1.e-3 * zetar * ( _lv - vhalfr ) * 9.648e4 / ( 8.315 * ( 273.16 + celsius ) ) ) ;
   
return _lalpr;
 }
 
static void _hoc_alpr(void) {
  double _r;
   _r =  alpr (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double betr (  double _lv ) {
   double _lbetr;
 _lbetr = exp ( 1.e-3 * zetar * gmr * ( _lv - vhalfr ) * 9.648e4 / ( 8.315 * ( 273.16 + celsius ) ) ) ;
   
return _lbetr;
 }
 
static void _hoc_betr(void) {
  double _r;
   _r =  betr (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double alpl (  double _lv , double _lvl ) {
   double _lalpl;
 _lalpl = exp ( 1.e-3 * zetal * ( _lv - _lvl ) * 9.648e4 / ( 8.315 * ( 273.16 + celsius ) ) ) ;
   
return _lalpl;
 }
 
static void _hoc_alpl(void) {
  double _r;
   _r =  alpl (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double betl (  double _lv , double _lvl ) {
   double _lbetl;
 _lbetl = exp ( 1.e-3 * zetal * gml * ( _lv - _lvl ) * 9.648e4 / ( 8.315 * ( 273.16 + celsius ) ) ) ;
   
return _lbetl;
 }
 
static void _hoc_betl(void) {
  double _r;
   _r =  betl (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
static int  states (  ) {
   rates ( _threadargscomma_ v , vhalf , vhalfn , vhalfl , vvh , rmax ) ;
   n = n + _zfacn * ( ninf - n ) ;
   l = l + _zfacl * ( linf - l ) ;
   r = r + _zfacr * ( rinf - r ) ;
   fr = inf ;
   
/*VERBATIM*/
        return 0;
  return 0; }
 
static void _hoc_states(void) {
  double _r;
   _r = 1.;
 states (  );
 hoc_retpushx(_r);
}
 
static int  rates (  double _lv , double _lvf , double _lvn , double _lvl , double _lvh , double _lrm ) {
   double _la , _lq10 ;
 _lq10 = pow( 1.5 , ( ( celsius - 22.0 ) / 10.0 ) ) ;
   inf = 1.0 / ( 1.0 + alp ( _threadargscomma_ _lv , _lvf ) ) ;
   _la = alpn ( _threadargscomma_ _lv , _lvn ) ;
   ninf = 1.0 / ( 1.0 + _la ) ;
   taun = betn ( _threadargscomma_ _lv , _lvn ) / ( _lq10 * ( a0n + b0n * _la ) ) ;
   if ( taun < nmax ) {
     taun = nmax ;
     }
   _zfacn = ( 1.0 - exp ( - dt / taun ) ) ;
   _la = alpl ( _threadargscomma_ _lv , _lvl ) ;
   linf = 1.0 / ( 1.0 + _la ) ;
   taul = betl ( _threadargscomma_ _lv , _lvl ) / ( _lq10 * ( a0l + b0l * _la ) ) ;
   if ( taul < lmax ) {
     taul = lmax ;
     }
   _zfacl = ( 1.0 - exp ( - dt / taul ) ) ;
   rinf = alpv ( _threadargscomma_ _lv , _lvh ) ;
   taur = betr ( _threadargscomma_ _lv ) / ( a0r + b0r * alpr ( _threadargscomma_ _lv ) ) ;
   if ( taur < _lrm ) {
     taur = _lrm ;
     }
   _zfacr = ( 1.0 - exp ( - dt / taur ) ) ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   _r = 1.;
 rates (  *getarg(1) , *getarg(2) , *getarg(3) , *getarg(4) , *getarg(5) , *getarg(6) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ hoc_execerror("namr", "cannot be used with CVODE"); return 0;}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_na_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_na_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_na_sym, _ppvar, 2, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  fr = fr0;
  l = l0;
  n = n0;
  r = r0;
 {
   rates ( _threadargscomma_ v , vhalf , vhalfn , vhalfl , vvh , rmax ) ;
   n = ninf ;
   l = linf ;
   r = rinf ;
   fr = inf ;
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
  ena = _ion_ena;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   ina = gnabar * n * l * r * ( v - ena ) ;
   Ir = - ina * fr ;
   }
 _current += ina;
 _current += Ir;

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
  ena = _ion_ena;
 _g = _nrn_current(_v + .001);
 	{ double _dina;
  _dina = ina;
 _rhs = _nrn_current(_v);
  _ion_dinadv += (_dina - ina)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ina += ina ;
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
  ena = _ion_ena;
 { error =  states();
 if(error){fprintf(stderr,"at line 80 in file namir.mod:\n	SOLVE states\n"); nrn_complain(_p); abort_run(error);}
 } }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "namir.mod";
static const char* nmodl_file_text = 
  "TITLE  Na channel\n"
  ": with recovery from inactivation and ir. M.Migliore BJ 1996\n"
  "NEURON {\n"
  "	SUFFIX namr\n"
  "	USEION na READ ena WRITE ina\n"
  "	NONSPECIFIC_CURRENT Ir\n"
  "        RANGE  gnabar,vhalf,vhalfn,vhalfl,vvh,rmax, r, b\n"
  "        GLOBAL ninf,linf,taul,taun,rinf,taur,inf\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(mA) = (milliamp)\n"
  "	(mV) = (millivolt)\n"
  "\n"
  "}\n"
  "\n"
  "INDEPENDENT {t FROM 0 TO 1 WITH 100 (ms)}\n"
  "\n"
  "PARAMETER {\n"
  "        dt (ms)\n"
  "	v (mV)\n"
  "        ena=50 (mV)\n"
  "	celsius = 22	(degC)\n"
  "	gnabar=.01 (mho/cm2)\n"
  "	vhalf=-50    (mV)\n"
  "	zeta=20       (1)\n"
  "        vhalfn=-30   (mV)\n"
  "        vhalfr=-60   (mV)\n"
  "        vhalfl=-57   (mV)\n"
  "        a0l=0.1      (/ms)\n"
  "        b0l=0.1      (/ms)\n"
  "        a0n=2    (/ms)\n"
  "        b0n=2    (/ms)\n"
  "        a0r=0.0003    (ms)\n"
  "        b0r=0.0003    (ms)\n"
  "        zetan=-4    (1)\n"
  "        zetar=12    (1)\n"
  "        zetal=4    (1)\n"
  "        gmn=0.9   (1)\n"
  "        gmr=0.2   (1)\n"
  "        gml=0.65   (1)\n"
  "	lmax=1   (1)\n"
  "	nmax=0.01   (1)\n"
  "	rmax=3 (1)\n"
  "	vvh=-58     (mv)\n"
  "	vvs=2 	(1)\n"
  "	b=0\n"
  "}\n"
  "\n"
  "\n"
  "STATE {\n"
  "	n\n"
  "        l\n"
  "	r\n"
  "	fr\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	Ir  (mA/cm2)\n"
  "	ina (mA/cm2)\n"
  "	inf\n"
  "        ninf\n"
  "        linf      \n"
  "	rinf\n"
  "        taul\n"
  "        taun\n"
  "	taur\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	rates(v,vhalf,vhalfn,vhalfl,vvh,rmax)\n"
  "	n=ninf\n"
  "	l=linf\n"
  "	r=rinf\n"
  "	fr=inf\n"
  "}\n"
  "\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE states\n"
  "	ina = gnabar*n*l*r*(v-ena)\n"
  "	Ir = -ina*fr\n"
  ":	Ir = 0\n"
  "}\n"
  "\n"
  "FUNCTION alp(v(mV),vf) { \n"
  "  alp = exp( 1.e-3*zeta*(v-vf)*9.648e4/(8.315*(273.16+celsius)))\n"
  "}       \n"
  "\n"
  "FUNCTION alpn(v(mV),vn) {\n"
  "  alpn = exp(1.e-3*zetan*(v-vn)*9.648e4/(8.315*(273.16+celsius))) \n"
  "}\n"
  "\n"
  "FUNCTION betn(v(mV),vn) {\n"
  "  betn = exp(1.e-3*zetan*gmn*(v-vn)*9.648e4/(8.315*(273.16+celsius))) \n"
  "}\n"
  "\n"
  "FUNCTION alpv(v(mV),vh) {\n"
  "         alpv = (1+b*exp((v-vh)/vvs))/(1+exp((v-vh)/vvs))\n"
  "}\n"
  "\n"
  "FUNCTION alpr(v(mV)) {\n"
  "  alpr = exp(1.e-3*zetar*(v-vhalfr)*9.648e4/(8.315*(273.16+celsius))) \n"
  "}\n"
  "\n"
  "FUNCTION betr(v(mV)) {\n"
  "  betr = exp(1.e-3*zetar*gmr*(v-vhalfr)*9.648e4/(8.315*(273.16+celsius))) \n"
  "}\n"
  "\n"
  "FUNCTION alpl(v(mV),vl) {\n"
  "  alpl = exp(1.e-3*zetal*(v-vl)*9.648e4/(8.315*(273.16+celsius))) \n"
  "}\n"
  "\n"
  "FUNCTION betl(v(mV),vl) {\n"
  "  betl = exp(1.e-3*zetal*gml*(v-vl)*9.648e4/(8.315*(273.16+celsius))) \n"
  "}\n"
  "\n"
  "LOCAL facn,facl,facr\n"
  "\n"
  ":if state_borgka is called from hoc, garbage or segmentation violation will\n"
  ":result because range variables won't have correct pointer.  This is because\n"
  ": only BREAKPOINT sets up the correct pointers to range variables.\n"
  "PROCEDURE states() {     : exact when v held constant; integrates over dt step\n"
  "        rates(v,vhalf,vhalfn,vhalfl,vvh,rmax)\n"
  "        n = n + facn*(ninf - n)\n"
  "        l = l + facl*(linf - l)\n"
  "        r = r + facr*(rinf - r)\n"
  "	fr = inf\n"
  "        VERBATIM\n"
  "        return 0;\n"
  "        ENDVERBATIM\n"
  "}\n"
  "\n"
  "PROCEDURE rates(v (mV),vf,vn,vl,vh,rm) { :callable from hoc\n"
  "        LOCAL a,q10\n"
  "        q10=1.5^((celsius-22)/10)\n"
  "        inf = 1/(1 + alp(v,vf))\n"
  "        a = alpn(v,vn)\n"
  "        ninf = 1/(1 + a)\n"
  "        taun = betn(v,vn)/(q10*(a0n+b0n*a))\n"
  "	if (taun<nmax) {taun=nmax}\n"
  "        facn = (1 - exp(-dt/taun))\n"
  "        a = alpl(v,vl)\n"
  "        linf = 1/(1+ a)\n"
  "        taul = betl(v,vl)/(q10*(a0l + b0l*a))\n"
  "	if (taul<lmax) {taul=lmax}\n"
  "        facl = (1 - exp(-dt/taul))\n"
  "        rinf = alpv(v,vh)\n"
  "        taur = betr(v)/(a0r+b0r*alpr(v))\n"
  "	if (taur<rm) {taur=rm}\n"
  "        facr = (1 - exp(-dt/taur))\n"
  "}\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  ;
#endif
