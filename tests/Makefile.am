# this file is part of gnucap
# (c) 2013 Felix Salfelder
# License: GPLv3
#
GNUCAP = $(top_builddir)/src/gnucap

# needed for gnucap-run.sh
AM_TESTS_ENVIRONMENT = \
    export PATH='../src:$(PATH)'\
           GNUCAP='gnucap'\
           REDIRECT='exec 2>&9'\
           MAKEFLAGS=--no-print-directory \
           srcdir='$(srcdir)';
AM_TESTS_FD_REDIRECT = 9>&2
TEST_EXTENSIONS = .gc .ckt

GC_LOG_COMPILER = $(top_srcdir)/tests/gc_log_compiler
CKT_LOG_COMPILER = $(top_srcdir)/tests/ckt_log_compiler
SH_LOG_COMPILER = $(SHELL)

AM_CKT_LOG_FLAGS = $(srcdir)

TESTS = \
    gnucap_run.sh \
    $(GC_TESTS) \
    $(LOGIC_TESTS) \
    $(UF_TESTS) \
    $(CKT_TESTS)

REF = ==out

patchlev_fake_%.h: $(top_builddir)/include/patchlev.h
	sed 's/e VERSION_REVISION.*/e VERSION_REVISION $*/' < $< > $@
patchlev_other.h: $(top_builddir)/include/patchlev.h
	sed 's/e VERSION_CURRENT .*/e VERSION_CURRENT 99999/' < $< > $@

BUILT_SOURCES = patchlev_fake_999.h patchlev_fake_0.h patchlev_other.h

transform=s&\.real$$&&;$(program_transform_name)
plugindir_transformed = ${libdir}/$(shell echo gnucap | sed -e "$(transform)")
plugindir:= $(plugindir_transformed)
plugin_PROGRAMS = gnucap-testfilter
gnucap_testfilter_SOURCES = filter.c
empty=
$(empty)install-pluginPROGRAMS: transform=

CKT_TESTS = \
    1lin.ckt \
    ap_ctof.2.ckt \
    ap_ctof.ckt \
    basic.ckt \
    bm_complex.1a.ckt \
    bm_complex.1.ckt \
    bm_complex.2.ckt \
    bm_complex.3.ckt \
    bm_complex.4.ckt \
    bm_cond.1.ckt \
    bm_cond.2.ckt \
    bm_cond.3.ckt \
    bm_exp.1.ckt \
    bm_exp.1p.ckt \
    bm_exp.2.ckt \
    bm_exp.3.ckt \
    bm_exp.3p.ckt \
    bm_exp.4p.ckt \
    bm_exp.5p.ckt \
    bm_fit.0.1.ckt \
    bm_fit.1.1.ckt \
    bm_fit.2.1.ckt \
    bm_fit.2.2.ckt \
    bm_fit.2.3.ckt \
    bm_fit.2.4.ckt \
    bm_fit.3.1.ckt \
    bm_fit.3.2.ckt \
    bm_fit.3.3.ckt \
    bmm_cap.1.ckt \
    bmm_cap.2.ckt \
    bmm_cap.3.ckt \
    bmm_cap.4.ckt \
    bmm_cap.5.ckt \
    bmm_res.1.ckt \
    bmm_res.1p.ckt \
    bmm_res.2.ckt \
    bmm_res.2p.ckt \
    bmm_res.3.ckt \
    bmm_res.3p.ckt \
    bmm_res.4.ckt \
    bmm_res.4p.ckt \
    bmm_res.5.ckt \
    bmm_res.5p.ckt \
    bm_poly.0.ckt \
    bm_poly.1.ckt \
    bm_poly.2.ckt \
    bm_posy.1.ckt \
    bm_pulse.1a.ckt \
    bm_pulse.1.ckt \
    bm_pulse.2a.ckt \
    bm_pulse.2.ckt \
    bm_pulse.3a.ckt \
    bm_pulse.3.ckt \
    bm_pulse.4.ckt \
    bm_pwl.1.ckt \
    bm_pwl.2.ckt \
    bm_pwl.3.ckt \
    bm_pwl.4.ckt \
    bm_pwl.5.ckt \
    bm_sffm.1.ckt \
    bm_sffm.2.ckt \
    bm_sffm.3.ckt \
    bm_sffm.4.ckt \
    bm_sffm.5.ckt \
    bm_sffm.6.ckt \
    bm_sin.1.ckt \
    bm_sin.2.ckt \
    bm_sin.3.ckt \
    bm_sin.4.ckt \
    bm_sin.5.ckt \
    bm_sin.6.ckt \
    bm_sin.7.ckt \
    bm_sin.8.ckt \
    bm_sin.9.ckt \
    bm_table.1.1.ckt \
    bm_table.2.1.ckt \
    bm_table.2.2.ckt \
    bm_table.2.3.ckt \
    bm_table.2.4.ckt \
    bm_table.3.1.ckt \
    bm_table.3.2.ckt \
    bm_table.3.3.ckt \
    bm_table.error1.ckt \
    bm_tanh.1.ckt \
    bm_tanh.2.ckt \
    cc.ckt \
    c.ckt \
    c_delete.2.ckt \
    c_delete.3.ckt \
    c_delete.4.ckt \
    c_delete.5.ckt \
    c_delete.6.ckt \
    c_delete.7.ckt \
    c_delete.ckt \
    c_eval.1.ckt \
    c_getckt.ckt \
    charge-cons2.ckt \
    charge-cons3.ckt \
    charge-cons.ckt \
    c_help.1.ckt \
    c_lib.1.ckt \
    c_list.1.ckt \
    c_list.ckt \
    c_measure.1.ckt \
    c_measure.2.ckt \
    c_modify.ckt \
    cont.1.ckt \
    cont.2.ckt \
    cont.3.ckt \
    current.ckt \
    d_bjt.2.ckt \
    d_bjt.ckt \
    d_bjt-diffpair-ccs.ckt \
    d_bjt-diffpair-cjc.ckt \
    d_bjt-diffpair-dc.ckt \
    d_bjt-diffpair-nocap.ckt \
    d_bjt-diffpair-tf.ckt \
    d_bjt-diffpair-tr2.ckt \
    d_bjt-diffpair-tran.ckt \
    d_bjt-diffpair-tr.ckt \
    d_bjt.error1.ckt \
    d_bjt-pnp.ckt \
    d_bjt-rca3040.ckt \
    d_bjt-rtlinv.ckt \
    d_bjt-schmitt-bypass.ckt \
    d_bjt-schmitt-nobypass.ckt \
    d_cap.1.ckt \
    d_cap.2.ckt \
    d_cap.3a.ckt \
    d_cap.3.ckt \
    d_cap.4.ckt \
    d_cap.5.ckt \
    d_cap.6.ckt \
    d_cap.ic1.ckt \
    d_cap.ic2.ckt \
    d_cap.ic4.ckt \
    d_cccs.1.ckt \
    d_cccs.2.ckt \
    d_cccs.3.ckt \
    d_cccs.4.ckt \
    d_cccs.5.ckt \
    d_cccs.6.ckt \
    d_ccvs.1.ckt \
    d_ccvs.2.ckt \
    d_ccvs.3.ckt \
    d_coil.1a.ckt \
    d_coil.1b.ckt \
    d_coil.1c.ckt \
    d_coil.1.ckt \
    d_coil.1d.ckt \
    d_coil.1e.ckt \
    d_coil.2.ckt \
    d_coil.3.ckt \
    d_coil.4.ckt \
    d_coil.5.ckt \
    d_coil.error1.ckt \
    d_coil.error2.ckt \
    d_coil.error3.ckt \
    d_coil.error4.ckt \
    d_coil.error5.ckt \
    d_coil.m3a.ckt \
    d_coil.m3.ckt \
    d_coil.m3r.ckt \
    dcsweep.ckt \
    dc_temp.ckt \
    d_diode.10a.ckt \
    d_diode.10.ckt \
    d_diode.11a.ckt \
    d_diode.11b.ckt \
    d_diode.11c.ckt \
    d_diode.11.ckt \
    d_diode.11d.ckt \
    d_diode.12.ckt \
    d_diode.1.ckt \
    d_diode.2.ckt \
    d_diode.3.ckt \
    d_diode.4.ckt \
    d_diode.5.ckt \
    d_diode.6.ckt \
    d_diode.7.ckt \
    d_diode.8.ckt \
    d_diode.9.ckt \
    d_diode.error1.ckt \
    d_logic.1a.ckt \
    d_logic.1.ckt \
    d_logic.2.ckt \
    d_logic.3.ckt \
    d_logic.4.ckt \
    d_logic.5.ckt \
    d_logic.6.ckt \
    d_logic.7.ckt \
    d_logic.error1.ckt \
    d_macro.1.ckt \
    d_mos1.bin1.ckt \
    d_mos1.bin2.ckt \
    d_mos1.bin3.ckt \
    d_mos1.inv1.ckt \
    d_mos1.inv2.ckt \
    d_mos1.inv3.ckt \
    d_mos1.lin1.ckt \
    d_mos1.lin1.late.ckt \
    d_mos1.lin1.none.ckt \
    d_mos1.lin2.ckt \
    d_mos1.lin3.ckt \
    d_mos1.lin4.ckt \
    d_mos1.lin5.ckt \
    d_mos1.n0.ckt \
    d_mos1.n-1.ckt \
    d_mos1.n+1.ckt \
    d_mos1.nand1.ckt \
    d_mos1.nand2.ckt \
    d_mos1.nand3.ckt \
    d_mos1.p0.ckt \
    d_mos1.p-1.ckt \
    d_mos1.p+1.ckt \
    d_mos1.plin.ckt \
    d_mos1.psat1.ckt \
    d_mos1.psat2.ckt \
    d_mos1.psat3.ckt \
    d_mos1.sat1.ckt \
    d_mos1.sat2.ckt \
    d_mos1.sat3.ckt \
    d_mos1.sat4.ckt \
    d_mos1.sat5.ckt \
    d_mos1.sat6.ckt \
    d_mos1.sat7.ckt \
    d_mos1.sts.ckt \
    d_mos1.zero.ckt \
    d_mos2.inv1.ckt \
    d_mos2.lin1.ckt \
    d_mos2.lin1.late.ckt \
    d_mos2.lin1r.ckt \
    d_mos2.lin1z.ckt \
    d_mos2.lin2.ckt \
    d_mos2.lin3.ckt \
    d_mos2.lin4.ckt \
    d_mos2.lin5.ckt \
    d_mos2.n0.ckt \
    d_mos2.n-1.ckt \
    d_mos2.n+1.ckt \
    d_mos2.nand1.ckt \
    d_mos2.p0.ckt \
    d_mos2.p-1.ckt \
    d_mos2.p+1.ckt \
    d_mos2.plin.ckt \
    d_mos2.psat1.ckt \
    d_mos2.psat2.ckt \
    d_mos2.psat3.ckt \
    d_mos2.sat1.ckt \
    d_mos2.sat2.ckt \
    d_mos2.sat3.ckt \
    d_mos2.sat4.ckt \
    d_mos2.sat5.ckt \
    d_mos2.sat6.ckt \
    d_mos2.sat7.ckt \
    d_mos2.sts.ckt \
    d_mos2.zero.ckt \
    d_mos3.inv1.ckt \
    d_mos3.lin1.ckt \
    d_mos3.lin1r.ckt \
    d_mos3.lin1z.ckt \
    d_mos3.lin2.ckt \
    d_mos3.lin3.ckt \
    d_mos3.lin4.ckt \
    d_mos3.lin5.ckt \
    d_mos3.n0.ckt \
    d_mos3.n-1.ckt \
    d_mos3.n+1.ckt \
    d_mos3.nand1.ckt \
    d_mos3.p0.ckt \
    d_mos3.p-1.ckt \
    d_mos3.p+1.ckt \
    d_mos3.plin.ckt \
    d_mos3.psat1.ckt \
    d_mos3.psat2.ckt \
    d_mos3.psat3.ckt \
    d_mos3.sat000.ckt \
    d_mos3.sat007.ckt \
    d_mos3.sat1.ckt \
    d_mos3.sat2.ckt \
    d_mos3.sat3.ckt \
    d_mos3.sat4.ckt \
    d_mos3.sat5.ckt \
    d_mos3.sat6.ckt \
    d_mos3.sat7.ckt \
    d_mos3.sts000.ckt \
    d_mos3.sts.ckt \
    d_mos3.zero.ckt \
    d_mos4.1.ckt \
    d_mos4.2.ckt \
    d_mos4.3.ckt \
    d_mos49.nand1.ckt \
    d_mos5.1.ckt \
    d_mos6.inv1.ckt \
    d_mos6.lin1.ckt \
    d_mos6.lin1r.ckt \
    d_mos6.lin1z.ckt \
    d_mos6.lin2.ckt \
    d_mos6.lin2r.ckt \
    d_mos6.lin2z.ckt \
    d_mos6.n0.ckt \
    d_mos6.n-1.ckt \
    d_mos6.n+1.ckt \
    d_mos6.nand1.ckt \
    d_mos6.p0.ckt \
    d_mos6.p-1.ckt \
    d_mos6.p+1.ckt \
    d_mos6.plin.ckt \
    d_mos6.psat1.ckt \
    d_mos6.psat2.ckt \
    d_mos6.psat3.ckt \
    d_mos6.sat1.ckt \
    d_mos6.sat1r.ckt \
    d_mos6.sat1s.ckt \
    d_mos6.sat1v.ckt \
    d_mos6.sat2.ckt \
    d_mos6.sat2r.ckt \
    d_mos6.sat3.ckt \
    d_mos6.sat3r.ckt \
    d_mos6.sts000.ckt \
    d_mos6.sts.ckt \
    d_mos6.zero.ckt \
    d_mos7.alpha.ckt \
    d_mos7.inv1.ckt \
    d_mos7.lin1.ckt \
    d_mos7.lin2.ckt \
    d_mos7.nand1.ckt \
    d_mos7.nand1.nobypass.ckt \
    d_mos7.plin.ckt \
    d_mos7.psat.ckt \
    d_mos7.sat1.ckt \
    d_mos7.sat2.ckt \
    d_mos7.sat3.ckt \
    d_mos7.sts.ckt \
    d_mos7.zero.ckt \
    d_mos8c2.nand1.ckt \
    d_mos8c2.nand1.nobypass.ckt \
    d_mos8.nand1.ckt \
    d_mos8.sat2.ckt \
    d_mos8.sat3.ckt \
    d_mos8.sat4.ckt \
    d_mos.error1.ckt \
    d_mos.error2.ckt \
    d_mos.error3.ckt \
    d_mos.error4.ckt \
    d_subckt.1.ckt \
    d_subckt.2.ckt \
    d_subckt.3.ckt \
    d_subckt.4.ckt \
    d_subckt.5.ckt \
    d_subckt.6.ckt \
    d_subckt.7.ckt \
    d_subckt.8.ckt \
    d_subckt.9.ckt \
    d_subckt.error1.ckt \
    d_subckt.error2.ckt \
    d_subckt.error3b.ckt \
    d_subckt.error3.ckt \
    d_subckt.error4.ckt \
    d_switch.c0.ckt \
    d_switch.c1a.ckt \
    d_switch.c1.ckt \
    d_switch.c1.inc.ckt \
    d_switch.c1.noinc.ckt \
    d_switch.c1r.ckt \
    d_switch.c2a.ckt \
    d_switch.c2.ckt \
    d_switch.c2r.ckt \
    d_switch.c3.ckt \
    d_switch.error1.ckt \
    d_switch.error2.ckt \
    d_switch.error3.ckt \
    d_switch.nro.1.ckt \
    d_switch.nro.1e.ckt \
    d_switch.nro.1g.ckt \
    d_switch.nro.1gr.ckt \
    d_switch.nro.2.auto.ckt \
    d_switch.nro.2.ckt \
    d_switch.nro.3.ckt \
    d_switch.nro.3h.ckt \
    d_switch.nro.4.ckt \
    d_switch.v.ckt \
    d_switch.v.cont.ckt \
    d_switch.vr.ckt \
    d_tcap.1.ckt \
    d_tcap.2.ckt \
    d_tcap.3.ckt \
    d_tcap.4.ckt \
    d_tcap.5.ckt \
    d_trln.ac.ckt \
    d_trln.tr.100.im.ckt \
    d_trln.tr.100.nim.ckt \
    d_trln.tr.300.im.ckt \
    d_trln.tr.300.nim.ckt \
    d_trln.tr.im.ckt \
    d_trln.tr.m.ckt \
    d_trln.tr.nim.ckt \
    d_vccap.1.ckt \
    d_vccap.2.ckt \
    d_vcg.1.ckt \
    d_vcg.2.ckt \
    d_vcg.3.ckt \
    d_vcg.4.ckt \
    d_vcg.5.ckt \
    d_vcr.1.ckt \
    d_vcr.2.ckt \
    d_vcr.3.ckt \
    d_vcr.4.ckt \
    d_vcr.5.ckt \
    d_vcvs.1.ckt \
    e_ccsrc.1.ckt \
    e_ccsrc.2.ckt \
    e_node_probes.ckt \
    eq3-1153.ckt \
    eq3-1153.fg2.ckt \
    eq3-1153.fg.ckt \
    eq3-1153.tr-b.ckt \
    eq3-1153.tr-nb.ckt \
    eq3-1153.tr-ni.ckt \
    eq4-9217.tran.ckt \
    eq4-9217.tran.nopreq.ckt \
    eqb1.ckt \
    eqboost.ckt \
    eqflat.ckt \
    eqmodify.ckt \
    e_storag.1.ckt \
    l2.ckt \
    l.ckt \
    list.1.ckt \
    ll1.ckt \
    ll.ckt \
    lll.ckt \
    llll.ckt \
    m_expression.1.ckt \
    m_expression.2.ckt \
    m_expression.3.ckt \
    m_expression.error.1.ckt \
    named_nodes.ckt \
    nmos100.1.ckt \
    nmos100.ckt \
    nmos15a.ckt \
    nmos15.ckt \
    nmp100.ckt \
    nothing.ckt \
    opamp-ol.ckt \
    opamp-ol-disto.2.ckt \
    opamp-ol-disto.ckt \
    opamp-vf.1.ckt \
    opamp-vf.1d.ckt \
    opamp-vf.2.ckt \
    opamp-vf.3.ckt \
    opamp-vf.ckt \
    oscillator.1.ckt \
    oscillator.2.ckt \
    oscillator.3.ckt \
    oscillator.4.ckt \
    oscillator.5.ckt \
    oscillator.6.ckt \
    oscillator.7.ckt \
    param.0a.ckt \
    param.0.ckt \
    param.10.ckt \
    param.11.ckt \
    param.12.ckt \
    param.1.ckt \
    param.2a.ckt \
    param.2b.ckt \
    param.2c.ckt \
    param.2.ckt \
    param.3.ckt \
    param.4.ckt \
    param.5.ckt \
    param.6.ckt \
    param.7.ckt \
    param.8.ckt \
    param.9.ckt \
    parse.1.ckt \
    sl.ckt \
    s_tr.1.ckt \
    trcurve2.ckt \
    trcurve.ckt \
    trcurves.ckt

# bm_cond stuff
XFAIL_TESTS = list.1.ckt

# recheck
XFAIL_TESTS += lang_verilog.3.gc

# typename issue, some other bugs
XFAIL_TESTS += lang_spectre.2.gc

# something wrong with .op command
XFAIL_TESTS += lang_spice.1b.gc

# needs work
XFAIL_TESTS += d_bjt-diffpair-tr2.ckt

# interactive mode broken?
XFAIL_TESTS += c_attach.1.gc

# cond parser broken
XFAIL_TESTS+= param.2c.ckt

# logic is broken. needs work
LOGIC_TESTS = \
    d_logic_nand-dc.ckt \
    d_logic_tr_ab1.ckt \
    d_logic_tr_ab2.ckt \
    d_logic_tr_an1.ckt \
    d_logic_tr_an2.ckt \
    d_logic_tr_db1.ckt \
    d_logic_tr_db2.ckt \
    d_logic_tr_dn1.ckt \
    d_logic_tr_dn2.ckt \
    d_logic_tr_mb1a.ckt \
    d_logic_tr_mb1b.ckt \
    d_logic_tr_mb2a.ckt \
    d_logic_tr_mb2b.ckt \
    d_logic_tr_mn1a.ckt \
    d_logic_tr_mn1b.ckt \
    d_logic_tr_mn2a.ckt \
    d_logic_tr_mn2b.ckt

# these are really broken...
XFAIL_TESTS += \
    c_delete.6.ckt \
    c_delete.7.ckt \
    m_expression.3.ckt \
    param.9.ckt \
    param.8.ckt

# functions need to return expressions, not doubles
XFAIL_TESTS+= c_measure.1.ckt

# something wrong with probes
XFAIL_TESTS+= print_v.sh

# these tests require gsl
GSL_TESTS = \
    d_ttcap.1.gc \
    pfet_shutdown.gc \
    s_ttt.rc.gc \
    s_ttt.lc.gc \
    s_ttt.lrc.gc

GSL_TESTS_LOG = $(GSL_TESTS:%.gc=%.log)

if USE_GSL
# don't skip GSL tests
else
$(GSL_TESTS_LOG): GC_LOG_FLAGS=SKIP
endif

# suspended.
d_bti.2.log: GC_LOG_FLAGS += SKIP

GC_TESTS = \
    $(GSL_TESTS) \
    ac_lp.gc \
    ac_noise.gc \
    alter.gc \
    bm_pulse.5.gc \
    bm_pulse.6.gc \
    bm_pulse.7.gc \
    bm_pulse.8.gc \
    bm_pulse.9.gc \
    bm_pulse.10.gc \
    bm_pulse.11.gc \
    bm_pulse.periodic.gc \
    c_attach.1.gc \
    c_attach.2.gc \
    c_debug.1.gc \
    c_debug.2.gc \
    c_echo.1.gc \
    c_for.1.gc \
    c_for.2.gc \
    c_for.3.gc \
    c_save.gc \
    c_printf.1.gc \
    c_printf.2.gc \
    c_wave.1.gc \
    c_wave.2.gc \
    c_wave.3.gc \
    c_wave.4.gc \
    c_wave.5.gc \
    c_wave.6.gc \
    c_wave.build.gc \
    c_wave.dump.gc \
    d_branchvs.1.gc \
    d_branchvs.2.gc \
    d_coil.brl.1.gc \
    d_coil.brl.2.gc \
    d_coil.brl.3.gc \
    d_cpoly_g.1.gc \
    d_cpoly_g.2.gc \
    d_cpoly_g.3.gc \
    d_cpoly_g.4.gc \
    d_mos.l0.gc \
    d_mos8.np.gc \
    d_res.gc \
    d_subckt.parm.gc \
    d_switch.disc.gc \
    d_switch.ic.gc \
    d_switch.tt.gc \
    dc_loop.gc \
    dc_itl.gc \
    dc_reverse.gc \
    dc_sweep.gc \
    dc_sweeps.gc \
    dc_trace.gc \
    dc_twice.gc \
    dc_twice.nopreq.gc \
    ddc_brl.gc \
    ddc_bug.gc \
    ddc_cc.gc \
    ddc_coil.gc \
    ddc_osc.gc \
    ddc_rc.gc \
    ddc_simple.gc \
    ddc_param.gc \
    ddc_sweep_twice.gc \
    ddc_v.1.gc \
    d_vs.1.gc \
    eval_param.gc \
    gnd.gc \
    hangbug.gc \
    lang_verilog.1a.gc \
    lang_verilog.1.gc \
    lang_verilog.2a.gc \
    lang_verilog.2b.gc \
    lang_verilog.2.gc \
    lang_verilog.3.gc \
    lang_verilog.comm.gc \
    lang_spectre.1.gc \
    lang_spectre.1a.gc \
    lang_spectre.2.gc \
    lang_spice.1.gc \
    lang_spice.1a.gc \
    lang_spice.1b.gc \
    lang_spice.comm.gc \
    LogicAND.gc \
    lring.gc \
    mark_dc_cont.gc \
    models-uf.gc \
    mosbug.2.gc \
    nodeorder.gc \
    node.gc \
    nodis.gc \
    param.13.gc \
    param_simple.gc \
    paramtest.gc \
    poly1_parse.gc \
    pulse.gc \
    ring.gc \
    s_sock.1.gc \
    s_sock.1b.gc \
    s_sock.1.preq.gc \
    s_sock.1b.preq.gc \
    s_sock.2.gc \
    s_sock.3.gc \
    s_sock.4.gc \
    s_sock.load.gc \
    s_tr.2.gc \
    s_tr.3.gc \
    s_tr.4.gc \
    s_tr.5.gc \
    s_tr.6.gc \
    s_tr.7.gc \
    s_tr.cont.gc \
    s_tr.lrc.gc \
    s_ttt.1.gc \
    s_ttt.label.gc \
    s_ttt.vib.gc \
    sens_rr.gc \
    spice_comment_continuation.gc \
    spice_mos.gc \
    temp.1.gc \
    tran_temp.gc \
    tr_bug.gc \
    tr_cont.gc \
    tr_cont.1.gc \
    tr_cont.2.gc \
    tr_cont.3.gc \
    tr_cont.4.gc \
    tr_cont_dc.gc \
    tr_mark.1.gc \
    tr_mark.2.gc \
    tr_mark.3.gc \
    tr_mark.gc \
    tr_nothing.gc \
    tt_twice.gc \
    ttr_cont.gc \
    tw_steps.gc \
    uic_cap.gc \
    v_uic.gc \
    vrc_cont.1.gc \
    vrc_evt.gc \
    Xnodeorder.gc \
    Xprobes.gc

#this is not ready yet.
XFAIL_TESTS+= s_tr.cont.gc
XFAIL_TESTS+= s_tr.lrc.gc

#reltol problem.
XFAIL_TESTS+= d_coil.brl.3.gc

# not implemented
XFAIL_TESTS+= c_for.3.gc

# abother bug
XFAIL_TESTS+= ddc_v.1.gc

# hmmm
XFAIL_TESTS+= ac_noise.gc

# somethings wrong
XFAIL_TESTS+= d_mos.l0.gc

# does not converge (why?)
XFAIL_TESTS+= mosbug.2.gc

CUSTOM_TESTS = c_attach.1.sh nosuchfile.sh

AUTO_TESTS = $(CKT_TESTS) $(CUSTOM_TESTS)

gsl_tests: $(GSL_TESTS:%.gc=%.log)
gc_tests: $(GC_TESTS:%.gc=%.log)
uf_tests: $(UF_TESTS:%.gc=%.log)
ckt_tests: $(CKT_TESTS:%.ckt=%.log)
logic_tests: $(LOGIC_TESTS:%.ckt=%.log)

.PHONY: .P

.P:
	@:

# $(AUTO_TESTS:%.sh=%.out): .P $(FILTER)
$(GC_TESTS:%.gc=%.gc.out): .P $(FILTER)
$(GC_TESTS:%.gc=%.log): .P $(FILTER)
$(CKT_TESTS:%.ckt=%.log): .P $(FILTER)
$(CKT_TESTS:%.ckt=%.ckt.out): .P $(FILTER)
# $(LOGIC_TESTS:%.ckt=%.log): .P $(FILTER)
$(LOGIC_TESTS:%=%.out): .P $(FILTER)
$(UF_TESTS:%.gc=%.gc.out): .P $(FILTER)
$(UF_TESTS:%.gc=%.log): .P $(FILTER)

# $(CKT_TESTS): %.sh: %.out

# try this first
$(CKT_TESTS:%.ckt:%.log): gnucap_run.log

# this will trigger a warning in autogen.sh
# (apparently this warning is s bug)
c_getckt.ckt.out: | $(PWD)/d_cccs.1.ckt
c_lib.1.ckt.out: | $(PWD)/c_lib.1.ckt

lang_spice.1a.log lang_verilog.1a.log: soclean | $(PWD)/bm_wrapper.cc
lang_spice.1a.gc.out lang_verilog.1a.gc.out: POSTDIFF+= '/^CXX/d'

$(GSL_TESTS:%=%.out) $(UF_TESTS:%=%.out): POSTDIFF+='/already·installed/{N;d}'

UF_TESTS = \
    d_bti.1.gc \
    d_bti.2.gc \
    d_bti.3.gc \
    d_bti.steps.gc \
    d_bti.symm.gc \
    d_hci.1.gc \
    d_hci.2.gc \
    d_hci.3.gc \
    d_hci.4.gc \
    d_hci.5.gc \
    d_hci.6.gc \
    d_hci.7.gc \
    d_hci_pd.gc \
    d_rcd.1.gc \
    d_rcd.2.gc \
    d_rcd.3.gc \
    d_rcd.neg.gc \
    hci_dumpload.gc \
    hci_mos.gc \
    hci_stepping.gc \
    no_hci.gc \
    param.14.gc \
    param_scopes.gc \
    pfet_age.gc \
    pfet_btisum.gc \
    pfet_age_pd.gc \
    pfet_age_pd.1.gc \
    pfet_age_tr.gc \
    pfet_rcd_param.gc \
    rcd_param.gc \
    rcd_tr.gc \
    rcd_tr.1.gc \
    rcd_tr.2.gc \
    rcd_ueff.gc \
    rcd_tw.gc \
    rcd_tw_dumpload.gc \
    sens_inv.gc \
    sens_mos8.gc \
    sens_uic.gc \
    spectre_param.gc

soclean:
	rm -f bm_wrapper.so

PREDIFF= \
   '1,/^core-lib·version/d' \
   '/^@/d'

## newer gnucap
#PREDIFF+='1,/^default·plugins/d'

# permit different control values
d_logic_tr_%.out: POSTDIFF+= '3·s/[01]\.········$$/X.········/'

# just zero noise
trcurve.ckt.out: POSTDIFF+= '54·s/-196\.E-24/·0.······/g' \
                        '411·s/727\.2E-18/0.·······/g' \
                        '411·s/317\.6E-18/0.·······/g' \
                        '462·s/795\.6E-18/0.·······/g' \
                        '462·s/378\.2E-18/0.·······/g' \
                        '513·s/864\.E-18/0.······/g' \
                        '513·s/438\.4E-18/0.·······/g'

$(CKT_TESTS:%=%.out) $(LOGIC_TESTS:%=%.out): %.ckt.out: %.ckt
	@LD_LIBRARY_PATH='../plugins/.libs$(LD_LIBRARY_PATH:%=:$(LD_LIBRARY_PATH))' \
	    ../src/gnucap $(GNUCAP_PLUGINS:%=-a %) $(GNUCAP_ARGS) -b $< | \
	    sed $(subst ·, ,$(PREDIFF:%=-e %)) \
	    $(if $(POSTDIFF),| sed) $(subst ·, ,$(POSTDIFF:%=-e %)) > $*.ckt.out

# close to zero...
d_vcg.%.out: POSTDIFF+= '27·s/975\.[0-9][0-9]f/0.·····/'

# close enough
bm_poly.%.out: POSTDIFF+= '30·s/-.\.f/·0.·/g'
bm_poly.%.out: POSTDIFF+= '30·s/-..\.f/·0.··/g'

# output seems to depend on lapack/architecture.
ddc_cc.gc.out: POSTDIFF+= 's/NaN/0.·/g'

# same
trcurve2.ckt.out: POSTDIFF+= 's/·[0-9]\.E-21/·0.····/'
trcurve2.ckt.out: POSTDIFF+= 's/·[0-9][0-9]\.E-21/·0.·····/'
trcurve2.ckt.out: POSTDIFF+= '11·s/.....E-30/0.·······/g'
trcurve2.ckt.out: POSTDIFF+= '17·s/3\.628E-21/0.·······/g'
trcurve2.ckt.out: POSTDIFF+= '17·s/110\.4E-30/0.·······/g'
trcurve2.ckt.out: POSTDIFF+= '23·s/5\.094E-2./0.·······/g'
trcurve2.ckt.out: POSTDIFF+= '23·s/73\.45E-21/0.·······/g'
trcurve2.ckt.out: POSTDIFF+= '29·s/630\.3E-24/0.·······/'
trcurve2.ckt.out: POSTDIFF+= '35·s/7\.08E-21/0.······/'
trcurve2.ckt.out: POSTDIFF+= '35·s/13\.28E-21/0.·······/'
trcurve2.ckt.out: POSTDIFF+= 's/[0-9]....E-30/0.·······/g'

bm_pulse.5.gc.out: POSTDIFF+= 's/[0-9]\....E-18/0.·······/g'

# negative zeroes? nice, but should not be machine dependent
bm_sffm.%.out: POSTDIFF+= 's/-0\.000\>/·0.000/g'

bm_table.2.4.ckt.out: POSTDIFF+= '54·s/[- ][0-9]\.f/·0.·/g'

c_attach.1.out: POSTDIFF+= 's\#`.*/tests\#`.../tests\#g'

FILTER = ./gnucap-testfilter

define testgen
	@echo "#!/bin/sh" > $@
	@echo "set -e" >> $@
	@echo "exec 2>&9" >> $@
	@echo "$(FILTER) $*.out < $< | diff -up - $< --label $< || exit 1" >> $@
	@echo "rm $*.out" >> $@
	@chmod +x $@
endef

$(CKT_TESTS:%.ckt=%.sh): %.sh: %.ckt.ref
	$(AM_V_GEN)$(testgen)

$(UF_TESTS:%=%.out) $(GC_TESTS:%=%.out): %.gc.out: %.gc
	$(AM_V_GEN)GNUCAP_PREQ=1 $(GNUCAP) $(GNUCAP_PLUGINS:%=-a ../modules/.libs/%) $(GNUCAP_ARGS) < $< | \
	    sed $(subst ·, ,$(PREDIFF:%=-e %)) \
	    $(if $(POSTDIFF),| sed) $(subst ·, ,$(POSTDIFF:%=-e %)) > $*.gc.out

# $(UF_TESTS): GNUCAP_PLUGINS = 'public models-uf.so'

# nosuchfile.out c_attach.1.out: %.out: %.gc
# 	-$(GNUCAP) < $< | \
# 	    sed $(subst ·, ,$(PREDIFF:%=-e %)) \
# 	    $(if $(POSTDIFF),| sed) $(subst ·, ,$(POSTDIFF:%=-e %)) > $*.out

LT=<
nosuchfile.sh: LT=<

# make sure we don't chdir...
lang_verilog.1a.sh: LT=<
lang_spice.1.sh: LT=<
c_attach.1.sh: I=-i

nosuchfile.sh: %.sh: %.gc Makefile %.ref
	@echo "#!/bin/sh" > $@
	@echo "set -e" >> $@
	@echo "GNUCAP=\$${GNUCAP-gnucap}" >> $@
	@echo "eval \$$REDIRECT" >> $@
	@echo "\$$GNUCAP \$$GNUCAP_ARGS $(LT) $I $< | sed $(subst ·, ,$(PREDIFF:%=-e %)) \\" >> $@
	@echo "$(if $(POSTDIFF),| sed) $(subst ·, ,$(POSTDIFF:%=-e %)) > $*.out" >> $@
	@echo "$(FILTER) $*.out < $(lastword $+) | diff -up - $(lastword $+) --label $(lastword $+) || exit 1" >> $@
	@echo "rm $*.out" >> $@
	@chmod +x $@

CLEANFILES = *.out $(AUTO_TESTS) core \
             patchlev_fake_0.h \
             patchlev_fake_999.h \
             patchlev_other.h

# hmmm.
if VPATH_BUILD
$(PWD)/%.ckt: %.ckt
	$(LN_S) $< $@

$(PWD)/%.cc: %.cc
	$(LN_S) $< $@

CLEANFILES+= d_cccs.1.ckt c_lib.1.ckt bm_wrapper.so
endif

EXTRA_DIST = gnucap_run.sh \
             gc_log_compiler \
             ckt_log_compiler \
             $(CKT_TESTS) $(CKT_TESTS:%=%.ref) \
             $(LOGIC_TESTS) $(LOGIC_TESTS:%=%.ref) \
             $(GC_TESTS) $(GC_TESTS:%.gc=%.ref) \
             $(UF_TESTS) $(UF_TESTS:%.gc=%.ref) \
             c_attach.1.gc \
             bm_wrapper.cc \
             disp_test.cc \
             bm_generator.doc bm_pulse.doc bm_pwl.doc bm_sin.doc d_cap.doc d_diode.doc \
             d_.doc d_subckt.doc u_.doc \
             d_coil.5.spice \
             d_mos1.lin1.spice \
             d_mos1.sat2.spice \
             d_mos6.lin1r.spice \
             d_mos6.lin1.spice \
             d_mos6.lin1z.spice \
             d_mos6.lin2r.spice \
             d_mos6.lin2.spice \
             d_mos6.lin2z.spice \
             d_mos6.sat1r.spice \
             d_mos6.sat1.spice \
             d_mos6.sat1s.spice \
             d_mos6.sat1v.spice \
             d_mos6.sat2r.spice \
             d_mos6.sat2.spice \
             d_mos6.sat3r.spice \
             d_mos6.sat3.spice \
             pwl.spice \
				 ngspice \
             nobypass \
             nolubypass \
             vs.spice

.PHONY: check

AM_CPPFLAGS = -I$(top_builddir)/include -I$(top_srcdir)/include \
              -I$(top_builddir)/src -I$(top_srcdir)/src
AM_DEFAULT_SOURCE_EXT = .cc
AM_LDFLAGS = -module -avoid-version -rpath /dev/null
check_LTLIBRARIES = \
	disp_test.la \
	bogus.la bogus1.la bogus2.la bogus3.la bogus4.la

# $(CUSTOM_TESTS:%.sh=%.log): %.log: %.out
	
c_attach.1.out: empty.cc bogus.so d_vs.cc loadtestdir loadtestdirro attachclean
c_attach.2.out c_attach.2.log: bogus.so bogus1.so bogus2.so bogus3.so bogus4.so

c_debug.2.gc.out: disp_test.so

attachclean:
	rm -f d_vs.so load*/*.so empty.so

loadtestdir:
	mkdir $@
	cd $@; $(LN_S) $(top_srcdir)/src/d_cccs.cc .
	cd $@; $(LN_S) $(top_srcdir)/src/d_ccvs.cc .

loadtestdirro:
	mkdir $@
	cp $(top_srcdir)/src/d_ccvs.cc $@

empty.cc:
	touch $@

d_vs.cc: | $(top_srcdir)/src/d_vs.cc
	$(LN_S) $| .

check_SO = \
	disp_test.so \
	bogus.so bogus1.so bogus2.so bogus3.so bogus4.so

$(check_SO): %.so: %.la
	[ -f $@ ] || $(LN_S) .libs/$@ .

.PHONY: c_attach.1.out attachclean soclean

# bug: @@@ also catches incomplete
# bug: case labels don't end with {
sedscript:
	grep -e ^@@@: -e ^@i@ *.log | sort -u | tr : " " | \
		$(AWK) '{printf "sed -i %s -e '"'"'%s s/\\({\\|:\\)\\(\\ \\|\\t\\)\\?\\(un\\|i\\)tested();$$/\\1/'"'"' # %s +%s\n", $$3, $$4, $$3, $$4}' > $@

# checked in files. sort out later.
EXTRA_DIST += \
bm_pulse.1a.ckt \
bogus.cc \
bogus2.cc \
bti_dumpload.gc \
powerdown.gc \
print_v.ckt \
pulse_freq.gc \
pulse_freq.ref \
rcd_cons_tw.sim \
rcd_cons_tw_more.sim \
rcd_tr.sim \
rcd_tw.sim \
res_cond.gc \
saveload.sim \
sweep.sim \
tt_inv.gc \
tw_trace.gc
