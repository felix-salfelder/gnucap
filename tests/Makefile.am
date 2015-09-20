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
TEST_EXTENSIONS = .sh

TESTS = \
    gnucap_run.sh \
    $(GC_TESTS) \
    $(CKT_TESTS)

REF = ==out

transform=s&\.real$$&&;$(program_transform_name)
plugindir_transformed = ${libdir}/$(shell echo gnucap | sed -e "$(transform)")
plugindir:= $(plugindir_transformed)
plugin_PROGRAMS = gnucap-testfilter
gnucap_testfilter_SOURCES = filter.c
empty=
$(empty)install-pluginPROGRAMS: transform=

CKT_TESTS = \
    1lin.sh \
    ap_ctof.2.sh \
    ap_ctof.sh \
    basic.sh \
    bm_complex.1a.sh \
    bm_complex.1.sh \
    bm_complex.2.sh \
    bm_complex.3.sh \
    bm_complex.4.sh \
    bm_cond.1.sh \
    bm_cond.2.sh \
    bm_cond.3.sh \
    bm_exp.1.sh \
    bm_exp.1p.sh \
    bm_exp.2.sh \
    bm_exp.3.sh \
    bm_exp.3p.sh \
    bm_exp.4p.sh \
    bm_exp.5p.sh \
    bm_fit.0.1.sh \
    bm_fit.1.1.sh \
    bm_fit.2.1.sh \
    bm_fit.2.2.sh \
    bm_fit.2.3.sh \
    bm_fit.2.4.sh \
    bm_fit.3.1.sh \
    bm_fit.3.2.sh \
    bm_fit.3.3.sh \
    bmm_cap.1.sh \
    bmm_cap.2.sh \
    bmm_cap.3.sh \
    bmm_cap.4.sh \
    bmm_cap.5.sh \
    bmm_res.1.sh \
    bmm_res.1p.sh \
    bmm_res.2.sh \
    bmm_res.2p.sh \
    bmm_res.3.sh \
    bmm_res.3p.sh \
    bmm_res.4.sh \
    bmm_res.4p.sh \
    bmm_res.5.sh \
    bmm_res.5p.sh \
    bm_poly.1.sh \
    bm_poly.2.sh \
    bm_posy.1.sh \
    bm_pulse.1.sh \
    bm_pulse.2.sh \
    bm_pulse.3.sh \
    bm_pwl.1.sh \
    bm_pwl.2.sh \
    bm_pwl.3.sh \
    bm_pwl.4.sh \
    bm_pwl.5.sh \
    bm_sffm.1.sh \
    bm_sffm.2.sh \
    bm_sffm.3.sh \
    bm_sffm.4.sh \
    bm_sffm.5.sh \
    bm_sffm.6.sh \
    bm_sin.1.sh \
    bm_sin.2.sh \
    bm_sin.3.sh \
    bm_sin.4.sh \
    bm_sin.5.sh \
    bm_sin.6.sh \
    bm_sin.7.sh \
    bm_sin.8.sh \
    bm_sin.9.sh \
    bm_table.1.1.sh \
    bm_table.2.1.sh \
    bm_table.2.2.sh \
    bm_table.2.3.sh \
    bm_table.2.4.sh \
    bm_table.3.1.sh \
    bm_table.3.2.sh \
    bm_table.3.3.sh \
    bm_table.error1.sh \
    bm_tanh.1.sh \
    bm_tanh.2.sh \
    cc.sh \
    c.sh \
    c_delete.2.sh \
    c_delete.3.sh \
    c_delete.4.sh \
    c_delete.5.sh \
    c_delete.6.sh \
    c_delete.7.sh \
    c_delete.sh \
    c_eval.1.sh \
    c_getckt.sh \
    charge-cons2.sh \
    charge-cons3.sh \
    charge-cons.sh \
    c_help.1.sh \
    c_lib.1.sh \
    c_list.sh \
    c_list.1.sh \
    c_measure.1.sh \
    c_measure.2.sh \
    c_modify.sh \
    cont.1.sh \
    cont.2.sh \
    cont.3.sh \
    current.sh \
    d_bjt.2.sh \
    d_bjt.sh \
    d_bjt-diffpair-ccs.sh \
    d_bjt-diffpair-cjc.sh \
    d_bjt-diffpair-dc.sh \
    d_bjt-diffpair-nocap.sh \
    d_bjt-diffpair-tf.sh \
    d_bjt-diffpair-tran.sh \
    d_bjt-diffpair-tr.sh \
    d_bjt-diffpair-tr2.sh \
    d_bjt.error1.sh \
    d_bjt-pnp.sh \
    d_bjt-rca3040.sh \
    d_bjt-rtlinv.sh \
    d_bjt-schmitt-bypass.sh \
    d_bjt-schmitt-nobypass.sh \
    d_cap.1.sh \
    d_cap.2.sh \
    d_cap.3a.sh \
    d_cap.3.sh \
    d_cap.4.sh \
    d_cap.5.sh \
    d_cap.6.sh \
    d_cap.ic1.sh \
    d_cap.ic2.sh \
    d_cap.ic4.sh \
    d_cccs.1.sh \
    d_cccs.2.sh \
    d_cccs.3.sh \
    d_cccs.4.sh \
    d_cccs.5.sh \
    d_cccs.6.sh \
    d_ccvs.1.sh \
    d_ccvs.2.sh \
    d_ccvs.3.sh \
    d_coil.1a.sh \
    d_coil.1b.sh \
    d_coil.1c.sh \
    d_coil.1.sh \
    d_coil.1d.sh \
    d_coil.1e.sh \
    d_coil.2.sh \
    d_coil.3.sh \
    d_coil.4.sh \
    d_coil.5.sh \
    d_coil.error1.sh \
    d_coil.error2.sh \
    d_coil.error3.sh \
    d_coil.error4.sh \
    d_coil.error5.sh \
    d_coil.m3a.sh \
    d_coil.m3.sh \
    d_diode.10a.sh \
    d_diode.10.sh \
    d_diode.11a.sh \
    d_diode.11b.sh \
    d_diode.11c.sh \
    d_diode.11.sh \
    d_diode.11d.sh \
    d_diode.12.sh \
    d_diode.1.sh \
    d_diode.2.sh \
    d_diode.3.sh \
    d_diode.4.sh \
    d_diode.5.sh \
    d_diode.6.sh \
    d_diode.7.sh \
    d_diode.8.sh \
    d_diode.9.sh \
    d_diode.error1.sh \
    d_logic.1.sh \
    d_logic.1a.sh \
    d_logic.2.sh \
    d_logic.3.sh \
    d_logic.4.sh \
    d_logic.5.sh \
    d_logic.6.sh \
    d_logic.7.sh \
    d_logic.error1.sh \
    d_macro.1.sh \
    d_mos1.bin1.sh \
    d_mos1.bin2.sh \
    d_mos1.bin3.sh \
    d_mos1.inv1.sh \
    d_mos1.inv2.sh \
    d_mos1.inv3.sh \
    d_mos1.lin1.sh \
    d_mos1.lin1.late.sh \
    d_mos1.lin1.none.sh \
    d_mos1.lin2.sh \
    d_mos1.lin3.sh \
    d_mos1.lin4.sh \
    d_mos1.lin5.sh \
    d_mos1.n0.sh \
    d_mos1.n-1.sh \
    d_mos1.n+1.sh \
    d_mos1.nand1.sh \
    d_mos1.nand2.sh \
    d_mos1.p0.sh \
    d_mos1.p-1.sh \
    d_mos1.p+1.sh \
    d_mos1.plin.sh \
    d_mos1.psat1.sh \
    d_mos1.psat2.sh \
    d_mos1.psat3.sh \
    d_mos1.sat1.sh \
    d_mos1.sat2.sh \
    d_mos1.sat3.sh \
    d_mos1.sat4.sh \
    d_mos1.sat5.sh \
    d_mos1.sat6.sh \
    d_mos1.sat7.sh \
    d_mos1.sts.sh \
    d_mos1.zero.sh \
    d_mos2.inv1.sh \
    d_mos2.lin1.sh \
    d_mos2.lin1.late.sh \
    d_mos2.lin1r.sh \
    d_mos2.lin1z.sh \
    d_mos2.lin2.sh \
    d_mos2.lin3.sh \
    d_mos2.lin4.sh \
    d_mos2.lin5.sh \
    d_mos2.n0.sh \
    d_mos2.n-1.sh \
    d_mos2.n+1.sh \
    d_mos2.nand1.sh \
    d_mos2.p0.sh \
    d_mos2.p-1.sh \
    d_mos2.p+1.sh \
    d_mos2.plin.sh \
    d_mos2.psat1.sh \
    d_mos2.psat2.sh \
    d_mos2.psat3.sh \
    d_mos2.sat1.sh \
    d_mos2.sat2.sh \
    d_mos2.sat3.sh \
    d_mos2.sat4.sh \
    d_mos2.sat5.sh \
    d_mos2.sat6.sh \
    d_mos2.sat7.sh \
    d_mos2.sts.sh \
    d_mos2.zero.sh \
    d_mos3.inv1.sh \
    d_mos3.lin1.sh \
    d_mos3.lin1r.sh \
    d_mos3.lin1z.sh \
    d_mos3.lin2.sh \
    d_mos3.lin3.sh \
    d_mos3.lin4.sh \
    d_mos3.lin5.sh \
    d_mos3.n0.sh \
    d_mos3.n-1.sh \
    d_mos3.n+1.sh \
    d_mos3.nand1.sh \
    d_mos3.p0.sh \
    d_mos3.p-1.sh \
    d_mos3.p+1.sh \
    d_mos3.plin.sh \
    d_mos3.psat1.sh \
    d_mos3.psat2.sh \
    d_mos3.psat3.sh \
    d_mos3.sat000.sh \
    d_mos3.sat007.sh \
    d_mos3.sat1.sh \
    d_mos3.sat2.sh \
    d_mos3.sat3.sh \
    d_mos3.sat4.sh \
    d_mos3.sat5.sh \
    d_mos3.sat6.sh \
    d_mos3.sat7.sh \
    d_mos3.sts000.sh \
    d_mos3.sts.sh \
    d_mos3.zero.sh \
    d_mos4.1.sh \
    d_mos4.2.sh \
    d_mos4.3.sh \
    d_mos49.nand1.sh \
    d_mos5.1.sh \
    d_mos6.inv1.sh \
    d_mos6.lin1.sh \
    d_mos6.lin1r.sh \
    d_mos6.lin1z.sh \
    d_mos6.lin2.sh \
    d_mos6.lin2r.sh \
    d_mos6.lin2z.sh \
    d_mos6.n0.sh \
    d_mos6.n-1.sh \
    d_mos6.n+1.sh \
    d_mos6.nand1.sh \
    d_mos6.p0.sh \
    d_mos6.p-1.sh \
    d_mos6.p+1.sh \
    d_mos6.plin.sh \
    d_mos6.psat1.sh \
    d_mos6.psat2.sh \
    d_mos6.psat3.sh \
    d_mos6.sat1.sh \
    d_mos6.sat1r.sh \
    d_mos6.sat1s.sh \
    d_mos6.sat1v.sh \
    d_mos6.sat2.sh \
    d_mos6.sat2r.sh \
    d_mos6.sat3.sh \
    d_mos6.sat3r.sh \
    d_mos6.sts000.sh \
    d_mos6.sts.sh \
    d_mos6.zero.sh \
    d_mos7.alpha.sh \
    d_mos7.inv1.sh \
    d_mos7.lin1.sh \
    d_mos7.lin2.sh \
    d_mos7.nand1.sh \
    d_mos7.nand1.nobypass.sh \
    d_mos7.plin.sh \
    d_mos7.psat.sh \
    d_mos7.sat1.sh \
    d_mos7.sat2.sh \
    d_mos7.sat3.sh \
    d_mos7.sts.sh \
    d_mos7.zero.sh \
    d_mos8c2.nand1.sh \
    d_mos8c2.nand1.nobypass.sh \
    d_mos8.nand1.sh \
    d_mos8.sat2.sh \
    d_mos8.sat3.sh \
    d_mos8.sat4.sh \
    d_mos.error1.sh \
    d_mos.error2.sh \
    d_mos.error3.sh \
    d_mos.error4.sh \
    d_subckt.1.sh \
    d_subckt.2.sh \
    d_subckt.3.sh \
    d_subckt.4.sh \
    d_subckt.5.sh \
    d_subckt.6.sh \
    d_subckt.7.sh \
    d_subckt.8.sh \
    d_subckt.9.sh \
    d_subckt.error1.sh \
    d_subckt.error2.sh \
    d_subckt.error3.sh \
    d_subckt.error4.sh \
    d_switch.c0.sh \
    d_switch.c1a.sh \
    d_switch.c1.sh \
    d_switch.c1.inc.sh \
    d_switch.c1.noinc.sh \
    d_switch.c1r.sh \
    d_switch.c2a.sh \
    d_switch.c2.sh \
    d_switch.c2r.sh \
    d_switch.c3.sh \
    d_switch.error1.sh \
    d_switch.error2.sh \
    d_switch.error3.sh \
    d_switch.nro.1.sh \
    d_switch.nro.1e.sh \
    d_switch.nro.1g.sh \
    d_switch.nro.1gr.sh \
    d_switch.nro.2.sh \
    d_switch.nro.3.sh \
    d_switch.nro.3h.sh \
    d_switch.nro.4.sh \
    d_switch.v.cont.sh \
    d_switch.v.sh \
    d_switch.vr.sh \
    d_tcap.1.sh \
    d_tcap.2.sh \
    d_tcap.3.sh \
    d_tcap.4.sh \
    d_trln.ac.sh \
    d_trln.tr.100.im.sh \
    d_trln.tr.100.nim.sh \
    d_trln.tr.300.im.sh \
    d_trln.tr.300.nim.sh \
    d_trln.tr.im.sh \
    d_trln.tr.m.sh \
    d_trln.tr.nim.sh \
    d_vccap.1.sh \
    d_vccap.2.sh \
    d_vcg.1.sh \
    d_vcg.2.sh \
    d_vcg.3.sh \
    d_vcg.4.sh \
    d_vcg.5.sh \
    d_vcr.1.sh \
    d_vcr.2.sh \
    d_vcr.3.sh \
    d_vcr.4.sh \
    d_vcr.5.sh \
    d_vcvs.1.sh \
    dc_temp.sh \
    dcsweep.sh \
    e_ccsrc.1.sh \
    e_ccsrc.2.sh \
    e_node_probes.sh \
    eq3-1153.sh \
    eq3-1153.fg2.sh \
    eq3-1153.fg.sh \
    eq3-1153.tr-b.sh \
    eq3-1153.tr-nb.sh \
    eq3-1153.tr-ni.sh \
    eq4-9217.tran.sh \
    eqb1.sh \
    eqboost.sh \
    eqflat.sh \
    eqmodify.sh \
    e_storag.1.sh \
    list.1.sh \
    l2.sh \
    l.sh \
    ll1.sh \
    ll.sh \
    lll.sh \
    llll.sh \
    m_expression.1.sh \
    m_expression.2.sh \
    m_expression.3.sh \
    m_expression.error.1.sh \
    named_nodes.sh \
    nmos100.1.sh \
    nmos100.sh \
    nmos15a.sh \
    nmos15.sh \
    nmp100.sh \
    nothing.sh \
    opamp-ol.sh \
    opamp-ol-disto.2.sh \
    opamp-ol-disto.sh \
    opamp-vf.1.sh \
    opamp-vf.1d.sh \
    opamp-vf.2.sh \
    opamp-vf.3.sh \
    opamp-vf.sh \
    oscillator.1.sh \
    oscillator.2.sh \
    oscillator.3.sh \
    oscillator.4.sh \
    oscillator.5.sh \
    oscillator.6.sh \
    oscillator.7.sh \
    param.0a.sh \
    param.0.sh \
    param.10.sh \
    param.11.sh \
    param.12.sh \
    param.1.sh \
    param.2a.sh \
    param.2b.sh \
    param.2c.sh \
    param.2.sh \
    param.3.sh \
    param.4.sh \
    param.5.sh \
    param.6.sh \
    param.7.sh \
    param.8.sh \
    param.9.sh \
    parse.1.sh \
    sl.sh \
    trcurve2.sh \
    trcurve.sh \
    trcurves.sh

XFAIL_TESTS = list.1.sh

# recheck
XFAIL_TESTS += lang_verilog.3.sh

# typename issue, some other bugs
XFAIL_TESTS += lang_spectre.2.sh

# something wrong with .op command
XFAIL_TESTS += lang_spice.1b.sh

# needs work
XFAIL_TESTS += d_bjt-diffpair-tr2.sh

# interactive mode broken?
XFAIL_TESTS += c_attach.1.sh

# cond parser broken
XFAIL_TESTS+= param.2c.sh

# logic is broken. needs work
LOGIC_TESTS = \
    d_logic_nand-dc.sh \
    d_logic_tr_ab1.sh \
    d_logic_tr_ab2.sh \
    d_logic_tr_an1.sh \
    d_logic_tr_an2.sh \
    d_logic_tr_db1.sh \
    d_logic_tr_db2.sh \
    d_logic_tr_dn1.sh \
    d_logic_tr_dn2.sh \
    d_logic_tr_mb1a.sh \
    d_logic_tr_mb1b.sh \
    d_logic_tr_mb2a.sh \
    d_logic_tr_mb2b.sh \
    d_logic_tr_mn1a.sh \
    d_logic_tr_mn1b.sh \
    d_logic_tr_mn2a.sh \
    d_logic_tr_mn2b.sh

# these are really broken...
XFAIL_TESTS += \
    c_delete.6.sh \
    c_delete.7.sh \
    m_expression.3.sh \
    param.9.sh \
    param.8.sh

# functions need to return expressions, not doubles
XFAIL_TESTS+= c_measure.1.sh

# something wrong with probes
XFAIL_TESTS+= print_v.sh

if USE_GSL
TESTS_GSL = \
    d_ttcap.1.sh \
    pfet_shutdown.sh
endif

SUSPENDED_TESTS = \
    d_bti.2.sh

GC_TESTS = \
    $(TESTS_GSL) \
    ac_lp.sh \
    ac_noise.sh \
    alter.sh \
    bm_pulse.4.sh \
    bm_pulse.5.sh \
    bm_pulse.6.sh \
    bm_pulse.7.sh \
    bm_pulse.8.sh \
    bm_pulse.9.sh \
    bm_pulse.10.sh \
    bm_pulse.11.sh \
    bm_pulse.periodic.sh \
    c_attach.1.sh \
    c_echo.1.sh \
    c_for.1.sh \
    c_for.2.sh \
    c_for.3.sh \
    c_save.sh \
    c_printf.1.sh \
    c_printf.2.sh \
    c_wave.1.sh \
    c_wave.2.sh \
    c_wave.3.sh \
    c_wave.4.sh \
    c_wave.5.sh \
    c_wave.6.sh \
    c_wave.build.sh \
    c_wave.dump.sh \
    d_branchvs.1.sh \
    d_branchvs.2.sh \
    d_bti.1.sh \
    d_bti.3.sh \
    d_bti.steps.sh \
    d_bti.symm.sh \
    d_coil.brl.1.sh \
    d_coil.brl.2.sh \
    d_coil.brl.3.sh \
    d_hci.1.sh \
    d_hci.2.sh \
    d_hci.3.sh \
    d_hci.4.sh \
    d_hci.5.sh \
    d_hci.6.sh \
    d_hci.7.sh \
    d_hci_pd.sh \
    d_mos8.np.sh \
    d_res.sh \
    d_switch.disc.sh \
    d_switch.ic.sh \
    d_switch.tt.sh \
    d_rcd.1.sh \
    d_rcd.2.sh \
    d_rcd.3.sh \
    d_rcd.neg.sh \
    dc_reverse.sh \
    dc_sweep.sh \
    dc_trace.sh \
    dc_twice.sh \
    ddc_brl.sh \
    ddc_bug.sh \
    ddc_cc.sh \
    ddc_coil.sh \
    ddc_osc.sh \
    ddc_rc.sh \
    ddc_simple.sh \
    ddc_param.sh \
    ddc_sweep_twice.sh \
    ddc_v.1.sh \
    d_vs.1.sh \
    eval_param.sh \
    gnd.sh \
    hangbug.sh \
    hci_dumpload.sh \
    hci_mos.sh \
    hci_stepping.sh \
    lang_verilog.1a.sh \
    lang_verilog.1.sh \
    lang_verilog.2a.sh \
    lang_verilog.2b.sh \
    lang_verilog.2.sh \
    lang_verilog.3.sh \
    lang_verilog.comm.sh \
    lang_spectre.1.sh \
    lang_spectre.1a.sh \
    lang_spectre.2.sh \
    lang_spice.1.sh \
    lang_spice.1a.sh \
    lang_spice.1b.sh \
    lang_spice.comm.sh \
    LogicAND.sh \
    lring.sh \
    mark_dc_cont.sh \
    mosbug.2.sh \
    no_hci.sh \
    nodeorder.sh \
    node.sh \
    nodis.sh \
    param.13.sh \
    param.14.sh \
    param_simple.sh \
    param_scopes.sh \
    paramtest.sh \
    pfet_age.sh \
    pfet_btisum.sh \
    pfet_age_pd.sh \
    pfet_age_pd.1.sh \
    pfet_age_tr.sh \
    pfet_rcd_param.sh \
    poly1_parse.sh \
    pulse.sh \
    rcd_param.sh \
    rcd_tr.sh \
    rcd_tr.1.sh \
    rcd_tr.2.sh \
    rcd_tw.sh \
    rcd_tw_dumpload.sh \
    rcd_ueff.sh \
    ring.sh \
    s_sock.1.sh \
    s_sock.1b.sh \
    s_sock.2.sh \
    s_tr.cont.sh \
    s_tr.lrc.sh \
    s_ttt.1.sh \
    s_ttt.label.sh \
    s_ttt.rc.sh \
    s_ttt.lc.sh \
    s_ttt.lrc.sh \
    s_ttt.vib.sh \
    sckt_bug.sh \
    sens_inv.sh \
    sens_mos8.sh \
    sens_rr.sh \
    sens_uic.sh \
    spectre_param.sh \
    spice_comment_continuation.sh \
    spice_mos.sh \
    tran_temp.sh \
    tr_bug.sh \
    tr_cont.sh \
    tr_cont.1.sh \
    tr_cont.2.sh \
    tr_cont.3.sh \
    tr_cont.4.sh \
    tr_cont_dc.sh \
    tr_mark.1.sh \
    tr_mark.2.sh \
    tr_mark.3.sh \
    tr_mark.sh \
    tr_nothing.sh \
    tt_twice.sh \
    ttr_cont.sh \
    tw_steps.sh \
    uic_cap.sh \
    v_uic.sh \
    vrc_cont.1.sh \
    vrc_evt.sh \
    Xnodeorder.sh \
    Xprobes.sh

#this is not ready yet.
XFAIL_TESTS+= s_tr.cont.sh
XFAIL_TESTS+= s_tr.lrc.sh

#reltol problem.
XFAIL_TESTS+= d_coil.brl.3.sh

# not implemented
XFAIL_TESTS+= c_for.3.sh

# abother bug
XFAIL_TESTS+= ddc_v.1.sh

# hmmm
XFAIL_TESTS+= ac_noise.sh

# somethings wrong
XFAIL_TESTS+= sckt_bug.sh

# does not converge (why?)
XFAIL_TESTS+= mosbug.2.sh

CUSTOM_TESTS = c_attach.1.sh nosuchfile.sh

AUTO_TESTS = $(CKT_TESTS) $(GC_TESTS) $(CUSTOM_TESTS)

gc_tests: $(GC_TESTS:%.sh=%.log)
ckt_tests: $(CKT_TESTS:%.sh=%.log)

.PHONY: .P

.P:
	@:

$(AUTO_TESTS:%.sh=%.out) $(GC_TESTS:%.sh=%.log): .P $(FILTER)

$(CKT_TESTS): %.sh: %.out

# try this first
$(CKT_TESTS:%.sh:%.log): gnucap_run.log

# this will trigger a warning in autogen.sh
# (apparently this warning is s bug)
c_getckt.out: | $(PWD)/d_cccs.1.ckt
c_lib.1.out: | $(PWD)/c_lib.1.ckt

lang_spice.1a.log lang_verilog.1a.log: soclean | $(PWD)/bm_wrapper.cc
lang_spice.1a.sh lang_verilog.1a.sh: POSTDIFF+= '/^CXX/d'

soclean:
	rm -f *.so

PREDIFF= \
   '1,/^core-lib·version/d' \
   '/^@/d'

## newer gnucap
#PREDIFF+='1,/^default·plugins/d'

# permit different control values
d_logic_tr_%.out: POSTDIFF+= '3·s/[01]\.········$$/X.········/'

# just zero noise
trcurve.out: POSTDIFF+= '54·s/-196\.E-24/·0.······/g' \
                        '411·s/727\.2E-18/0.·······/g' \
                        '411·s/317\.6E-18/0.·······/g' \
                        '462·s/795\.6E-18/0.·······/g' \
                        '462·s/378\.2E-18/0.·······/g' \
                        '513·s/864\.E-18/0.······/g' \
                        '513·s/438\.4E-18/0.·······/g'

$(CKT_TESTS:%.sh=%.out): %.out: %.ckt
	$(AM_V_at)LD_LIBRARY_PATH='../plugins/.libs$(LD_LIBRARY_PATH:%=:$(LD_LIBRARY_PATH))' \
	    ../src/gnucap $(GNUCAP_PLUGINS:%=-a %) $(GNUCAP_ARGS) -b $< | \
	    sed $(subst ·, ,$(PREDIFF:%=-e %)) \
	    $(if $(POSTDIFF),| sed) $(subst ·, ,$(POSTDIFF:%=-e %)) > $*.out

# close to zero...
d_vcg.%.out: POSTDIFF+= '27·s/975\.[0-9][0-9]f/0.·····/'

# close enough
bm_poly.%.out: POSTDIFF+= '30·s/-.\.f/·0.·/g'
bm_poly.%.out: POSTDIFF+= '30·s/-..\.f/·0.··/g'

# output seems to depend on lapack/architecture.
ddc_cc.sh: POSTDIFF+= 's/NaN/0.·/g'

# same
trcurve2.out: POSTDIFF+= 's/·[0-9]\.E-21/·0.····/'
trcurve2.out: POSTDIFF+= 's/·[0-9][0-9]\.E-21/·0.·····/'
trcurve2.out: POSTDIFF+= '11·s/.....E-30/0.·······/g'
trcurve2.out: POSTDIFF+= '17·s/3\.628E-21/0.·······/g'
trcurve2.out: POSTDIFF+= '17·s/110\.4E-30/0.·······/g'
trcurve2.out: POSTDIFF+= '23·s/5\.094E-2./0.·······/g'
trcurve2.out: POSTDIFF+= '23·s/73\.45E-21/0.·······/g'
trcurve2.out: POSTDIFF+= '29·s/630\.3E-24/0.·······/'
trcurve2.out: POSTDIFF+= '35·s/7\.08E-21/0.······/'
trcurve2.out: POSTDIFF+= '35·s/13\.28E-21/0.·······/'
trcurve2.out: POSTDIFF+= 's/[0-9]....E-30/0.·······/g'

# negative zeroes? nice, but should not be machine dependent
bm_sffm.%.out: POSTDIFF+= 's/-0\.000\>/·0.000/g'

bm_table.2.4.out: POSTDIFF+= '54·s/[- ][0-9]\.f/·0.·/g'

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

$(CKT_TESTS): %.sh: $(REF)/%.ckt.ref
	$(AM_V_GEN)$(testgen)

$(GC_TESTS:%.sh=%.out): %.out: %.gc
	$(AM_V_GEN)$(GNUCAP) $(GNUCAP_PLUGINS:%=-a ../modules/.libs/%) $(GNUCAP_ARGS) < $< | \
	    sed $(subst ·, ,$(PREDIFF:%=-e %)) \
	    $(if $(POSTDIFF),| sed) $(subst ·, ,$(POSTDIFF:%=-e %)) > $*.out

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

nosuchfile.sh $(GC_TESTS): %.sh: %.gc Makefile %.ref
	@echo "#!/bin/sh" > $@
	@echo "set -e" >> $@
	@echo "GNUCAP=\$${GNUCAP-gnucap}" >> $@
	@echo "eval \$$REDIRECT" >> $@
	@echo "\$$GNUCAP \$$GNUCAP_ARGS $(LT) $I $< | sed $(subst ·, ,$(PREDIFF:%=-e %)) \\" >> $@
	@echo "$(if $(POSTDIFF),| sed) $(subst ·, ,$(POSTDIFF:%=-e %)) > $*.out" >> $@
	@echo "$(FILTER) $*.out < $(lastword $+) | diff -up - $(lastword $+) --label $(lastword $+) || exit 1" >> $@
	@echo "rm $*.out" >> $@
	@chmod +x $@

CLEANFILES = *.out $(AUTO_TESTS) core

# hmmm.
if VPATH_BUILD
$(PWD)/%.ckt: %.ckt
	$(LN_S) $< $@

$(PWD)/%.cc: %.cc
	$(LN_S) $< $@

CLEANFILES+= d_cccs.1.ckt c_lib.1.ckt bm_wrapper.so
endif

EXTRA_DIST = gnucap_run.sh \
             $(CKT_TESTS:%.sh=%.ckt) $(CKT_TESTS:%.sh===out/%.ckt.ref) \
             $(LOGIC_TESTS:%.sh=%.ckt) $(LOGIC_TESTS:%.sh===out/%.ckt.ref) \
             $(GC_TESTS:%.sh=%.gc) $(GC_TESTS:%.sh=%.ref) \
             c_attach.1.gc \
             bm_wrapper.cc \
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

AM_CPPFLAGS = -I$(top_builddir)/include -I$(top_srcdir)/include
AM_DEFAULT_SOURCE_EXT = .cc
AM_LDFLAGS = -module -avoid-version -rpath /dev/null
check_LTLIBRARIES = bogus.la

# $(CUSTOM_TESTS:%.sh=%.log): %.log: %.out
	
c_attach.1.out: empty.cc bogus.so d_vs.cc loadtestdir loadtestdirro attachclean

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

bogus.so: %.so: %.la
	[ -f $@ ] || $(LN_S) .libs/$@ .

.PHONY: c_attach.1.out attachclean soclean

# bug: @@@ also catches incomplete
# bug: case labels don't end with {
sedscript: 2
	grep -e ^@@@ -e ^@i@ 2 | sort -u | tr : " " | \
		$(AWK) '{printf "sed -i %s -e '"'"'%s s/\\({\\|:\\)\\(\\ \\|\\t\\)\\?\\(un\\|i\\)tested();$$/\\1/'"'"' # %s +%s\n", $$2, $$3, $$2, $$3}' > $@

# checked in files. sort out later.
EXTRA_DIST += \
bm_pulse.1a.ckt \
bogus.cc \
bti_dumpload.gc \
dc_sweep.ckt \
dc_sweeps.gc \
dc_trace.ckt \
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
